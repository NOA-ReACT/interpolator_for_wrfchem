import shutil
from pathlib import Path
from typing import Optional

import click
import numpy as np
import xarray as xr

from interpolator_for_wrfchem import utils
from interpolator_for_wrfchem.global_models import GLOBAL_MODELS
from interpolator_for_wrfchem.interpolation import interpolate_to_wrf
from interpolator_for_wrfchem.species_map import SpeciesMap, convert_si
from interpolator_for_wrfchem.wrf import WRFBoundary, WRFInput


def _parse_hoz_shift(ctx, param, value: Optional[str]) -> tuple[int, int]:
    """Parse the --hoz-shift option, 'LON,LAT' cell counts, into a tuple."""
    if value is None:
        return (0, 0)
    try:
        shift_x, shift_y = (int(part) for part in value.split(","))
    except ValueError:
        raise click.BadParameter(
            f"expected two comma-separated integers 'LON,LAT', got {value!r}"
        )
    return (shift_x, shift_y)


def _parse_multiplier(ctx, param, value: tuple[str, ...]) -> dict[str, float]:
    """Parse the --multiplier options, 'VAR=COEFF', into a {var: coeff} dict."""
    scale_factors: dict[str, float] = {}
    for entry in value:
        if "=" not in entry:
            raise click.BadParameter(
                f"expected 'VAR=COEFF', got {entry!r}"
            )
        var, _, coeff = entry.partition("=")
        var = var.strip()
        try:
            coeff_value = float(coeff)
        except ValueError:
            raise click.BadParameter(
                f"coefficient must be a number, got {coeff!r} in {entry!r}"
            )
        if var in scale_factors:
            raise click.BadParameter(
                f"variable {var!r} given more than once"
            )
        scale_factors[var] = coeff_value
    return scale_factors


def backup_files(wrfinput: Path, wrfbdy: Optional[Path]):
    """
    Make a copy of wrfinput and wrfbdy in the same directory, adding a .orig suffix.
    If the .orig file already exists, this function will not overwrite it.
    `wrfbdy` is skipped if None
    """

    print("Backing up wrfinput/wrfbdy to wrfinput.orig/wrfbdy.orig")
    if wrfinput.with_suffix(".orig").exists():
        print("Backup file already exists, will not overwrite")
    else:
        shutil.copy(wrfinput, wrfinput.with_suffix(".orig"))

    if wrfbdy is not None:
        if wrfbdy.with_suffix(".orig").exists():
            print("Backup file already exists, will not overwrite")
        else:
            shutil.copy(wrfbdy, wrfbdy.with_suffix(".orig"))


def do_mappings(
    mapping: SpeciesMap, variables: xr.Dataset, shape: tuple[int, ...], squeeze=False
):
    """
    Returns a dictionary of name -> DataArray pairs with the mapped arrays.
    For each variable in the mapping, the function will:
    - Convert the source variable to the target variable's units
    - Apply the coefficients
    - Apply the weight and offset

    Args:
        - mapping: SpeciesMap object
        - vars: xarray.Dataset with the source variables
        - shape: Shape of the output arrays
        - squeeze: Whether to squeeze the input arrays before applying the mapping. Use when doing the boundary.
    """

    if squeeze:
        # Remove any 1-sized dimensions from the shape
        shape = tuple(s for s in shape if s != 1)

    out = {}
    for name, spec in mapping.map.items():
        out[name] = np.zeros(shape)

        for component, coefficient in mapping.map[name].coeffs.items():
            alias = mapping.aliases_source.get(component, component)
            arr = variables[alias]

            if squeeze:
                arr = arr.squeeze()

            assert arr.shape == shape, (
                f"Shape mismatch for {alias}: {arr.shape} != {shape}"
            )

            out[name] += (
                convert_si(
                    arr,
                    spec.units.global_model,
                    spec.units.regional_model,
                )
                * coefficient
            )
        out[name] = (out[name] * spec.weight) + spec.offset

    return out


def do_initial_conditions(
    wrf: WRFInput,
    wrf_ds: xr.Dataset,
    global_model_ds: xr.Dataset,
    mappings: SpeciesMap,
    write_diagnostics=False,
    below_surface: str = "clamp",
    scale_factors: Optional[dict[str, float]] = None,
):
    """Interpolate global model fields to WRF-Chem domain for initial conditions.

    Args:
        wrf: WRFInput object, used to write the interpolated fields to the `wrfinput` file
        wrf_ds: xarray.Dataset containing the WRF-Chem domain's basic coordinates and pressure field
        global_model_ds: xarray.Dataset containing the global model fields
        mappings: SpeciesMap object to map from global to WRF-CHEM species
        write_diagnostics: Whether to write out a diagnostic file for debugging purposes.
        below_surface: Policy for WRF half-levels below the global-model surface;
            see ``interpolate_to_wrf``.
        scale_factors: Optional {variable: coefficient} mapping. Matching output
            variables are multiplied by their coefficient before being written and
            the coefficient is recorded in the `interpolator_scale_factor` attribute.

    Returns:
        The set of `scale_factors` keys that matched a written variable.
    """
    scale_factors = scale_factors or {}
    matched: set[str] = set()

    interp_ds = interpolate_to_wrf(
        wrf_ds, global_model_ds, below_surface=below_surface
    )

    # Write interpolated source field to diagnostic file, if enabled.
    if write_diagnostics:
        interp_diag_ds = interp_ds.copy()
        interp_diag_ds["pres"] = wrf_ds["pres"]
        interp_diag_ds.to_netcdf("diag_interp.nc", "w")

    # Compute mappings
    wrf_vars = do_mappings(mappings, interp_ds, wrf_ds.pres.shape)

    # Write to WRF
    for name, arr in wrf_vars.items():
        alias = mappings.aliases_target.get(name, name)

        coeff = scale_factors.get(name)
        if coeff is not None:
            arr = arr * coeff
            matched.add(name)

        if alias not in wrf.nc_file.variables:
            wrf.nc_file.createVariable(
                alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
            )
        wrf.nc_file.variables[alias][0, :, :, :] = arr

        if coeff is not None:
            wrf.nc_file.variables[alias].setncattr("interpolator_scale_factor", coeff)

    return matched


def do_boundary_conditions(
    wrfbdy_path: Path,
    wrf_ds: xr.Dataset,
    global_model,
    mapping: SpeciesMap,
    skip_vertical: bool = False,
    below_surface: str = "clamp",
    hoz_shift: tuple[int, int] = (0, 0),
    scale_factors: Optional[dict[str, float]] = None,
):
    """Interpolate global model fields to the boundary of the WRF-Chem file and compute tendencies.

    Args:
        wrfbdy_path: Path to the WRF boundary file (wrfbdy_d01)
        wrf_ds: Dataset containing the WRF-CHEM domain's basic coordinates and pressure field
        global_model: GlobalModel object to fetch global model fields for each time
        mapping: SpeciesMap object to map from global to WRF-CHEM species
        hoz_shift: Horizontal shift (x, y) applied by `global_model`, recorded
            in the output file's root attributes for accounting
        scale_factors: Optional {variable: coefficient} mapping. Matching boundary
            variables are multiplied by their coefficient before being written
            (tendencies follow implicitly) and the coefficient is recorded in the
            `interpolator_scale_factor` attribute.

    Returns:
        The set of `scale_factors` keys that matched a written variable.
    """
    scale_factors = scale_factors or {}
    matched: set[str] = set()

    wrfbdy = WRFBoundary(wrfbdy_path)
    wrfbdy.nc_file.setncattr("interpolator_hoz_shift", f"{hoz_shift[0]},{hoz_shift[1]}")

    interp_last_t = {}
    for t_idx, t in enumerate(wrfbdy.times):
        # Fetch global model data for this specific time (temporally interpolated
        # if t falls between two available global-model times)
        global_model_ds = global_model.get_dataset_interpolated(t)
        if skip_vertical:
            global_model_ds.attrs["skip_vertical"] = True

        wrf_pres = wrf_ds["pres"]

        for bdy in ["BXS", "BXE", "BYS", "BYE"]:
            # Get WRF profile and replace pressure from met_em
            wrf_bdy = utils.get_boundary_profile(wrf_ds, bdy)
            wrf_bdy["pres"] = utils.get_boundary_profile(wrf_pres, bdy)

            # Interpolate to CAMS
            cams_bdy = interpolate_to_wrf(
                wrf_bdy, global_model_ds, below_surface=below_surface
            )

            # Do mappings
            # We use squeeze here because up until this point, all arrays keep their
            # original dimensions, even if they are of size 1.
            # This helps share the interpolation code between the initial and
            # boundary conditions, but now we don't need the extra dimension.
            wrf_vars = do_mappings(mapping, cams_bdy, wrf_bdy.pres.shape, squeeze=True)

            # Write to wrfbdy, except the last time step, which is only used for tendencies
            for name, arr in wrf_vars.items():
                alias = mapping.aliases_target.get(name, name) + f"_{bdy}"

                coeff = scale_factors.get(name)
                if coeff is not None:
                    arr = arr * coeff
                    matched.add(name)

                if t == wrfbdy.times[-1]:
                    interp_last_t[alias] = arr.to_numpy()
                    continue

                if alias not in wrfbdy.nc_file.variables:
                    wrfbdy.nc_file.createVariable(alias, "f4", ("Time", *arr.dims))
                wrfbdy.nc_file.variables[alias][t_idx, ...] = arr

                if coeff is not None:
                    wrfbdy.nc_file.variables[alias].setncattr(
                        "interpolator_scale_factor", coeff
                    )

    # Compute tendencies
    # For each boundary, store the difference between the current and previous value
    # inside the {var}_BT{bdy} variable, where {bdy} is one of XS, XE, YS, YE.
    # At this point, the `interp_last_t` variable contains the interpolated values for the last
    # timestep, which should be used as the "next" value for the last timestep.
    for t_idx, t in enumerate(wrfbdy.times):
        if t_idx == len(wrfbdy.times) - 1:
            continue  # No tendency for last time step

        dt = (wrfbdy.times[t_idx + 1] - t).total_seconds()

        for name in wrf_vars.keys():
            for bdy in ["XS", "XE", "YS", "YE"]:
                bdy_var = f"{name}_B{bdy}"
                bdy_t_var = f"{name}_BT{bdy}"

                if bdy_t_var not in wrfbdy.nc_file.variables:
                    wrfbdy.nc_file.createVariable(
                        bdy_t_var, "f4", wrfbdy.nc_file[bdy_var].dimensions
                    )

                curr_var = wrfbdy.nc_file.variables[bdy_var][t_idx, ...]
                next_var = (
                    wrfbdy.nc_file.variables[bdy_var][t_idx + 1, ...]
                    if t != wrfbdy.times[-2]
                    else interp_last_t[bdy_var]
                )

                wrfbdy.nc_file.variables[bdy_t_var][t_idx, ...] = (
                    next_var - curr_var
                ) / dt

    return matched


@click.command(name="interpolator-for-wrf")
@click.argument("global_model", type=click.Choice(list(GLOBAL_MODELS.keys())))
@click.argument(
    "input_files",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
)
@click.argument(
    "mapping",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.argument(
    "wrfinput",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.option(
    "--wrfbdy",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.option("--no-ic", is_flag=True)
@click.option("--copy-icbc", is_flag=True)
@click.option("--diagnostics", is_flag=True)
@click.option(
    "--below-surface-fill",
    type=click.Choice(["clamp", "extend"]),
    default="clamp",
    help=(
        "How to fill WRF half-levels that sit below the global model's surface "
        "(deepest WRF layer in low-elevation cells). 'clamp' is strictly "
        "mass-conservative (no extra mass below the source surface); 'extend' "
        "holds the near-surface mixing ratio constant down to the WRF surface."
    ),
)
@click.option(
    "--hoz-shift",
    default=None,
    callback=_parse_hoz_shift,
    help=(
        "Perturbation: shift the global model fields by whole grid cells "
        "before interpolation, as 'LON,LAT' (e.g. 2,-3). Positive values move "
        "the field east/north. Applied as a periodic roll of the data; the "
        "coordinate axes are unchanged."
    ),
)
@click.option(
    "--multiplier",
    multiple=True,
    callback=_parse_multiplier,
    help=(
        "Scale an output variable by a constant before writing it, as "
        "'VAR=COEFF' (e.g. DUST_1=0.98). May be given multiple times for "
        "different variables. The coefficient is recorded in the variable's "
        "'interpolator_scale_factor' netCDF attribute."
    ),
)
def main(
    global_model: str,
    input_files: Path,
    mapping: Path,
    wrfinput: Path,
    wrfbdy: Optional[Path],
    no_ic: bool,
    copy_icbc: bool,
    diagnostics: bool,
    below_surface_fill: str,
    hoz_shift: tuple[int, int],
    multiplier: dict[str, float],
):
    """
    Interpolate global model fields to WRF-Chem input files.

    Args:
        GLOBAL_MODEL: Global model to use
        INPUT_FILES: Path to directory containing the global model files. They must be
                     in subdirectories per day, in YYYY-MM-DD format.
        MAPPING: Path to the mapping file (toml) to use
        WRFINPUT: Path to the WRF input file

    Options:
        --wrfbdy: Path to the WRF boundary file to update, optional. If omitted, no boundary file is updated.
        --no-ic: Do not update initial conditions, wrfinput file still required for reading the grid.
        --copy-icbc: Copy wrfinput/wrfbdy to wrfinput.orig/wrfbdy.orig before modifying, will not overwrite any existing .orig files.
        --diagnostics: Write out diagnostic file for debugging. It will be written to the current directory with the name `diag_cams_interp.nc`. The diagnositc file contains the raw fields interpolated to the WRF grid, before the
        species mapping.
        --hoz-shift: Shift the global model fields by whole grid cells ('LON,LAT') before interpolation, as a perturbation. The applied shift is recorded in the `interpolator_hoz_shift` root attribute of the output files.
        --multiplier: Scale an output variable by a constant before writing, as 'VAR=COEFF' (e.g. DUST_1=0.98). May be repeated for different variables. The coefficient is recorded in the variable's `interpolator_scale_factor` attribute.
    """

    scale_factors = multiplier
    matched_scale_factors: set[str] = set()

    # Backup wrfinput/wrfbdy if requested
    if copy_icbc:
        backup_files(wrfinput, wrfbdy)

    # Read input files: mapping def., global model, WRF
    mapping = SpeciesMap(Path(mapping))
    print(mapping)

    wrf = WRFInput(wrfinput, read_only=no_ic)
    wrf_ds = wrf.get_dataset()
    print(wrf)

    global_model = GLOBAL_MODELS[global_model](
        Path(input_files), mapping.required_source_species
    )
    global_model.hoz_shift = hoz_shift
    if hoz_shift != (0, 0):
        print(f"Applying horizontal shift perturbation: {hoz_shift} cells (x, y)")
    print(global_model)

    # If the global model needs coordinate transformation (e.g., wrfout projected grid)
    if hasattr(global_model, "transform_target_coords"):
        wrf_ds = global_model.transform_target_coords(wrf_ds)

    # Detect matching vertical levels: if source and target share the same eta levels
    # (ZNU), vertical interpolation is unnecessary and skipped entirely.
    skip_vertical = False
    if hasattr(global_model, "znu") and "ZNU" in wrf_ds.coords:
        src_znu = global_model.znu
        dst_znu = wrf_ds.ZNU.values
        if len(src_znu) == len(dst_znu) and np.allclose(src_znu, dst_znu):
            skip_vertical = True
            print("Source and target WRF vertical levels match — skipping vertical interpolation")

    # Interpolate initial conditions (temporally interpolated if wrf.time falls
    # between two available global-model times)
    global_model_ds = global_model.get_dataset_interpolated(wrf.time)
    if skip_vertical:
        global_model_ds.attrs["skip_vertical"] = True

    # Initial conditions
    if not no_ic:
        matched_scale_factors |= do_initial_conditions(
            wrf,
            wrf_ds,
            global_model_ds,
            mapping,
            diagnostics,
            below_surface=below_surface_fill,
            scale_factors=scale_factors,
        )

        # Record the applied perturbation for accounting
        wrf.nc_file.setncattr(
            "interpolator_hoz_shift", f"{hoz_shift[0]},{hoz_shift[1]}"
        )

    # Compute boundary
    if wrfbdy:
        print(f"Doing boundary conditions ({wrfbdy})")
        matched_scale_factors |= do_boundary_conditions(
            wrfbdy,
            wrf_ds,
            global_model,
            mapping,
            skip_vertical=skip_vertical,
            below_surface=below_surface_fill,
            hoz_shift=hoz_shift,
            scale_factors=scale_factors,
        )

    # Warn about any --multiplier variables that never matched a written field
    unmatched = set(scale_factors) - matched_scale_factors
    if unmatched:
        print(
            "Warning: these --multiplier variables did not match any written "
            f"variable: {', '.join(sorted(unmatched))}"
        )

    wrf.close()
