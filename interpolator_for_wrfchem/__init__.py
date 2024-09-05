import shutil
from pathlib import Path
from typing import Optional

import click
import numpy as np
import xarray as xr

from interpolator_for_wrfchem import utils
from interpolator_for_wrfchem.global_models import GLOBAL_MODELS
from interpolator_for_wrfchem.interpolation import interpolate_to_wrf
from interpolator_for_wrfchem.met_em import MetEm
from interpolator_for_wrfchem.species_map import SpeciesMap, convert_si
from interpolator_for_wrfchem.wrf import WRFBoundary, WRFInput


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

            assert (
                arr.shape == shape
            ), f"Shape mismatch for {alias}: {arr.shape} != {shape}"

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
):
    """Interpolate global model fields to WRF-Chem domain for initial conditions.

    Args:
        wrf: WRFInput object, used to write the interpolated fields to the `wrfinput` file
        wrf_ds: xarray.Dataset containing the WRF-Chem domain's basic coordinates and pressure field
        global_model_ds: xarray.Dataset containing the global model fields
        mappings: SpeciesMap object to map from global to WRF-CHEM species
        write_diagnostics: Whether to write out a diagnostic file for debugging purposes.
    """

    interp_ds = interpolate_to_wrf(wrf_ds, global_model_ds)

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

        if alias not in wrf.nc_file.variables:
            wrf.nc_file.createVariable(
                alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
            )
        wrf.nc_file.variables[alias][0, :, :, :] = arr


def do_boundary_conditions(
    wrfbdy_path: Path,
    met_em: MetEm,
    wrf_ds: xr.Dataset,
    global_model_ds: xr.Dataset,
    mapping: SpeciesMap,
):
    """Interpolate global model fields to the boundary of the WRF-Chem file and compute tendencies.

    Args:
        wrfbdy_path: Path to the WRF boundary file (wrfbdy_d01)
        met_em: MetEm object representing the met_em files
        wrf_ds: Dataset containing the WRF-CHEM domain's basic coordinates and pressure field
        global_model: Dataset containing the global model fields
        mapping: SpeciesMap object to map from global to WRF-CHEM species
    """

    wrfbdy = WRFBoundary(wrfbdy_path)

    for t_idx, t in enumerate(wrfbdy.times):
        wrf_pres = met_em.get_pres(t)

        for bdy in ["BXS", "BXE", "BYS", "BYE"]:
            # Get WRF profile and replace pressure from met_em
            wrf_bdy = utils.get_boundary_profile(wrf_ds, bdy)
            wrf_bdy["pres"] = utils.get_boundary_profile(wrf_pres, bdy)

            # Interpolate to CAMS
            cams_bdy = interpolate_to_wrf(wrf_bdy, global_model_ds)

            # Do mappings
            # We use squeeze here because up until this point, all arrays keep their
            # original dimensions, even if they are of size 1.
            # This helps share the interpolation code between the initial and
            # boundary conditions, but now we don't need the extra dimension.
            wrf_vars = do_mappings(mapping, cams_bdy, wrf_bdy.pres.shape, squeeze=True)

            # Write to wrfbdy
            for name, arr in wrf_vars.items():
                alias = mapping.aliases_target.get(name, name) + f"_{bdy}"

                if alias not in wrfbdy.nc_file.variables:
                    wrfbdy.nc_file.createVariable(alias, "f4", ("Time", *arr.dims))
                wrfbdy.nc_file.variables[alias][t_idx, ...] = arr

        # Compute tendencies
        # For each boundary, store the difference between the current and previous value
        # inside the {var}_BT{bdy} variable, where {bdy} is one of XS, XE, YS, YE.
        for t_idx, t in enumerate(wrfbdy.times):
            if t_idx == 0:
                continue
            dt = (t - wrfbdy.times[t_idx - 1]).total_seconds()

            for name in wrf_vars.keys():
                for bdy in ["XS", "XE", "YS", "YE"]:
                    bdy_var = f"{name}_B{bdy}"
                    bdy_t_var = f"{name}_BT{bdy}"

                    if bdy_t_var not in wrfbdy.nc_file.variables:
                        wrfbdy.nc_file.createVariable(
                            bdy_t_var, "f4", wrfbdy.nc_file[bdy_var].dimensions
                        )

                    wrfbdy.nc_file.variables[bdy_t_var][t_idx, ...] = (
                        wrfbdy.nc_file.variables[bdy_var][t_idx, ...]
                        - wrfbdy.nc_file.variables[bdy_var][t_idx - 1, ...]
                    ) / dt


@click.command(name="interpolator-for-wrf")
@click.argument("global_model", type=click.Choice(list(GLOBAL_MODELS.keys())))
@click.argument(
    "input_files",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
)
@click.argument(
    "met_em",
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
def main(
    global_model: str,
    input_files: Path,
    met_em: Path,
    mapping: Path,
    wrfinput: Path,
    wrfbdy: Optional[Path],
    no_ic: bool,
    copy_icbc: bool,
    diagnostics: bool,
):
    """
    Interpolate global model fields to WRF-Chem input files.

    Args:
        GLOBAL_MODEL: Global model to use
        INPUT_FILES: Path to directory containing the global model files. They must be
                     in subdirectories per day, in YYYY-MM-DD foramt.
        MET_EM: Path to the directory containing the met_em files.
        MAPPING: Path to the mapping file (toml) to use
        WRFINPUT: Path to the WRF input file

    Options:
        --wrfbdy: Path to the WRF boundary file to update, optional. If omitted, no boundary file is updated.
        --no-ic: Do not update initial conditions, wrfinput file still required for reading the grid.
        --copy-icbc: Copy wrfinput/wrfbdy to wrfinput.orig/wrfbdy.orig before modifying, will not overwrite any existing .orig files.
        --diagnostics: Write out diagnostic file for debugging. It will be written to the current directory with the name `diag_cams_interp.nc`. The diagnositc file contains the raw fields interpolated to the WRF grid, before the
        species mapping.
    """

    # Backup wrfinput/wrfbdy if requested
    if copy_icbc:
        backup_files(wrfinput, wrfbdy)

    # Read input files: mapping def., global model, WRF
    mapping = SpeciesMap(Path(mapping))
    print(mapping)

    wrf = WRFInput(wrfinput)
    wrf_ds = wrf.get_dataset()
    print(wrf)

    met_em = MetEm(Path(met_em), wrf_ds)
    print(met_em)

    global_model = GLOBAL_MODELS[global_model](
        Path(input_files), mapping.required_source_species
    )
    print(global_model)

    # Interpolate initial conditions
    if wrf.time not in global_model.times:
        raise RuntimeError(
            f"Could not find global model file for wrfinput time {wrf.nc_file_time}"
        )

    global_model_ds = global_model.get_dataset(wrf.time)

    # Initial conditions
    if not no_ic:
        do_initial_conditions(wrf, wrf_ds, global_model_ds, mapping, diagnostics)

    # Compute boundary
    if wrfbdy:
        print(f"Doing boundary conditions ({wrfbdy})")
        do_boundary_conditions(wrfbdy, met_em, wrf_ds, global_model_ds, mapping)

    wrf.close()
