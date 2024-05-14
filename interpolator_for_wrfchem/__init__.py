import argparse
import shutil
from collections import namedtuple
from pathlib import Path

import numpy as np

from interpolator_for_wrfchem import utils
from interpolator_for_wrfchem.global_model import CAMS_EAC4_ML
from interpolator_for_wrfchem.interpolation import interpolate_to_wrf
from interpolator_for_wrfchem.met_em import MetEm
from interpolator_for_wrfchem.species_map import SpeciesMap, convert_si
from interpolator_for_wrfchem.wrf import WRFBoundary, WRFInput


def main():
    """Entrypoint, parse arguments and run script"""
    args = cmd_args()

    # Backup wrfinput/wrfbdy if requested
    if args.copy_icbc:
        print("Backing up wrfinput/wrfbdy to wrfinput.orig/wrfbdy.orig")
        if args.wrfinput.with_suffix(".orig").exists():
            print("Backup file already exists, will not overwrite")
        else:
            shutil.copy(args.wrfinput, args.wrfinput.with_suffix(".orig"))

        if "wrfbdy" in args and args.wrfbdy is not None:
            if args.wrfbdy.with_suffix(".orig").exists():
                print("Backup file already exists, will not overwrite")
            else:
                shutil.copy(args.wrfbdy, args.wrfbdy.with_suffix(".orig"))

    # Read input files: mapping def., global model, WRF
    mapping = SpeciesMap(Path(args.mapping))
    print(mapping)

    wrf = WRFInput(args.wrfinput)
    wrf_ds = wrf.get_dataset()
    print(wrf)

    met_em = MetEm(Path(args.met_em), wrf_ds)
    print(met_em)

    cams = CAMS_EAC4_ML(Path(args.input_files), mapping.required_source_species)
    print(cams)

    # Interpolate initial conditions
    if wrf.time not in cams.times:
        raise RuntimeError(
            f"Could not find global model file for wrfinput time {wrf.nc_file_time}"
        )

    cams_ds = cams.get_dataset(wrf.time)

    # Initial conditions
    if not args.no_ic:
        cams_interp_ds = interpolate_to_wrf(wrf_ds, cams_ds)

        if args.diagnostics:
            cams_diag_ds = cams_interp_ds.copy()
            cams_diag_ds["pres"] = wrf_ds["pres"]
            cams_interp_ds.to_netcdf("diag_cams_interp.nc", "w")

        # Compute mappings
        wrf_vars = {}
        for name, spec in mapping.map.items():
            wrf_vars[name] = np.zeros(wrf_ds.pres.shape)

            for component, coefficient in mapping.map[name].coeffs.items():
                alias = mapping.aliases_source.get(component, component)
                wrf_vars[name] += (
                    convert_si(
                        cams_interp_ds[alias],
                        spec.units.global_model,
                        spec.units.regional_model,
                    )
                    * coefficient
                )
            wrf_vars[name] = (wrf_vars[name] * spec.weight) + spec.offset

        # Write to WRF
        for name, arr in wrf_vars.items():
            alias = mapping.aliases_target.get(name, name)

            if alias not in wrf.nc_file.variables:
                wrf.nc_file.createVariable(
                    alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
                )
            wrf.nc_file.variables[alias][0, :, :, :] = arr

    # Compute boundary
    if "wrfbdy" in args:
        wrfbdy = WRFBoundary(args.wrfbdy)
        for t_idx, t in enumerate(wrfbdy.times):
            wrf_pres = met_em.get_pres(t)

            for bdy in ["BXS", "BXE", "BYS", "BYE"]:
                # Get WRF profile and replace pressure from met_em
                wrf_bdy = utils.get_boundary_profile(wrf_ds, bdy)
                wrf_bdy["pres"] = utils.get_boundary_profile(wrf_pres, bdy)

                # Interpolate to CAMS
                cams_bdy = interpolate_to_wrf(wrf_bdy, cams_ds)

                # Do mappings
                # We use squeeze() at various points here because up until this point,
                # all arrays keep their original dimensions, even if they are of size 1.
                # This helps share the interpolation code between the initial and
                # boundary conditions, but now we don't need the extra dimension.
                wrf_vars = {}
                for name, spec in mapping.map.items():
                    wrf_vars[name] = np.zeros(wrf_bdy.pres.squeeze().shape)

                    for component, coefficient in mapping.map[name].coeffs.items():
                        alias = mapping.aliases_source.get(component, component)
                        wrf_vars[name] += (
                            convert_si(
                                cams_bdy[alias].squeeze(),
                                spec.units.global_model,
                                spec.units.regional_model,
                            )
                            * coefficient
                        )
                    wrf_vars[name] = (wrf_vars[name] * spec.weight) + spec.offset

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

    wrf.close()


def cmd_args():
    argparser = argparse.ArgumentParser(
        description="Interpolate WRF-Chem output to a regular grid"
    )
    argparser.add_argument(
        "input_files", type=Path, help="Global model fields to interpolate"
    )
    argparser.add_argument("met_em", type=Path, help="Directory of met_em files")
    argparser.add_argument("mapping", type=Path, help="Mapping file to use")

    argparser.add_argument("wrfinput", type=Path, help="WRF input file")
    argparser.add_argument("--wrfbdy", type=Path, help="WRF boundary file to update")

    argparser.add_argument(
        "--no-ic",
        action="store_true",
        help="Do not update initial conditions, wrfinput file still required.",
    )
    argparser.add_argument(
        "--copy-icbc",
        action="store_true",
        help="Copy wrfinput/wrfbdy to wrfinput.orig/wrfbdy.orig before modifying, will not overwrite any existing .orig files",
    )
    argparser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Write out diagnostic file for debugging. It will be written to the current directory with the name `diag_cams_interp.nc`.",
    )
    args = vars(argparser.parse_args())

    args_type = namedtuple(
        "Args",
        args.keys(),
    )
    return args_type(**args)
