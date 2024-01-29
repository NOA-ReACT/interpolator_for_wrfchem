import argparse
import shutil
from collections import namedtuple
from pathlib import Path

import numpy as np
from rich.progress import track

from interpolator_for_wrfchem.global_model_2 import CAMS_EAC4_ML
from interpolator_for_wrfchem.interpolation import interpolate_to_wrf
from interpolator_for_wrfchem.met_em import MetEm
from interpolator_for_wrfchem.species_map import SpeciesMap
from interpolator_for_wrfchem.wrf import WRF
from interpolator_for_wrfchem import utils


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

        if args.wrfbdy.with_suffix(".orig").exists():
            print("Backup file already exists, will not overwrite")
        else:
            shutil.copy(args.wrfbdy, args.wrfbdy.with_suffix(".orig"))

    # Read input files: mapping def., global model, WRF
    mapping = SpeciesMap(Path(args.mapping))
    print(mapping)

    wrf = WRF(Path(args.wrfinput), Path(args.wrfbdy))
    wrf_ds = wrf.get_dataset()
    print(wrf)

    met_em = MetEm(Path(args.met_em), wrf_ds)
    print(met_em)

    cams = CAMS_EAC4_ML(Path(args.input_files), mapping.required_source_species)
    print(cams)

    # Interpolate initial conditions
    if wrf.wrfinput_time not in cams.times:
        raise RuntimeError(
            f"Could not find global model file for wrfinput time {wrf.wrfinput_time}"
        )

    cams_ds = cams.get_dataset(
        wrf.wrfinput_time,
    )

    # Initial conditions
    if not args.no_ic:
        cams_interp_ds = interpolate_to_wrf(wrf_ds, cams_ds)

        if args.diagnostics:
            cams_diag_ds = cams_interp_ds.copy()
            cams_diag_ds["pres"] = wrf_ds["pres"]
            cams_interp_ds.to_netcdf("diag_cams_interp.nc", "w")

        # Compute mappings
        wrf_vars = {}
        for species in track(mapping.map.keys(), description="Mapping..."):
            wrf_vars[species] = np.zeros(wrf_ds.pres.shape)

            for component, coefficient in mapping.map[species].items():
                alias = mapping.aliases_source.get(component, component)
                wrf_vars[species] += (
                    cams_interp_ds[alias] * mapping.units.source * coefficient
                ) / mapping.units.target

        # Write to WRF
        for name, arr in track(wrf_vars.items(), description="Writing..."):
            alias = mapping.aliases_target.get(name, name)

            if alias not in wrf.wrfinput.variables:
                wrf.wrfinput.createVariable(
                    alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
                )
            wrf.wrfinput.variables[alias][0, :, :, :] = arr

    # Compute boundary
    if not args.no_bc:
        for t_idx, t in enumerate(wrf.wrfbdy_time):
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
                for species in mapping.map.keys():
                    wrf_vars[species] = np.zeros(wrf_bdy.pres.squeeze().shape)

                    for component, coefficient in mapping.map[species].items():
                        alias = mapping.aliases_source.get(component, component)
                        wrf_vars[species] += (
                            cams_bdy[alias].squeeze()
                            * mapping.units.source
                            * coefficient
                        ) / mapping.units.target

                # Write to wrfbdy
                for name, arr in wrf_vars.items():
                    alias = mapping.aliases_target.get(name, name) + f"_{bdy}"

                    if alias not in wrf.wrfbdy.variables:
                        wrf.wrfbdy.createVariable(alias, "f4", ("Time", *arr.dims))
                    wrf.wrfbdy.variables[alias][t_idx, ...] = arr

        # Compute tendencies
        # For each boundary, store the difference between the current and previous value
        # inside the {var}_BT{bdy} variable, where {bdy} is one of XS, XE, YS, YE.
        for t_idx, t in enumerate(wrf.wrfbdy_time):
            if t_idx == 0:
                continue
            dt = (t - wrf.wrfbdy_time[t_idx - 1]).total_seconds()

            for name in wrf_vars.keys():
                for bdy in ["XS", "XE", "YS", "YE"]:
                    bdy_var = f"{name}_B{bdy}"
                    bdy_t_var = f"{name}_BT{bdy}"

                    if bdy_t_var not in wrf.wrfbdy.variables:
                        wrf.wrfbdy.createVariable(
                            bdy_t_var, "f4", wrf.wrfbdy[bdy_var].dimensions
                        )

                    wrf.wrfbdy.variables[bdy_t_var][t_idx, ...] = (
                        wrf.wrfbdy.variables[bdy_var][t_idx, ...]
                        - wrf.wrfbdy.variables[bdy_var][t_idx - 1, ...]
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
    argparser.add_argument("wrfinput", type=Path, help="WRF input file to update")
    argparser.add_argument("wrfbdy", type=Path, help="WRF boundary file to update")
    argparser.add_argument("mapping", type=Path, help="Mapping file to use")
    argparser.add_argument(
        "--copy-icbc",
        "-c",
        action="store_true",
        help="Do not modify wrfinput/wrfbdy, instead copy them to wrfinput.orig/wrfbdy.orig",
    )
    argparser.add_argument(
        "--diagnostics",
        "-d",
        action="store_true",
        help="Write out diagnostic file for debugging. It will be written to the current directory with the name `diag.nc`.",
    )
    argparser.add_argument(
        "--no-ic",
        action="store_true",
        help="Do not interpolate initial conditions",
    )
    argparser.add_argument(
        "--no-bc",
        action="store_true",
        help="Do not interpolate boundary conditions",
    )
    args = vars(argparser.parse_args())

    args_type = namedtuple(
        "Args",
        args.keys(),
    )
    return args_type(**args)
