import argparse
import shutil
from collections import namedtuple
from pathlib import Path

import netCDF4 as nc
import numpy as np
from rich.progress import track
from scipy import interpolate

from interpolator_for_wrfchem.global_model import CAMS_EAC4
from interpolator_for_wrfchem.met_em import MetEm
from interpolator_for_wrfchem.species_map import SpeciesMap
from interpolator_for_wrfchem.wrf import WRF


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
    print(wrf)

    cams = CAMS_EAC4(Path(args.input_files))
    print(cams)

    # Open diagnostic file if requested
    diag_nc = None
    if args.diagnostics:
        diag_nc = nc.Dataset("diag.nc", "w")

    # Interpolate initial conditions
    if wrf.wrfinput_time not in cams.times:
        raise RuntimeError(
            f"Could not find global model file for wrfinput time {wrf.wrfinput_time}"
        )

    wrf_xlong = wrf.wrfinput.variables["XLONG"][0, :, :]
    wrf_xlat = wrf.wrfinput.variables["XLAT"][0, :, :]
    wrf_pres = (
        wrf.wrfinput.variables["P"][0, :, :, :]
        + wrf.wrfinput.variables["PB"][0, :, :, :]
    ) * 0.01
    cams_pres = cams.get_var(wrf.wrfinput_time, "level")

    cams_interp = {}
    for var in track(mapping.required_source_species, description="Interpolating..."):
        cams_interp[var] = np.empty(wrf_pres.shape)

        v = cams.interpolate_hoz(wrf.wrfinput_time, var, wrf_xlong, wrf_xlat)

        # Interpolate vertical levels
        for (x, y), _ in np.ndenumerate(wrf_xlong):
            cams_interp[var][:, x, y] = interpolate.interp1d(
                cams_pres,
                v[:, x, y],
                kind="linear",
                fill_value="extrapolate",
            )(wrf_pres[:, x, y])

    if diag_nc is not None:
        # Write interpolated fields to diagnostic file
        diag_nc.createDimension("south_north", wrf.size_south_north)
        diag_nc.createDimension("west_east", wrf.size_west_east)
        diag_nc.createDimension("bottom_top", wrf.size_bottom_top)

        for name, arr in cams_interp.items():
            diag_nc.createVariable(
                name, "f4", ("bottom_top", "south_north", "west_east")
            )
            diag_nc.variables[name][:, :, :] = arr

    # Compute mappings
    wrf_vars = {}
    for species in track(mapping.map.keys(), description="Mapping..."):
        wrf_vars[species] = np.zeros(wrf_pres.shape)

        for component, coefficient in mapping.map[species].items():
            alias = mapping.aliases_source.get(component, component)
            wrf_vars[species] += cams_interp[alias] * coefficient

    # Write to WRF
    for name, arr in track(wrf_vars.items(), description="Writing..."):
        alias = mapping.aliases_target.get(name, name)

        if alias not in wrf.wrfinput.variables:
            wrf.wrfinput.createVariable(
                alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
            )
        wrf.wrfinput.variables[name][0, :, :, :] = arr
        # TODO Units?

    # Compute boundary

    wrf.close()
    if diag_nc is not None:
        diag_nc.close()


def cmd_args():
    argparser = argparse.ArgumentParser(
        description="Interpolate WRF-Chem output to a regular grid"
    )
    argparser.add_argument(
        "input_files", type=Path, help="Global model fields to interpolate"
    )
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
    args = vars(argparser.parse_args())

    args_type = namedtuple(
        "Args",
        args.keys(),
    )
    return args_type(**args)
