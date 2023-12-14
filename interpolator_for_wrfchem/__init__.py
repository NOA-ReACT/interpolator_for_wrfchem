import argparse
from collections import namedtuple
from pathlib import Path
import shutil

import numpy as np

from interpolator_for_wrfchem.global_model import CAMS_EAC4
from interpolator_for_wrfchem.species_map import SpeciesMap
from interpolator_for_wrfchem.wrf import WRF
from rich.progress import track

from scipy import interpolate


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

        wrf.wrfinput.createVariable(
            alias, "f4", ("Time", "bottom_top", "south_north", "west_east")
        )
        wrf.wrfinput.variables[name][0, :, :, :] = arr
        # TODO Units?

    # Compute boundary

    wrf.close()


def cmd_args():
    argparser = argparse.ArgumentParser(
        description="Interpolate WRF-Chem output to a regular grid"
    )
    argparser.add_argument("input_files", help="Global model fields to interpolate")
    argparser.add_argument("wrfinput", help="WRF input file to update")
    argparser.add_argument("wrfbdy", help="WRF boundary file to update")
    argparser.add_argument("mapping", help="Mapping file to use")
    argparser.add_argument(
        "--copy-icbc",
        "-c",
        action="store_true",
        help="Do not modify wrfinput/wrfbdy, instead copy them to wrfinput_orig/wrfbdy_orig",
    )
    args = argparser.parse_args()

    args_type = namedtuple(
        "Args",
        ["input_files", "wrfinput", "wrfbdy", "mapping", "copy_icbc"],
    )
    return args_type(
        input_files=Path(args.input_files),
        wrfinput=Path(args.wrfinput),
        wrfbdy=Path(args.wrfbdy),
        mapping=Path(args.mapping),
        copy_icbc=args.copy_icbc,
    )
