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

    met_em = MetEm(Path(args.met_em), wrf.wrfinput)
    print(met_em)

    cams = CAMS_EAC4_ML(Path(args.input_files), mapping.required_source_species)
    print(cams)

    # Interpolate initial conditions
    if wrf.wrfinput_time not in cams.times:
        raise RuntimeError(
            f"Could not find global model file for wrfinput time {wrf.wrfinput_time}"
        )

    wrf_ds = wrf.get_dataset()
    cams_ds = cams.get_dataset(
        wrf.wrfinput_time,
    )
    cams_interp_ds = interpolate_to_wrf(wrf_ds, cams_ds)

    if args.diagnostics:
        cams_interp_ds.to_netcdf("diag.nc", "w")

    # Compute mappings
    wrf_vars = {}
    for species in track(mapping.map.keys(), description="Mapping..."):
        wrf_vars[species] = np.zeros(wrf_ds.pres.shape)

        for component, coefficient in mapping.map[species].items():
            alias = mapping.aliases_source.get(component, component)
            wrf_vars[species] += cams_interp_ds[alias] * coefficient

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
