"""
Time series of global-model fields on the WRF grid, written to a file of their own.

`interpolator-for-wrfchem` writes the initial and boundary conditions into the
wrfinput/wrfbdy pair. Some fields are needed at other times as well, for example the
background oxidants of the GOCART sulphur chemistry (BACKG_OH, BACKG_H2O2, BACKG_NO3):
WRF reads them once from wrfinput and keeps them, also through restarts, so a cycling
system has to overwrite them every cycle. This module interpolates the fields of a
species map to the WRF grid at a series of times and writes them all to one netCDF
file (dimensions Time, bottom_top, south_north, west_east, plus a WRF-style `Times`
variable), from which a cycling system picks the right record.

The vertical remap uses the pressure of the wrfinput for all times (its column mass at
the wrfinput time). That is good enough for slowly varying fields such as oxidant
daily means, but not for the transported tracers.
"""

import datetime as dt
from pathlib import Path
from typing import Optional

import click
import netCDF4 as nc
import numpy as np
import xarray as xr

from interpolator_for_wrfchem import do_mappings
from interpolator_for_wrfchem.global_models import GLOBAL_MODELS
from interpolator_for_wrfchem.interpolation import interpolate_to_wrf
from interpolator_for_wrfchem.species_map import SpeciesMap
from interpolator_for_wrfchem.wrf import WRFInput

WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"


def daily_mean_dataset(global_model, times: list[dt.datetime]) -> xr.Dataset:
    """Average the global-model datasets at `times`, all fields including pressure.

    The pressure fields are averaged too, so the vertical remap sees a consistent
    mean column. Every dataset of a model shares its grid, so a plain sum works.
    """
    total = None
    for t in times:
        ds = global_model.get_dataset(t)
        total = ds if total is None else total + ds
    with xr.set_options(keep_attrs=True):
        return total / len(times)


def group_by_day(
    times: list[dt.datetime],
) -> dict[dt.date, list[dt.datetime]]:
    """Group the available times by UTC day, [00:00, 24:00)."""
    days: dict[dt.date, list[dt.datetime]] = {}
    for t in sorted(times):
        days.setdefault(t.date(), []).append(t)
    return days


def select_days(
    available: list[dt.datetime],
    start: Optional[dt.date],
    end: Optional[dt.date],
    allow_partial: bool,
) -> dict[dt.date, list[dt.datetime]]:
    """
    The days to average, each with its source times.

    A day is complete when it has as many source times as the best covered day (e.g.
    4 for 6-hourly data). Incomplete days, usually the last one when the data stops at
    00 UTC, are an error inside [start, end] unless `allow_partial`; outside an
    explicit range they are dropped with a message.
    """
    days = group_by_day(available)
    if start is not None:
        days = {d: ts for d, ts in days.items() if d >= start}
    if end is not None:
        days = {d: ts for d, ts in days.items() if d <= end}
    if not days:
        raise click.ClickException("No global-model times in the requested days")

    n_full = max(len(ts) for ts in days.values())
    partial = {d: ts for d, ts in days.items() if len(ts) < n_full}
    for d, ts in partial.items():
        hours = ", ".join(f"{t:%H}" for t in ts)
        message = f"{d} has {len(ts)} of {n_full} source times ({hours} UTC)"
        if allow_partial:
            print(f"Warning: {message}, averaging what is there")
        elif start is None and end is None:
            print(f"Skipping {message}")
            del days[d]
        else:
            raise click.ClickException(
                f"{message}. Download the missing times, narrow --start/--end or pass "
                "--allow-partial"
            )

    if start is not None and end is not None:
        wanted = {start + dt.timedelta(days=i) for i in range((end - start).days + 1)}
        missing = sorted(wanted - set(days))
        if missing:
            raise click.ClickException(
                "No global-model data for " + ", ".join(str(d) for d in missing)
            )

    return dict(sorted(days.items()))


class FieldsFile:
    """Output file: WRF-like dimensions and grid attributes copied from wrfinput."""

    def __init__(self, path: Path, wrf: WRFInput, description: str):
        self.nc = nc.Dataset(str(path), "w", format="NETCDF4")
        src = wrf.nc_file
        for name in src.ncattrs():
            self.nc.setncattr(name, src.getncattr(name))
        self.nc.setncattr("TITLE", "interpolator-for-wrfchem fields")
        self.nc.setncattr("interpolator_description", description)

        self.nc.createDimension("Time", None)
        self.nc.createDimension("DateStrLen", 19)
        for dim in ("bottom_top", "south_north", "west_east"):
            self.nc.createDimension(dim, src.dimensions[dim].size)

        self.times = self.nc.createVariable("Times", "S1", ("Time", "DateStrLen"))
        for name in ("XLAT", "XLONG"):
            var = self.nc.createVariable(name, "f4", ("south_north", "west_east"))
            var[:] = src.variables[name][0]
            var.units = src.variables[name].units

        self.n = 0

    def write(self, time: dt.datetime, fields: dict[str, np.ndarray], attrs: dict):
        """Append one record. `attrs` holds per-variable attributes, keyed by name."""
        self.times[self.n] = np.frombuffer(
            time.strftime(WRF_TIME_FORMAT).encode("ascii"), dtype="S1"
        )
        for name, arr in fields.items():
            if name not in self.nc.variables:
                var = self.nc.createVariable(
                    name,
                    "f4",
                    ("Time", "bottom_top", "south_north", "west_east"),
                    zlib=True,
                )
                var.FieldType = np.int32(104)
                var.MemoryOrder = "XYZ"
                var.stagger = ""
                for key, value in attrs.get(name, {}).items():
                    var.setncattr(key, value)
            self.nc.variables[name][self.n] = np.asarray(arr, dtype="f4")
        self.n += 1

    def close(self):
        self.nc.close()


def _parse_date(ctx, param, value: Optional[str]) -> Optional[dt.date]:
    if value is None:
        return None
    try:
        return dt.date.fromisoformat(value)
    except ValueError:
        raise click.BadParameter(f"expected YYYY-MM-DD, got {value!r}")


@click.command(name="interpolator-for-wrfchem-fields")
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
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--daily-mean",
    "mode",
    flag_value="daily-mean",
    default=True,
    help="One record per UTC day, the mean of the source times in [00, 24) UTC "
    "(default).",
)
@click.option(
    "--instantaneous",
    "mode",
    flag_value="instantaneous",
    help="One record per source time, no averaging.",
)
@click.option("--start", callback=_parse_date, help="First day, YYYY-MM-DD.")
@click.option("--end", callback=_parse_date, help="Last day, YYYY-MM-DD (inclusive).")
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Average days that miss some source times instead of stopping.",
)
@click.option(
    "--below-surface-fill",
    type=click.Choice(["clamp", "extend"]),
    default="extend",
    help="See interpolator-for-wrfchem. Default 'extend': these are concentration-"
    "like fields, not tracer mass that has to be conserved.",
)
def main(
    global_model: str,
    input_files: Path,
    mapping: Path,
    wrfinput: Path,
    output: Path,
    mode: str,
    start: Optional[dt.date],
    end: Optional[dt.date],
    allow_partial: bool,
    below_surface_fill: str,
):
    """
    Interpolate the fields of MAPPING at many times and write them to OUTPUT.

    Arguments are as for interpolator-for-wrfchem, plus the OUTPUT netCDF path. The
    wrfinput only provides the grid and the pressure, it is not modified.

    With --daily-mean (default) each record is the mean of one UTC day and its `Times`
    is that day at 00:00:00. With --instantaneous each source time is one record.
    """
    species_map = SpeciesMap(mapping)
    print(species_map)

    wrf = WRFInput(wrfinput, read_only=True)
    wrf_ds = wrf.get_dataset()
    print(wrf)

    model = GLOBAL_MODELS[global_model](
        Path(input_files), species_map.required_source_species
    )
    print(model)

    if mode == "daily-mean":
        days = select_days(model.available_times, start, end, allow_partial)
        records = [
            (dt.datetime.combine(d, dt.time()), ts) for d, ts in days.items()
        ]
        description = (
            "Daily means: each record is the mean of the source times in "
            "[Times, Times + 1 day)"
        )
    else:
        times = [
            t
            for t in model.available_times
            if (start is None or t.date() >= start) and (end is None or t.date() <= end)
        ]
        if not times:
            raise click.ClickException("No global-model times in the requested days")
        records = [(t, [t]) for t in times]
        description = "Instantaneous fields at the source times"

    out = FieldsFile(output, wrf, description)
    out.nc.setncattr("interpolator_source", f"{global_model} {Path(input_files).resolve()}")
    out.nc.setncattr("interpolator_species_map", str(Path(mapping).resolve()))
    out.nc.setncattr("interpolator_pressure_from", str(Path(wrfinput).resolve()))

    attrs = {
        species_map.aliases_target.get(name, name): {
            "description": f"{name} from species map '{species_map.name}'",
            "interpolator_coeffs": ", ".join(
                f"{k}={v}" for k, v in spec.coeffs.items()
            ),
        }
        for name, spec in species_map.map.items()
    }

    for label, source_times in records:
        hours = ", ".join(f"{t:%Y-%m-%d %H}" for t in source_times)
        print(f"{label:%Y-%m-%d %H:%M}: {len(source_times)} source time(s) ({hours})")
        gm_ds = daily_mean_dataset(model, source_times)
        interp = interpolate_to_wrf(wrf_ds, gm_ds, below_surface=below_surface_fill)
        fields = do_mappings(species_map, interp, wrf_ds.pres.shape)
        fields = {
            species_map.aliases_target.get(name, name): arr
            for name, arr in fields.items()
        }
        for name in fields:
            attrs.setdefault(name, {})["interpolator_n_source_times"] = np.int32(
                len(source_times)
            )
        out.write(label, fields, attrs)

    out.close()
    wrf.close()
    print(f"Wrote {len(records)} record(s) to {output}")
