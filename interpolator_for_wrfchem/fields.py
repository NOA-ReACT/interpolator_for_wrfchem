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
WRF_TITLE = "OUTPUT FROM interpolator-for-wrfchem FOR WRF V4.x"


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
        src = wrf.nc_file
        self.nc = nc.Dataset(str(path), "w", format=src.data_model)
        for name in src.ncattrs():
            self.nc.setncattr(name, src.getncattr(name))
        # WRF only reads input files whose TITLE contains " V4." (share/input_wrf.F)
        self.nc.setncattr("TITLE", WRF_TITLE)
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
                    zlib=self.nc.data_model.startswith("NETCDF4"),
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


def _parse_time(ctx, param, value: Optional[str]) -> Optional[dt.datetime]:
    if value is None:
        return None
    for fmt in ("%Y-%m-%d_%H:%M:%S", "%Y-%m-%d_%H", "%Y-%m-%dT%H", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise click.BadParameter(f"expected YYYY-MM-DD or YYYY-MM-DD_HH, got {value!r}")


class InterpolatedFields:
    """Mapped fields on the WRF grid at any time, from the global model."""

    def __init__(self, model, species_map, wrf_ds, below_surface: str):
        self.model = model
        self.species_map = species_map
        self.wrf_ds = wrf_ds
        self.below_surface = below_surface

    def from_dataset(self, gm_ds) -> dict[str, np.ndarray]:
        interp = interpolate_to_wrf(self.wrf_ds, gm_ds, below_surface=self.below_surface)
        fields = do_mappings(self.species_map, interp, self.wrf_ds.pres.shape)
        return {
            self.species_map.aliases_target.get(name, name): np.asarray(arr)
            for name, arr in fields.items()
        }

    def at(self, t: dt.datetime) -> dict[str, np.ndarray]:
        """Instantaneous, linearly interpolated in time between source times."""
        return self.from_dataset(self.model.get_dataset_interpolated(t))


def write_reference_files(
    fields: InterpolatedFields,
    wrf: WRFInput,
    out_dir: Path,
    start: dt.datetime,
    end: dt.datetime,
    interval: int,
    file_template: str,
    attrs: dict,
    global_attrs: dict,
) -> int:
    """
    One WRF auxinput file per `interval` minutes from `start` to `end` (inclusive), each
    holding every field of the map twice: `<NAME>_NDG_OLD` at the file time and
    `<NAME>_NDG_NEW` one interval later, as WRF-ReACT's chem nudging (and grid FDDA)
    reads them. So the data has to reach `end` + `interval`.
    """
    step = dt.timedelta(minutes=interval)
    times = []
    t = start
    while t <= end:
        times.append(t)
        t += step
    if not times:
        raise click.ClickException("--end is before --start")

    out_dir.mkdir(parents=True, exist_ok=True)
    cache: dict[dt.datetime, dict[str, np.ndarray]] = {}

    def get(t):
        if t not in cache:
            print(f"  interpolating {t:%Y-%m-%d %H:%M}")
            cache[t] = fields.at(t)
            for old in [k for k in cache if k < t - step]:
                del cache[old]
        return cache[t]

    for t in times:
        old, new = get(t), get(t + step)
        name = file_template.format(date=t.strftime(WRF_TIME_FORMAT))
        out = FieldsFile(out_dir / name, wrf, f"Reference fields: *_NDG_OLD at Times, "
                         f"*_NDG_NEW {interval} min later")
        for k, v in global_attrs.items():
            out.nc.setncattr(k, v)
        record = {}
        for var, arr in old.items():
            record[f"{var}_NDG_OLD"] = arr
            record[f"{var}_NDG_NEW"] = new[var]
        rec_attrs = {}
        for var in old:
            for suffix, when in (("OLD", "the file time"), ("NEW", f"{interval} min after the file time")):
                rec_attrs[f"{var}_NDG_{suffix}"] = {
                    **attrs.get(var, {}),
                    "description": f"{var} reference at {when}",
                }
        out.write(t, record, rec_attrs)
        out.close()
        print(f"{name}: OLD {t:%Y-%m-%d %H:%M}, NEW {t + step:%Y-%m-%d %H:%M}")
    return len(times)


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
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--daily-mean",
    "mode",
    flag_value="daily-mean",
    default=True,
    help="One record per UTC day, the mean of the source times in [00, 24) UTC "
    "(default). OUTPUT is a file.",
)
@click.option(
    "--instantaneous",
    "mode",
    flag_value="instantaneous",
    help="One record per source time, no averaging. OUTPUT is a file.",
)
@click.option(
    "--wrf-reference",
    "mode",
    flag_value="wrf-reference",
    help="WRF auxinput reference files for chem nudging: one file per --interval, "
    "with <NAME>_NDG_OLD at the file time and <NAME>_NDG_NEW one interval later. "
    "OUTPUT is a directory; --start and --end are required.",
)
@click.option(
    "--start",
    callback=_parse_time,
    help="First day (YYYY-MM-DD), or first file time for --wrf-reference "
    "(YYYY-MM-DD_HH).",
)
@click.option(
    "--end",
    callback=_parse_time,
    help="Last day, or last file time for --wrf-reference (inclusive).",
)
@click.option(
    "--interval",
    default=360,
    show_default=True,
    help="--wrf-reference: minutes between files, WRF's auxinput19_interval_m.",
)
@click.option(
    "--file-template",
    default="auxinput19_d01_{date}",
    show_default=True,
    help="--wrf-reference: file name, {date} becomes the WRF date of the file. Must "
    "match WRF's auxinput19_inname.",
)
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Average days that miss some source times instead of stopping.",
)
@click.option(
    "--below-surface-fill",
    type=click.Choice(["clamp", "extend"]),
    default=None,
    help="See interpolator-for-wrfchem. Default 'extend' for --daily-mean and "
    "--instantaneous (WRF layers below the global model's surface get the lowest "
    "source value, not zero oxidants) and 'clamp' for --wrf-reference, the default "
    "of interpolator-for-wrfchem, so the reference matches the tracer's IC/BC.",
)
def main(
    global_model: str,
    input_files: Path,
    mapping: Path,
    wrfinput: Path,
    output: Path,
    mode: str,
    start: Optional[dt.datetime],
    end: Optional[dt.datetime],
    interval: int,
    file_template: str,
    allow_partial: bool,
    below_surface_fill: Optional[str],
):
    """
    Interpolate the fields of MAPPING at many times and write them to OUTPUT.

    Arguments are as for interpolator-for-wrfchem, plus OUTPUT. The wrfinput only
    provides the grid and the pressure (used for every time), it is not modified.

    --daily-mean (default): each record is the mean of one UTC day, `Times` is that
    day at 00:00:00. --instantaneous: one record per source time. Both write one file.

    --wrf-reference: one file per --interval from --start to --end, for WRF-ReACT's
    chem nudging (auxinput19), e.g. the no3f relaxation reference.
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
    if below_surface_fill is None:
        below_surface_fill = "clamp" if mode == "wrf-reference" else "extend"
    print(f"Below the global model's surface: {below_surface_fill}")
    fields = InterpolatedFields(model, species_map, wrf_ds, below_surface_fill)

    global_attrs = {
        "interpolator_source": f"{global_model} {Path(input_files).resolve()}",
        "interpolator_species_map": str(Path(mapping).resolve()),
        "interpolator_pressure_from": str(Path(wrfinput).resolve()),
        "interpolator_below_surface_fill": below_surface_fill,
    }
    attrs = {
        species_map.aliases_target.get(name, name): {
            "description": f"{name} from species map '{species_map.name}'",
            "interpolator_coeffs": ", ".join(
                f"{k}={v}" for k, v in spec.coeffs.items()
            ),
        }
        for name, spec in species_map.map.items()
    }

    if mode == "wrf-reference":
        if start is None or end is None:
            raise click.ClickException("--wrf-reference needs --start and --end")
        n = write_reference_files(
            fields, wrf, output, start, end, interval, file_template, attrs,
            global_attrs,
        )
        wrf.close()
        print(f"Wrote {n} file(s) to {output}")
        print(
            "WRF namelist: &chem chem_nudge_opt = 1; &time_control "
            f"auxinput19_interval_m = {interval}, io_form_auxinput19 = 2, "
            f"auxinput19_inname = '{file_template.replace('{date}', '<date>').replace('_d01_', '_d<domain>_')}'"
            ". Runs must start at one of the file times."
        )
        return

    start_day = start.date() if start else None
    end_day = end.date() if end else None
    if mode == "daily-mean":
        days = select_days(model.available_times, start_day, end_day, allow_partial)
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
            if (start_day is None or t.date() >= start_day)
            and (end_day is None or t.date() <= end_day)
        ]
        if not times:
            raise click.ClickException("No global-model times in the requested days")
        records = [(t, [t]) for t in times]
        description = "Instantaneous fields at the source times"

    out = FieldsFile(output, wrf, description)
    for k, v in global_attrs.items():
        out.nc.setncattr(k, v)

    for label, source_times in records:
        hours = ", ".join(f"{t:%Y-%m-%d %H}" for t in source_times)
        print(f"{label:%Y-%m-%d %H:%M}: {len(source_times)} source time(s) ({hours})")
        record = fields.from_dataset(daily_mean_dataset(model, source_times))
        for name in record:
            attrs.setdefault(name, {})["interpolator_n_source_times"] = np.int32(
                len(source_times)
            )
        out.write(label, record, attrs)

    out.close()
    wrf.close()
    print(f"Wrote {len(records)} record(s) to {output}")
