"""
Gridded emission inventories (e.g. CAMS-GLOB-ANT) to WRF-Chem `wrfchemi` files.

The inventories are surface fluxes in kg m-2 s-1 on a regular lat-lon grid, one
variable per sector, with a monthly (or daily) time axis. For each output time this
module:

1. regrids every sector to the WRF grid, conserving mass (see `ConservativeRegridder`),
2. spreads it over the lowest WRF layers with the sector's vertical profile,
3. sums the sectors and sources of each WRF emission field with the coefficients of an
   emission map (TOML), which also do the unit conversion,
4. writes one `wrfchemi_d01_<date>` file per output time, for `io_style_emissions = 2`.

Version 1 has no diurnal or weekly cycle: within a month (or day, for daily
inventories) the emissions are constant.
"""

import datetime as dt
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import cftime
import click
import netCDF4 as nc
import numpy as np
import pyproj
import scipy.sparse as sp
import tomli

WRF_EARTH_RADIUS = 6370000.0  # m, the sphere WRF's map projections use
WRF_TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"


# --------------------------------------------------------------------------------------
# WRF grid
# --------------------------------------------------------------------------------------


class WRFGrid:
    """Horizontal grid and layer heights of a wrfinput file."""

    def __init__(self, wrfinput: Path):
        with nc.Dataset(str(wrfinput)) as ds:
            ds.set_auto_mask(False)
            self.global_attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
            self.data_model = ds.data_model
            self.xlat = np.asarray(ds["XLAT"][0], dtype=float)
            self.xlong = np.asarray(ds["XLONG"][0], dtype=float)
            self.dx = float(ds.DX)
            self.dy = float(ds.DY)
            map_proj = int(ds.MAP_PROJ)
            truelat1, truelat2 = float(ds.TRUELAT1), float(ds.TRUELAT2)
            stand_lon = float(ds.STAND_LON)
            cen_lat = float(ds.MOAD_CEN_LAT)
            # Height of the layer interfaces above ground, (bottom_top_stag, sn, we)
            z_w = (ds["PH"][0] + ds["PHB"][0]) / 9.81
            self.z_w_agl = np.asarray(z_w - ds["HGT"][0][np.newaxis], dtype=float)

        self.ny, self.nx = self.xlat.shape

        sphere = f"+R={WRF_EARTH_RADIUS} +units=m +no_defs"
        if map_proj == 1:
            proj = (
                f"+proj=lcc +lat_1={truelat1} +lat_2={truelat2} +lat_0={cen_lat} "
                f"+lon_0={stand_lon} {sphere}"
            )
        elif map_proj == 3:
            proj = f"+proj=merc +lat_ts={truelat1} +lon_0={stand_lon} {sphere}"
        else:
            raise NotImplementedError(
                f"MAP_PROJ={map_proj} is not supported, only Lambert conformal (1) "
                "and Mercator (3)"
            )
        self.transformer = pyproj.Transformer.from_crs(
            pyproj.CRS.from_proj4("+proj=longlat " + sphere),
            pyproj.CRS.from_proj4(proj),
            always_xy=True,
        )

        # The mass points are a regular grid in projected coordinates. Take the origin
        # from cell (0, 0) and check the others agree, which also checks the projection.
        x, y = self.transformer.transform(self.xlong, self.xlat)
        self.x0, self.y0 = float(x[0, 0]), float(y[0, 0])
        jj, ii = np.mgrid[0 : self.ny, 0 : self.nx]
        err = max(
            np.abs(x - (self.x0 + ii * self.dx)).max() / self.dx,
            np.abs(y - (self.y0 + jj * self.dy)).max() / self.dy,
        )
        if err > 0.02:
            raise RuntimeError(
                f"The WRF mass points are up to {err:.3f} cells off a regular grid in "
                "the reconstructed projection; the projection parameters do not match"
            )
        self.projection_error = err

    def locate(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Flat WRF cell index (j * nx + i) of each point, -1 outside the domain."""
        x, y = self.transformer.transform(lon, lat)
        i = np.floor((x - self.x0) / self.dx + 0.5).astype(np.int64)
        j = np.floor((y - self.y0) / self.dy + 0.5).astype(np.int64)
        inside = (i >= 0) & (i < self.nx) & (j >= 0) & (j < self.ny)
        return np.where(inside, j * self.nx + i, -1)


# --------------------------------------------------------------------------------------
# Conservative regridding
# --------------------------------------------------------------------------------------


def _regular_edges(centres: np.ndarray, name: str) -> tuple[float, float]:
    """First centre and spacing of a regular axis."""
    centres = np.asarray(centres, dtype=float)
    # Mean spacing; the axes are often float32, so single steps wobble by ~1e-5 deg
    step = (centres[-1] - centres[0]) / (len(centres) - 1)
    if np.abs(np.diff(centres) - step).max() > 1e-3 * abs(step):
        raise RuntimeError(f"The {name} axis of the inventory is not regular")
    return float(centres[0]), float(step)


class ConservativeRegridder:
    """
    Area-weighted regridding of fluxes from a regular lat-lon grid to the WRF grid.

    Each source cell is split into `n_sub` x `n_sub` sub-cells, and each sub-cell is
    assigned whole to the WRF cell its centre falls in. That gives a sparse matrix W
    (WRF cells x source cells) of overlap areas. A flux field f (per m2) becomes
        F_wrf = (W f) / (W 1),
    the source mass rate that lands in each WRF cell over the source area that lands
    there. The total mass rate over the domain is conserved exactly (up to the
    sub-cell resolution at the domain edge), and the WRF cell area comes out of the
    same sampling, so no map factors are needed.
    """

    def __init__(
        self,
        grid: WRFGrid,
        src_lat: np.ndarray,
        src_lon: np.ndarray,
        n_sub: int = 5,
        margin: float = 1.0,
    ):
        lat0, dlat = _regular_edges(src_lat, "latitude")
        lon0, dlon = _regular_edges(src_lon, "longitude")
        if dlat <= 0 or dlon <= 0:
            raise RuntimeError("The inventory axes must be ascending")

        # Source cells around the domain only
        lon_min, lon_max = grid.xlong.min() - margin, grid.xlong.max() + margin
        lat_min, lat_max = grid.xlat.min() - margin, grid.xlat.max() + margin
        if lon_min < -180 or lon_max > 180:
            raise NotImplementedError("Domains across the date line are not supported")
        self.lat_slice = slice(
            max(0, int(math.floor((lat_min - lat0) / dlat))),
            min(len(src_lat), int(math.ceil((lat_max - lat0) / dlat)) + 1),
        )
        self.lon_slice = slice(
            max(0, int(math.floor((lon_min - lon0) / dlon))),
            min(len(src_lon), int(math.ceil((lon_max - lon0) / dlon)) + 1),
        )
        lats = src_lat[self.lat_slice]
        lons = src_lon[self.lon_slice]
        self.shape = (len(lats), len(lons))

        # Sub-cell centres, as offsets from the cell centre
        off = ((np.arange(n_sub) + 0.5) / n_sub - 0.5)
        sub_lon = (lons[:, None] + off[None, :] * dlon).ravel()  # (nlon * n_sub)

        rows, cols, vals = [], [], []
        r = WRF_EARTH_RADIUS
        for row, lat in enumerate(lats):
            sub_lat = lat + off * dlat  # (n_sub,)
            # Exact area of each sub-cell band on the sphere
            band_area = (
                r**2
                * math.radians(dlon / n_sub)
                * (
                    np.sin(np.radians(sub_lat + dlat / n_sub / 2))
                    - np.sin(np.radians(sub_lat - dlat / n_sub / 2))
                )
            )  # (n_sub,)
            plon, plat = np.meshgrid(sub_lon, sub_lat)  # (n_sub, nlon * n_sub)
            area = np.broadcast_to(band_area[:, None], plon.shape)
            idx = grid.locate(plon.ravel(), plat.ravel())
            src = row * len(lons) + np.broadcast_to(
                np.repeat(np.arange(len(lons)), n_sub)[None, :], plon.shape
            ).ravel()
            inside = idx >= 0
            rows.append(idx[inside])
            cols.append(src[inside])
            vals.append(area.ravel()[inside])

        n_wrf = grid.ny * grid.nx
        self.W = sp.coo_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
            shape=(n_wrf, self.shape[0] * self.shape[1]),
        ).tocsr()
        self.wrf_area = np.asarray(self.W.sum(axis=1)).ravel()  # m2
        if (self.wrf_area <= 0).any():
            raise RuntimeError(
                "Some WRF cells received no source area; the inventory does not cover "
                "the domain"
            )
        self.wrf_shape = (grid.ny, grid.nx)
        # Source cell areas, for the mass check
        self.src_area = np.asarray(self.W.sum(axis=0)).ravel()

    def subset(self, field2d: np.ndarray) -> np.ndarray:
        return field2d[self.lat_slice, self.lon_slice]

    def __call__(self, flux_subset: np.ndarray) -> np.ndarray:
        """Regrid a flux (per m2) already cut with `subset`, to (south_north, west_east)."""
        mass = self.W @ np.asarray(flux_subset, dtype=float).ravel()
        return (mass / self.wrf_area).reshape(self.wrf_shape)


# --------------------------------------------------------------------------------------
# Vertical distribution
# --------------------------------------------------------------------------------------


@dataclass
class VerticalProfile:
    edges: list[float]  # m above ground, len(fractions) + 1, starting at 0
    fractions: list[float]

    def __post_init__(self):
        if len(self.edges) != len(self.fractions) + 1:
            raise ValueError("A profile needs one more edge than fractions")
        if self.edges[0] != 0 or any(np.diff(self.edges) <= 0):
            raise ValueError("Profile edges must start at 0 and increase")
        if not math.isclose(sum(self.fractions), 1.0, abs_tol=1e-6):
            raise ValueError(f"Profile fractions sum to {sum(self.fractions)}, not 1")

    @property
    def top(self) -> float:
        """Top of the highest layer that has emissions"""
        last = max(i for i, f in enumerate(self.fractions) if f > 0)
        return self.edges[last + 1]

    def wrf_fractions(self, z_w_agl: np.ndarray, kemit: int) -> np.ndarray:
        """
        Fraction of the emissions in each of the lowest `kemit` WRF layers, per column,
        shape (kemit, sn, we). Each profile layer is spread evenly in height and split
        by overlap with the WRF layers. Anything above WRF layer `kemit` goes into the
        top one, so the columns always sum to 1.
        """
        out = np.zeros((kemit,) + z_w_agl.shape[1:])
        for p, frac in enumerate(self.fractions):
            if frac == 0:
                continue
            lo, hi = self.edges[p], self.edges[p + 1]
            for k in range(kemit):
                zb = z_w_agl[k]
                zt = z_w_agl[k + 1] if k < kemit - 1 else np.inf
                if k == 0:
                    zb = np.full_like(zb, -np.inf)
                overlap = np.clip(np.minimum(zt, hi) - np.maximum(zb, lo), 0, None)
                out[k] += frac * overlap / (hi - lo)
        return out


# --------------------------------------------------------------------------------------
# Emission map
# --------------------------------------------------------------------------------------


@dataclass
class Source:
    name: str
    file: str  # may contain {year}
    sectors: list[str] | str = "all"  # "all" = every (time, lat, lon) variable
    sum_variable: Optional[str] = None  # excluded from "all", used as a check
    year: Optional[int] = None  # use this year's file and records for every date
    # "auto": read the time axis ("months since ..." or a CF time unit).
    # "day_of_year": record n is day n + 1 of the year, whatever the units say (for
    # files with a broken time axis, e.g. CAMS-GLOB-OCE v3.1: 365 daily records
    # labelled "months since"). 29 February uses 28 February in non-leap years.
    time_axis: str = "auto"


@dataclass
class Target:
    name: str
    units: str
    description: str
    coeffs: dict[str, float] = field(default_factory=dict)


class EmissionMap:
    def __init__(self, path: Path):
        self.path = path
        with open(path, "rb") as fin:
            data = tomli.load(fin)

        meta = data.get("meta", {})
        if meta.get("version") != 1:
            raise ValueError(f"{path}: only version 1 emission maps are supported")
        self.name = meta.get("name", str(path))

        self.sources = {
            name: Source(
                name=name,
                file=raw["file"],
                sectors=raw.get("sectors", "all"),
                sum_variable=raw.get("sum_variable"),
                year=raw.get("year"),
                time_axis=raw.get("time_axis", "auto"),
            )
            for name, raw in data.get("sources", {}).items()
        }
        for src in self.sources.values():
            if src.time_axis not in ("auto", "day_of_year"):
                raise ValueError(
                    f"{path}: time_axis of {src.name} must be 'auto' or 'day_of_year'"
                )

        self.targets = {}
        for name, raw in data.get("targets", {}).items():
            coeffs = {k: float(v) for k, v in raw.get("coeffs", {}).items()}
            unknown = set(coeffs) - set(self.sources)
            if unknown:
                raise ValueError(f"{path}: {name} uses unknown sources {unknown}")
            self.targets[name] = Target(
                name=name,
                units=raw["units"],
                description=raw.get("description", name),
                coeffs=coeffs,
            )

        vertical = data.get("vertical", {})
        self.profiles = {
            name: VerticalProfile(edges=raw["edges"], fractions=raw["fractions"])
            for name, raw in vertical.get("profiles", {}).items()
        }
        self.sector_profiles: dict[str, str] = vertical.get("sectors", {})
        self.default_profile: str = vertical.get("default", "surface")
        for prof in [*self.sector_profiles.values(), self.default_profile]:
            if prof not in self.profiles:
                raise ValueError(f"{path}: unknown vertical profile '{prof}'")

    def profile_for(self, sector: str) -> str:
        return self.sector_profiles.get(sector, self.default_profile)


# --------------------------------------------------------------------------------------
# Reading inventories
# --------------------------------------------------------------------------------------


class InventoryReader:
    """Reads sector fluxes of one source at a given date, regridded to WRF."""

    def __init__(self, source: Source, input_dir: Path, grid: WRFGrid, n_sub: int):
        self.source = source
        self.input_dir = input_dir
        self.grid = grid
        self.n_sub = n_sub
        self.regridder: Optional[ConservativeRegridder] = None
        self._cache: dict[tuple, dict[str, np.ndarray]] = {}
        self.warned_unmapped: set[str] = set()

    def path_for(self, year: int) -> Path:
        return self.input_dir / self.source.file.format(year=year)

    def record_for(self, ds: nc.Dataset, t: dt.datetime) -> tuple[int, str]:
        """Index of the record that covers `t`, by month or day, and a label for it."""
        time = ds["time"]
        year = self.source.year or t.year
        if self.source.time_axis == "day_of_year":
            leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
            day = 28 if (t.month, t.day) == (2, 29) and not leap else t.day
            date = dt.date(year, t.month, day)
            idx = date.timetuple().tm_yday - 1
            if idx >= len(time):
                raise RuntimeError(
                    f"{ds.filepath()} has {len(time)} records, no day {idx + 1}"
                )
            return idx, date.isoformat()
        if time.units.startswith("months since"):
            # cftime only allows "months since" for 360-day calendars, so count the
            # whole months from the base date by hand
            base = dt.date.fromisoformat(time.units.split()[2][:10])
            months = [base.month - 1 + int(round(float(v))) for v in time[:]]
            key = [(base.year + m // 12, m % 12 + 1) for m in months]
            want = (year, t.month)
            label = f"{year}-{t.month:02d}"
        elif "month" in time.units:
            raise RuntimeError(f"Unsupported time units '{time.units}'")
        else:
            dates = cftime.num2date(
                time[:], time.units, getattr(time, "calendar", "standard")
            )
            key = [(d.year, d.month, d.day) for d in dates]
            want = (year, t.month, t.day)
            label = f"{year}-{t.month:02d}-{t.day:02d}"
        if want not in key:
            raise RuntimeError(f"No record for {label} in {ds.filepath()}")
        return key.index(want), label

    def sectors(self, ds: nc.Dataset) -> list[str]:
        if self.source.sectors == "all":
            return [
                v
                for v, var in ds.variables.items()
                if var.dimensions == ("time", "lat", "lon")
                and v != self.source.sum_variable
            ]
        return list(self.source.sectors)

    def read(self, t: dt.datetime) -> tuple[str, dict[str, np.ndarray]]:
        """Regridded flux of every sector at `t`, kg m-2 s-1 (or the file's unit)."""
        year = self.source.year or t.year
        path = self.path_for(year)
        with nc.Dataset(str(path)) as ds:
            ds.set_auto_mask(False)
            rec, label = self.record_for(ds, t)
            key = (str(path), rec)
            if key in self._cache:
                return label, self._cache[key]

            if self.regridder is None:
                print(f"  Building the regridding weights for {self.source.name}")
                self.regridder = ConservativeRegridder(
                    self.grid, ds["lat"][:], ds["lon"][:], n_sub=self.n_sub
                )

            out = {}
            for sector in self.sectors(ds):
                flux = self.regridder.subset(ds[sector][rec])
                flux = np.where(np.isfinite(flux), flux, 0.0)
                out[sector] = self.regridder(flux)

            if self.source.sum_variable:
                total = self.regridder(
                    self.regridder.subset(ds[self.source.sum_variable][rec])
                )
                parts = sum(out.values())
                scale = max(float(np.abs(total).max()), 1e-30)
                err = float(np.abs(parts - total).max()) / scale
                if err > 1e-3:
                    print(
                        f"  Warning: the sectors of {self.source.name} differ from "
                        f"'{self.source.sum_variable}' by up to {err:.2%} of its maximum"
                    )

        self._cache = {key: out}  # keep only the current record
        return label, out


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------


def write_wrfchemi(
    path: Path,
    grid: WRFGrid,
    times: list[dt.datetime],
    fields: list[dict[str, np.ndarray]],
    targets: dict[str, Target],
    attrs: dict[str, str],
):
    """One wrfchemi file with `len(times)` frames."""
    kemit = next(iter(fields[0].values())).shape[0]
    with nc.Dataset(str(path), "w", format=grid.data_model) as ds:
        for k, v in grid.global_attrs.items():
            ds.setncattr(k, v)
        # WRF accepts input files whose TITLE contains " V4." or one of a few markers,
        # EMISSIONS among them (share/input_wrf.F), and stops on anything else
        ds.setncattr("TITLE", "EMISSIONS FROM interpolator-for-wrfchem-emissions")
        for k, v in attrs.items():
            ds.setncattr(k, v)
        ds.createDimension("Time", None)
        ds.createDimension("DateStrLen", 19)
        ds.createDimension("west_east", grid.nx)
        ds.createDimension("south_north", grid.ny)
        ds.createDimension("emissions_zdim", kemit)

        times_var = ds.createVariable("Times", "S1", ("Time", "DateStrLen"))
        for n, t in enumerate(times):
            times_var[n] = np.frombuffer(
                t.strftime(WRF_TIME_FORMAT).encode("ascii"), dtype="S1"
            )

        compress = grid.data_model.startswith("NETCDF4")
        for name, target in targets.items():
            var = ds.createVariable(
                name,
                "f4",
                ("Time", "emissions_zdim", "south_north", "west_east"),
                zlib=compress,
            )
            var.FieldType = np.int32(104)
            var.MemoryOrder = "XYZ"
            var.description = target.description
            var.units = target.units
            var.stagger = ""
            var.coordinates = "XLONG XLAT"
            for n, frame in enumerate(fields):
                var[n] = frame[name].astype("f4")


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def _parse_datetime(ctx, param, value: str) -> dt.datetime:
    for fmt in ("%Y-%m-%d_%H", "%Y-%m-%d_%H:%M:%S", "%Y-%m-%dT%H", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(value, fmt)
        except ValueError:
            pass
    raise click.BadParameter(f"expected YYYY-MM-DD_HH, got {value!r}")


@click.command(name="interpolator-for-wrfchem-emissions")
@click.argument(
    "input_dir",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
)
@click.argument(
    "emission_map",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.argument(
    "wrfinput",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
)
@click.argument(
    "output_dir", type=click.Path(file_okay=False, path_type=Path)
)
@click.option("--start", required=True, callback=_parse_datetime, help="YYYY-MM-DD_HH")
@click.option(
    "--end",
    required=True,
    callback=_parse_datetime,
    help="YYYY-MM-DD_HH, inclusive. Include the end of the last forecast.",
)
@click.option(
    "--interval",
    default=60,
    show_default=True,
    help="Minutes between frames, WRF's auxinput5_interval_m.",
)
@click.option(
    "--frames-per-file",
    default=1,
    show_default=True,
    help="WRF's frames_per_auxinput5. WRF opens the file named after the time it "
    "needs a new one, so the files must line up with the start of every run.",
)
@click.option(
    "--kemit",
    type=int,
    default=None,
    help="Number of emission levels (WRF's kemit). Default: the fewest WRF layers that "
    "hold every profile everywhere. Emissions above the last level go into it.",
)
@click.option(
    "--subcells",
    default=5,
    show_default=True,
    help="Sub-cells per inventory cell side for the conservative regridding.",
)
@click.option("--domain", default=1, show_default=True)
def main(
    input_dir: Path,
    emission_map: Path,
    wrfinput: Path,
    output_dir: Path,
    start: dt.datetime,
    end: dt.datetime,
    interval: int,
    frames_per_file: int,
    kemit: Optional[int],
    subcells: int,
    domain: int,
):
    """
    Regrid the emission inventories of EMISSION_MAP, found in INPUT_DIR, to the grid of
    WRFINPUT and write `wrfchemi_d<domain>_<date>` files to OUTPUT_DIR, for
    `io_style_emissions = 2`. The wrfinput is not modified.

    Set in the WRF namelist: auxinput5_interval_m = --interval, frames_per_auxinput5 =
    --frames-per-file, io_form_auxinput5 = 2, kemit = the value printed at the end.
    """
    emap = EmissionMap(emission_map)
    print(f"Emission map: {emap.name} ({emission_map})")
    print(f"Targets: {', '.join(emap.targets)}")

    grid = WRFGrid(wrfinput)
    print(
        f"WRF grid: {grid.nx} x {grid.ny}, dx = {grid.dx / 1000:g} km, projection "
        f"check within {grid.projection_error:.1e} cells"
    )

    # Vertical levels
    tops = {name: prof.top for name, prof in emap.profiles.items()}
    needed_top = max(tops[emap.profile_for(s)] for s in [*emap.sector_profiles, "*"])
    k_needed = int(np.max(np.argmax(grid.z_w_agl >= needed_top - 1e-6, axis=0)))
    if kemit is None:
        kemit = max(1, k_needed)
    elif kemit < k_needed:
        print(
            f"Warning: kemit = {kemit} is below the {k_needed} layers needed for "
            f"profiles up to {needed_top:g} m; the mass above goes into layer {kemit}"
        )
    fractions = {
        name: prof.wrf_fractions(grid.z_w_agl, kemit)
        for name, prof in emap.profiles.items()
    }

    readers = {
        name: InventoryReader(src, input_dir, grid, subcells)
        for name, src in emap.sources.items()
    }

    times = []
    t = start
    while t <= end:
        times.append(t)
        t += dt.timedelta(minutes=interval)
    if not times:
        raise click.ClickException("--end is before --start")

    output_dir.mkdir(parents=True, exist_ok=True)
    attrs = {
        "interpolator_emission_map": str(emission_map.resolve()),
        "interpolator_input_dir": str(input_dir.resolve()),
        "interpolator_temporal_profile": "none: monthly or daily means of the inventory",
        "interpolator_kemit": str(kemit),
    }

    frames, frame_times, labels_seen = [], [], {}
    n_files = 0
    for t in times:
        frame = {
            name: np.zeros((kemit, grid.ny, grid.nx)) for name in emap.targets
        }
        for src_name, reader in readers.items():
            users = {
                tname: tgt.coeffs[src_name]
                for tname, tgt in emap.targets.items()
                if src_name in tgt.coeffs
            }
            if not users:
                continue
            label, sectors = reader.read(t)
            if labels_seen.get(src_name) != label:
                labels_seen[src_name] = label
                print(f"{t:%Y-%m-%d %H:%M}: {src_name} record {label}")
                _report_mass(src_name, reader, sectors)
            for sector, flux in sectors.items():
                prof = emap.profile_for(sector)
                if sector not in emap.sector_profiles and sector not in reader.warned_unmapped:
                    reader.warned_unmapped.add(sector)
                    print(
                        f"  Note: sector '{sector}' of {src_name} has no vertical "
                        f"profile, using '{prof}'"
                    )
                spread = fractions[prof] * flux[np.newaxis]
                for tname, coeff in users.items():
                    frame[tname] += coeff * spread

        frames.append(frame)
        frame_times.append(t)
        if len(frames) == frames_per_file or t == times[-1]:
            path = output_dir / (
                f"wrfchemi_d{domain:02d}_{frame_times[0].strftime(WRF_TIME_FORMAT)}"
            )
            write_wrfchemi(path, grid, frame_times, frames, emap.targets, attrs)
            n_files += 1
            frames, frame_times = [], []

    print(f"Wrote {n_files} file(s) to {output_dir}")
    print(
        "WRF namelist (&chem / &time_control): emiss_opt = 6, io_style_emissions = 2, "
        f"kemit = {kemit}, auxinput5_interval_m = {interval}, "
        f"frames_per_auxinput5 = {frames_per_file}, io_form_auxinput5 = 2, "
        "auxinput5_inname = 'wrfchemi_d<domain>_<date>'"
    )


def _report_mass(name: str, reader: InventoryReader, sectors: dict[str, np.ndarray]):
    """Domain total of the regridded fluxes, to compare with the inventory."""
    rg = reader.regridder
    total = sum(float((f.ravel() * rg.wrf_area).sum()) for f in sectors.values())
    print(f"  {name}: domain total {total * 86400 * 365.25 / 1e9:.4g} Tg/yr-equivalent")
