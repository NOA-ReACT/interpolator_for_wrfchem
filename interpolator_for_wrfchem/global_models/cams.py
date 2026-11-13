import datetime as dt
import importlib.resources
from pathlib import Path

import cftime
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr

import interpolator_for_wrfchem.res as res
from interpolator_for_wrfchem.global_models.prototype import GlobalModel


class CAMS_Base(GlobalModel):
    """
    Common base for all CAMS global model implementations.

    Provides shared attributes and utilities used by both model-level and
    pressure-level subclasses.

    The data directory contains one subdirectory per download (usually one per day,
    named `YYYY-MM-DD`). Inside it, any number of netCDF files can hold the fields:
    they are told apart by their content, not by their name. Files that have the
    subclass' level dimension are multi-level files, and the file that has `sp` is the
    surface file. A large request can therefore be downloaded in several pieces (e.g.
    a few variables each) and dropped into the same directory as is. All multi-level
    files of a directory must share the same grid, levels and times.
    """

    dir: Path
    """Where the files are located"""

    times: dict[dt.datetime, dict]
    """
    For each available time, the files to read from: `sfc` is the surface file and
    `lev` maps each multi-level variable to the file that holds it.
    """

    required_vars: list[str]
    """The list of variables to be read from the files"""

    level_dim: str
    """Name of the vertical dimension that identifies the multi-level files"""

    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        if isinstance(dir, str):
            dir = Path(dir)
        self.dir = dir
        self.required_vars = required_vars

    def _get_dates_from_file(self, filepath: Path) -> list[dt.datetime]:
        """Returns the date from a single file, as datetime object"""

        with nc.Dataset(filepath, "r") as ds:
            time_var = ds.variables["valid_time"][:]
            dates = time_var.flatten().tolist()
            units = ds.variables["valid_time"].units
            calendar = ds.variables["valid_time"].calendar

        return cftime.num2pydate(dates, units, calendar)

    def _explore_directory(self):
        """Traverses the directory to find CAMS files, reads time and variable info"""

        self.times = {}
        for subdir in sorted(self.dir.iterdir()):
            if not subdir.is_dir():
                continue

            sfc = None
            lev: dict[str, Path] = {}
            for filepath in sorted(subdir.glob("*.nc")):
                with nc.Dataset(filepath, "r") as ds:
                    if "valid_time" not in ds.variables:
                        continue
                    if self.level_dim in ds.dimensions:
                        for var in ds.variables:
                            lev.setdefault(var, filepath)
                    elif "sp" in ds.variables and sfc is None:
                        sfc = filepath

            if sfc is None or len(lev) == 0:
                continue

            for filepath in set(lev.values()):
                for date in self._get_dates_from_file(filepath):
                    self.times[date] = {"sfc": sfc, "lev": lev}

        # Sort dictionary by key
        self.times = dict(sorted(self.times.items()))
        if len(self.times) == 0:
            raise RuntimeError(f"No files found in {self.dir}")

    @staticmethod
    def _read_at_time(ds: nc.Dataset, var: str, t: dt.datetime) -> np.ndarray:
        """
        Reads variable `var` at valid time `t`, without its time dimensions.

        CAMS files have two time dimensions (reference time and forecast period),
        with `valid_time` defined on both. Their order is not always the same
        between `valid_time` and the variables, so they are matched by name.
        """

        time = ds.variables["valid_time"]
        time_num = cftime.date2num(t, units=time.units, calendar=time.calendar)
        matches = np.where(time[:] == time_num)
        if len(matches[0]) == 0:
            raise RuntimeError(f"Time {t} not found in {ds.filepath()}")
        time_idx = {dim: int(idx[0]) for dim, idx in zip(time.dimensions, matches)}

        variable = ds.variables[var]
        return np.asarray(
            variable[tuple(time_idx.get(dim, slice(None)) for dim in variable.dimensions)]
        )

    def _read_surface_pressure(self, t: dt.datetime) -> np.ndarray:
        """Returns the surface pressure at time `t` in Pa, dimensions (lat, lon)"""

        with nc.Dataset(self.times[t]["sfc"], "r") as ds:
            ds.set_auto_mask(False)
            return self._read_at_time(ds, "sp", t)

    def _read_level_variables(self, t: dt.datetime) -> dict[str, np.ndarray]:
        """
        Reads all `required_vars` at time `t` from the multi-level files, each with
        dimensions (level, lat, lon). The level, `latitude` and `longitude` axes are
        included in the result, as they are in the files, under the names
        `self.level_dim`, `latitude` and `longitude`.
        """

        lev = self.times[t]["lev"]
        missing = [var for var in self.required_vars if var not in lev]
        if missing:
            raise RuntimeError(
                f"Variables {', '.join(missing)} not found in any {self.level_dim} "
                f"file of {self.times[t]['sfc'].parent}"
            )

        out: dict[str, np.ndarray] = {}
        axes = (self.level_dim, "latitude", "longitude")
        for filepath in sorted({lev[var] for var in self.required_vars}):
            with nc.Dataset(filepath, "r") as ds:
                ds.set_auto_mask(False)

                for axis in axes:
                    values = np.asarray(ds[axis][:])
                    if axis not in out:
                        out[axis] = values
                    elif not np.array_equal(out[axis], values):
                        raise RuntimeError(
                            f"{filepath} has a different {axis} axis than the "
                            f"other files of the same directory"
                        )

                for var in self.required_vars:
                    if lev[var] != filepath:
                        continue
                    assert ds[var].dimensions[-3:] == axes, (
                        f"The last three dimensions of {var} should be "
                        f"{', '.join(axes)} in this order"
                    )
                    out[var] = self._read_at_time(ds, var, t)

        return out

    @property
    def available_times(self) -> list[dt.datetime]:
        """Returns the list of available times"""

        return list(self.times.keys())


class CAMS_ModelLevel_Base(CAMS_Base):
    """
    Base class for CAMS data on native model levels (hybrid sigma-pressure coordinates).

    Handles interactions with a directory of CAMS Global Atmospheric Composition data
    on model levels, see `CAMS_Base` for the directory layout. The usual download is
    one `data_sfc.nc` and one `data_mlev.nc` per `YYYY-MM-DD` directory.

    The 3D pressure field is created using the surface pressure and the model
    definition [1], which is different between the re-analysis (e.g., EAC4) and the
    operational forecasts. Thus, we have two subclasses of this that are almost the
    same, but use different model level definitions.

    The files do not need to hold all model levels. The upper levels can be left out
    of the download, as long as the levels that are there are contiguous and reach
    the surface (e.g. levels 46 to 137 of L137 start at 42 hPa, enough for a WRF
    domain with its top at 50 hPa). WRF layers above the highest downloaded level
    receive no mass.

    [1] https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """

    level_dim = "model_level"

    level_def: pd.DataFrame
    """Model level definitions (L60 or L137)"""

    def get_dataset(self, t: dt.datetime) -> xr.Dataset:
        """
        Returns the CAMS dataset at time `t` for model-level data.

        Coordinates:
            - Longitude
            - Latitude,
            - Level (pressure)

        Vars:
            - All variables in `required_vars`, passed when creating the object
        """

        sp = self._read_surface_pressure(t)  # (lat, lon), in Pa
        fields = self._read_level_variables(t)

        # Order the levels top-first, like the level definition table
        levels = np.rint(fields.pop("model_level")).astype(int)
        order = np.argsort(levels)
        levels = levels[order]

        # `level_def` has n_levels+1 rows: each is a half-level edge, from model
        # top (row 0) to surface (row n_levels). Full level k lies between rows
        # k-1 and k. The files might only hold the lowest levels of the model.
        n_levels = len(self.level_def.index) - 1
        if not (np.diff(levels) == 1).all() or levels[-1] != n_levels:
            raise RuntimeError(
                f"The model levels of {self.times[t]['sfc'].parent} must be contiguous "
                f"and reach the surface (level {n_levels}), found {levels[0]} to "
                f"{levels[-1]} ({len(levels)} levels)"
            )
        first_hf = levels[0] - 1

        data = {}
        for var in self.required_vars:
            data[var] = (("level", "latitude", "longitude"), fields[var][order])

        data["longitude"] = ((fields["longitude"] - 180) % 360) - 180
        data["latitude"] = fields["latitude"]
        data["level"] = levels

        # Create the 3D pressure field. Full-level (centre) pressure is the
        # average of the bracketing half-level pressures.
        a = self.level_def["a [Pa]"].values[first_hf:].reshape(-1, 1, 1)
        b = self.level_def["b"].values[first_hf:].reshape(-1, 1, 1)

        pres_hf = a + b * sp  # Half-level (interface) pressure, Pa
        pres = (pres_hf[1:, :, :] + pres_hf[:-1, :, :]) / 2  # Full-level pressure
        pres = pres / 100  # Convert to hPa
        pres_hf = pres_hf / 100  # Convert to hPa
        data["pres"] = (("level", "latitude", "longitude"), pres)
        data["pres_hf"] = (("level_hf", "latitude", "longitude"), pres_hf)
        data["level_hf"] = np.arange(first_hf, n_levels + 1)

        # Create the xarray Dataset
        ds = xr.Dataset(data)

        # Create dataset
        ds = ds.set_coords(["longitude", "latitude", "level", "level_hf"])

        # Sort latitude and longitude, surface should be the first level
        ds = ds.sortby("latitude")
        ds = ds.sortby("longitude")
        ds = ds.sortby("level", ascending=False)
        ds = ds.sortby("level_hf", ascending=False)

        return ds


class CAMS_PressureLevel_Base(CAMS_Base):
    """
    Base class for CAMS data already on pressure levels.

    See `CAMS_Base` for the directory layout. The usual download is one
    `data_sfc.nc` and one `data_plev.nc` per `YYYY-MM-DD` directory. The surface
    pressure is needed for the bottom edge of the lowest layer; the pressure of the
    other levels is directly available as the level coordinate (in hPa).
    """

    level_dim = "pressure_level"

    def get_dataset(self, t: dt.datetime) -> xr.Dataset:
        """
        Returns the CAMS dataset at time `t` for pressure-level data.

        Coordinates:
            - Longitude
            - Latitude
            - Level (pressure, in hPa)

        Vars:
            - All variables in `required_vars`, passed when creating the object
        """

        # Surface pressure is needed for the bottom half-level edge of the
        # mass-conservative vertical remap
        sp = self._read_surface_pressure(t)  # (lat, lon), in Pa
        fields = self._read_level_variables(t)

        data = {}
        for var in self.required_vars:
            data[var] = (("level", "latitude", "longitude"), fields[var])

        data["longitude"] = ((fields["longitude"] - 180) % 360) - 180
        data["latitude"] = fields["latitude"]
        pressure_levels = fields["pressure_level"]

        # Sort levels top-first (ascending pressure) so the edge construction
        # below is independent of the file's native ordering, and apply the
        # same permutation to all variables.
        sort_idx = np.argsort(pressure_levels)
        pressure_levels = pressure_levels[sort_idx]
        for var in self.required_vars:
            dims, arr = data[var]
            data[var] = (dims, arr[sort_idx, :, :])

        # Drop the topmost full level — pressure-level CAMS has no surface
        # pressure analog at the top, and the highest level (typically <5 hPa)
        # carries negligible aerosol/chemistry mass. Its pressure is reused
        # as the top half-level edge of the kept column.
        p_top_original = float(pressure_levels[0])
        p_kept = pressure_levels[1:]
        for var in self.required_vars:
            dims, arr = data[var]
            data[var] = (dims, arr[1:, :, :])
        data["level"] = p_kept

        n_lat = len(data["latitude"])
        n_lon = len(data["longitude"])
        sp_hpa = sp / 100.0  # Pa -> hPa, shape (lat, lon)

        # Broadcast kept pressure levels (already in hPa) to a 3D full-level field
        pres_3d = np.broadcast_to(
            p_kept.reshape(-1, 1, 1),
            (len(p_kept), n_lat, n_lon),
        ).copy()
        data["pres"] = (("level", "latitude", "longitude"), pres_3d)

        # Build half-level (edge) pressures. For N-1 kept centroids we emit
        # N edges: top edge = dropped level's pressure; interior edges =
        # midpoints between kept centroids; bottom edge = surface pressure.
        n_edges = len(pressure_levels)  # = (N-1) + 1
        pres_hf = np.empty((n_edges, n_lat, n_lon))
        pres_hf[0, :, :] = p_top_original
        if n_edges > 2:
            mids = (p_kept[:-1] + p_kept[1:]) / 2.0
            pres_hf[1:-1, :, :] = mids.reshape(-1, 1, 1)
        pres_hf[-1, :, :] = sp_hpa
        # High-terrain columns: surface can sit above the deepest kept centroid.
        # Clamp so half-levels never exceed sp; the resulting zero-thickness
        # bottom layers are tolerated by _remap_column's dp>0 guard.
        pres_hf = np.minimum(pres_hf, sp_hpa[None, :, :])

        data["pres_hf"] = (("level_hf", "latitude", "longitude"), pres_hf)
        data["level_hf"] = np.arange(n_edges)

        ds = xr.Dataset(data)
        ds = ds.set_coords(["longitude", "latitude", "level", "level_hf"])

        ds = ds.sortby("latitude")
        ds = ds.sortby("longitude")
        ds = ds.sortby("level", ascending=False)
        ds = ds.sortby("level_hf", ascending=False)

        return ds


class CAMS_EAC4(CAMS_ModelLevel_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS EAC4 reanalysis product, in model levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        super().__init__(dir, required_vars)
        self.level_def = pd.read_csv(importlib.resources.files(res) / "cams_l60.csv")
        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS EAC4({self.dir})"


class CAMS_Global_Forecasts(CAMS_ModelLevel_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS Global Forecasts product, in model levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        super().__init__(dir, required_vars)
        self.level_def = pd.read_csv(importlib.resources.files(res) / "cams_l137.csv")
        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS Global Forecasts({self.dir})"


class CAMS_EAC4_Pressure(CAMS_PressureLevel_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS EAC4 reanalysis product, on pressure levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        super().__init__(dir, required_vars)
        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS EAC4 Pressure Levels({self.dir})"


class CAMS_Global_Forecasts_Pressure(CAMS_PressureLevel_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS Global Forecasts product, on pressure levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        super().__init__(dir, required_vars)
        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS Global Forecasts Pressure Levels({self.dir})"
