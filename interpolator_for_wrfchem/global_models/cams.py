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
    """

    dir: Path
    """Where the files are located"""

    times: dict[dt.datetime, dict[str, Path]]
    """Dictionary of times and file paths"""

    required_vars: list[str]
    """The list of variables to be read from the files"""

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

    @property
    def available_times(self) -> list[dt.datetime]:
        """Returns the list of available times"""

        return list(self.times.keys())


class CAMS_ModelLevel_Base(CAMS_Base):
    """
    Base class for CAMS data on native model levels (hybrid sigma-pressure coordinates).

    Handles interactions with a directory of CAMS Global Atmospheric Composition data
    on model levels. The given directory should contain a set of directories named
    `YYYY-MM-DD`, each containing a pair of `data_sfc.nc` and `data_mlev.nc` files.

    The 3D pressure field is created using the surface pressure and the model
    definition [1], which is different between the re-analysis (e.g., EAC4) and the
    operational forecasts. Thus, we have two subclasses of this that are almost the
    same, but use different model level definitions.

    [1] https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """

    level_def: pd.DataFrame
    """Model level definitions (L60 or L137)"""

    def _explore_directory(self):
        """Traverses the directory to find CAMS files, reads time and variable info"""

        self.times = {}
        for filepath in self.dir.iterdir():
            if not filepath.is_dir():
                continue

            sfc_lv = filepath / "data_sfc.nc"
            ml_lv = filepath / "data_mlev.nc"

            if not sfc_lv.exists() or not ml_lv.exists():
                continue

            for date in self._get_dates_from_file(ml_lv):
                self.times[date] = {"sfc": sfc_lv, "ml": ml_lv}

        # Sort dictionary by key
        self.times = dict(sorted(self.times.items()))
        if len(self.times) == 0:
            raise RuntimeError(f"No files found in {self.dir}")

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

        # Read surface pressure from sfc file
        with nc.Dataset(self.times[t]["sfc"], "r") as ds:
            ds.set_auto_mask(False)
            time = ds.variables["valid_time"]
            time_dims = time.dimensions
            cftime_seconds = cftime.date2num(
                t, units=time.units, calendar=time.calendar
            )
            time_idx = np.where(
                time == cftime_seconds
            )  # Should contain 2 values, one for period and one for reference_time
            time_idx = {dim: idx for dim, idx in zip(time_dims, time_idx)}

            sp_var = ds["sp"]

            # The time dimensions might not be in the same order between `valid_time` and `sp`, so
            # we make sure we index the correct things
            time_idx_4d = tuple(
                time_idx[dim] if dim in time_idx else slice(None)
                for dim in sp_var.dimensions
            )
            sp = sp_var[*time_idx_4d]

        # Drop one of the two redundant time dimensions, so now we are working on (level, lat, lon)
        # We pretty much name one of the time dimensions as "level"
        sp = np.squeeze(sp, axis=0)

        with nc.Dataset(self.times[t]["ml"], "r") as ds:
            ds.set_auto_mask(False)

            assert ds[self.required_vars[0]].dimensions[-2:] == (
                "latitude",
                "longitude",
            ), (
                "The last two dimensions of the variables should be latitude and longitude in this order"
            )

            # Read variables
            data = {}
            for var in self.required_vars:
                # Drop time dimensions w/ squeeze
                data[var] = (
                    ("level", "latitude", "longitude"),
                    np.squeeze(ds[var][*time_idx_4d][:]),
                )

            # Read coordinates
            data["longitude"] = ((ds["longitude"][:] - 180) % 360) - 180
            data["latitude"] = ds["latitude"][:]
            data["level"] = ds["model_level"][:]

        # Create the 3D pressure field
        n_levels = len(self.level_def.index)
        a = self.level_def["a [Pa]"].values.reshape(n_levels, 1, 1)
        b = self.level_def["b"].values.reshape(n_levels, 1, 1)

        pres_hf = a + b * sp  # psfc  # Half-level pressure
        pres = (pres_hf[1:, :, :] + pres_hf[:-1, :, :]) / 2  # Full-level pressure
        pres = pres / 100  # Convert to hPa
        data["pres"] = (("level", "latitude", "longitude"), pres)

        # Create the xarray Dataset
        ds = xr.Dataset(data)

        # Create dataset
        ds = ds.set_coords(["longitude", "latitude", "level"])

        # Sort latitude and longitude, surface should be the first level
        ds = ds.sortby("latitude")
        ds = ds.sortby("longitude")
        ds = ds.sortby("level", ascending=False)

        return ds


class CAMS_PressureLevel_Base(CAMS_Base):
    """
    Base class for CAMS data already on pressure levels.

    The given directory should contain a set of directories named `YYYY-MM-DD`, each
    containing a `data_plev.nc` file. No surface pressure file is needed since the
    pressure is directly available as the level coordinate (in hPa).
    """

    def _explore_directory(self):
        """Traverses the directory to find CAMS pressure-level files"""

        self.times = {}
        for filepath in self.dir.iterdir():
            if not filepath.is_dir():
                continue

            plev = filepath / "data_plev.nc"

            if not plev.exists():
                continue

            for date in self._get_dates_from_file(plev):
                self.times[date] = {"pl": plev}

        self.times = dict(sorted(self.times.items()))
        if len(self.times) == 0:
            raise RuntimeError(f"No files found in {self.dir}")

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

        with nc.Dataset(self.times[t]["pl"], "r") as ds:
            ds.set_auto_mask(False)

            time = ds.variables["valid_time"]
            time_dims = time.dimensions
            cftime_seconds = cftime.date2num(
                t, units=time.units, calendar=time.calendar
            )
            time_idx = np.where(time == cftime_seconds)
            time_idx = {dim: idx for dim, idx in zip(time_dims, time_idx)}

            first_var = ds[self.required_vars[0]]

            assert first_var.dimensions[-2:] == (
                "latitude",
                "longitude",
            ), (
                "The last two dimensions of the variables should be latitude and longitude in this order"
            )

            time_idx_4d = tuple(
                time_idx[dim] if dim in time_idx else slice(None)
                for dim in first_var.dimensions
            )

            data = {}
            for var in self.required_vars:
                data[var] = (
                    ("level", "latitude", "longitude"),
                    np.squeeze(ds[var][*time_idx_4d][:]),
                )

            data["longitude"] = ((ds["longitude"][:] - 180) % 360) - 180
            data["latitude"] = ds["latitude"][:]
            pressure_levels = ds["pressure_level"][:]
            data["level"] = pressure_levels

        # Broadcast pressure levels (already in hPa) to a 3D field
        n_lat = len(data["latitude"])
        n_lon = len(data["longitude"])
        pres_3d = np.broadcast_to(
            pressure_levels.reshape(-1, 1, 1),
            (len(pressure_levels), n_lat, n_lon),
        ).copy()
        data["pres"] = (("level", "latitude", "longitude"), pres_3d)

        ds = xr.Dataset(data)
        ds = ds.set_coords(["longitude", "latitude", "level"])

        ds = ds.sortby("latitude")
        ds = ds.sortby("longitude")
        ds = ds.sortby("level", ascending=False)

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
