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
    Specific for the ADS-Beta that launched in Summer of 2024.

    Handles interactions w/ a directory of CAMS Global Atmospheric Composition forecasts.
    The given directory should contain a set of directories named `YYYY-MM-DD`, each
    containing a pair of `levtype_ml.nc` and `levtype_sfc.nc` files. This class handles
    model level data.

    The 3D pressure field is created using the surface pressure and the model definition [1],
    which is different between the re-analysis (e.g., EAC4) and the operational forecasts.
    Thus, we have two subclasses of this that are almost the same, but use different model
    level definitions.

    [1] https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """

    dir: Path
    """Where the files are located"""

    times: dict[dt.datetime, dict[str, Path]]
    """Dictionary of times and file paths"""

    required_vars: list[str]
    """The list of variables to be read from the files"""

    level_def: pd.DataFrame
    """Model level definitions (L60 or L137)"""

    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS EAC4 reanalysis product, in model levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                           returned datasets
        """
        if isinstance(dir, str):
            dir = Path(dir)
        self.dir = dir
        self.required_vars = required_vars

        # Read information about the vertical levels, required for creating the
        # 3D pressure field
        self.level_def = pd.read_csv(
            importlib.resources.files(res)
            / "cams_l60.csv"  # Change filename in subclass!
        )

        self._explore_directory()

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

    def _get_dates_from_file(self, filepath: Path) -> list[dt.datetime]:
        """Returns the date from a single file, as datetime object"""

        with nc.Dataset(filepath, "r") as ds:

            time_var = ds.variables["valid_time"][:]
            dates = time_var.flatten().tolist()
            units = ds.variables["valid_time"].units
            calendar = ds.variables["valid_time"].calendar

        return cftime.num2pydate(dates, units, calendar)

    def get_dataset(self, t: dt.datetime) -> xr.Dataset:
        """
        Returns the CAMS EAC4 dataset at time `t`

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
            cftime_seconds = cftime.date2num(
                t, units=time.units, calendar=time.calendar
            )
            time_idx = np.where(
                time == cftime_seconds
            )  # Should contain 2 values, one for period and one for reference_time
            sp = ds["sp"][*time_idx, :, :]  # period, reference_time, lat, lon

        # Drop one of the two redundant time dimensions, so now we are working on (level, lat, lon)
        # We pretty much name one of the time dimensions as "level"
        sp = np.squeeze(sp, axis=0)

        with nc.Dataset(self.times[t]["ml"], "r") as ds:
            ds.set_auto_mask(False)

            # Read variables
            data = {}
            for var in self.required_vars:

                # Drop time dimensions w/ squeeze
                data[var] = (
                    ("level", "latitude", "longitude"),
                    np.squeeze(ds[var][*time_idx, ...][:]),
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

    @property
    def available_times(self) -> list[dt.datetime]:
        """Returns the list of available times"""

        return list(self.times.keys())


class CAMS_EAC4(CAMS_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS EAC4 reanalysis product, in model levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        if isinstance(dir, str):
            dir = Path(dir)
        self.dir = dir
        self.required_vars = required_vars

        # Read information about the vertical levels, required for creating the
        # 3D pressure field
        self.level_def = pd.read_csv(importlib.resources.files(res) / "cams_l60.csv")

        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS EAC4({self.dir})"


class CAMS_Global_Forecasts(CAMS_Base):
    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        """
        Prepare a CAMS EAC4 reanalysis product, in model levels

        Args:
            dir: The directory containing the files
            required_vars: A list of variables that are used and will be present in the
                            returned datasets
        """
        if isinstance(dir, str):
            dir = Path(dir)
        self.dir = dir
        self.required_vars = required_vars

        # Read information about the vertical levels, required for creating the
        # 3D pressure field
        self.level_def = pd.read_csv(importlib.resources.files(res) / "cams_l137.csv")

        self._explore_directory()

    def __str__(self) -> str:
        return f"CAMS Global Forecasts({self.dir})"
