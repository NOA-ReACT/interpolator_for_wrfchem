import datetime as dt
from pathlib import Path
import cftime
import netCDF4 as nc
import pandas as pd
import xarray as xr
import importlib.resources
import interpolator_for_wrfchem.res as res


class CAMS_EAC4_ML:
    """
    Handles interactions with a directory of EAC4-model level files. The given directory
    should contains a set of directories named `YYYY-MM-DD`, each containing a pair of
    `levtype_ml.nc` and `levtype_sfc.nc` files. This class handles model level data.

    The 3D pressure field is created using the surface pressure and the model definitions
    from [1], while the methodology is described in [2].

    [1] https://confluence.ecmwf.int/display/UDOC/L60+model+level+definitions
    [2] https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """

    dir: Path
    """Where the files are located"""

    times: dict[dt.datetime, dict[str, Path]]
    """Dictionary of times and file paths"""

    required_vars: list[str]
    """The list of variables to be read from the files"""

    level_def: pd.DataFrame
    """L60 model level definitions"""

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
            importlib.resources.files(res) / "cams_eac4_model_levels.csv"
        )

        self._explore_directory()

    def _explore_directory(self):
        """Traverses the directory to find CAMS files, reads time and variable info"""

        self.times = {}
        for filepath in self.dir.iterdir():
            if not filepath.is_dir():
                continue

            sfc_lv = filepath / "levtype_sfc.nc"
            ml_lv = filepath / "levtype_ml.nc"

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
            time_var = ds.variables["time"][:]
            units = ds.variables["time"].units
            calendar = ds.variables["time"].calendar
        return cftime.num2pydate(time_var, units, calendar)

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
            time = ds.variables["time"]
            time_idx = cftime.date2index(t, time, calendar=time.calendar)
            sp = ds["sp"][time_idx, :, :]  # time, lat, lon

        with nc.Dataset(self.times[t]["ml"], "r") as ds:
            ds.set_auto_mask(False)

            # Read variables
            data = {}
            for var in self.required_vars:
                data[var] = (("level", "latitude", "longitude"), ds[var][time_idx, ...])

            # Read coordinates
            data["longitude"] = ((ds["longitude"][:] - 180) % 360) - 180
            data["latitude"] = ds["latitude"][:]
            data["level"] = ds["level"][:]

        # Create the 3D pressure field
        psfc = sp.reshape(1, *sp.shape)
        a = self.level_def["a [Pa]"].values.reshape(61, 1, 1)
        b = self.level_def["b"].values.reshape(61, 1, 1)

        pres_hf = a + b * psfc  # Half-level pressure
        pres = (pres_hf[1:, :, :] + pres_hf[:-1, :, :]) / 2  # Full-level pressure
        pres = pres / 100  # Convert to hPa
        data["pres"] = (("level", "latitude", "longitude"), pres)

        # Create dataset
        ds = xr.Dataset(data)
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
