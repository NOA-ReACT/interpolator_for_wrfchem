import datetime as dt
from pathlib import Path
import netCDF4 as nc
import cftime
import numpy as np
from scipy import interpolate
import pandas as pd
import importlib.resources
import interpolator_for_wrfchem.res as res


class CAMS_EAC4_PL:
    """
    Handles interactions with a directory of EAC4 files, named `cams_YYYY-MM-DD.nc`, each
    file containing one or more time steps. This class handles pressure level data.
    """

    def __init__(self, dir: Path):
        self.dir = dir
        self.times = {}

        self.read_time_info()

    def read_time_info(self):
        """
        Read the time information from the files in the directory, fills the self.times
        dictionary
        """

        self.times = {}
        for filepath in self.dir.glob("cams_*.nc"):
            for date in self.get_dates_from_file(filepath):
                self.times[date] = filepath

        # Sort dictionary by key
        self.times = dict(sorted(self.times.items()))

    def get_dates_from_file(self, filepath: Path):
        """Yields the dates from a single file, as datetime objects"""

        with nc.Dataset(filepath, "r") as ds:
            time_var = ds.variables["time"][:]
            units = ds.variables["time"].units
            calendar = ds.variables["time"].calendar
        return cftime.num2pydate(time_var, units, calendar)

    def get_var(self, t: dt.datetime, name: str):
        """Returns the variable `name` at time `t`"""

        with nc.Dataset(self.times[t], "r") as ds:
            ds.set_auto_mask(False)

            var = ds.variables[name]
            if var.ndim == 4:
                times = ds.variables["time"][:]
                time_index = cftime.date2index(t, times, calendar=var.calendar)
                return var[time_index, :, :, :]
            else:
                return var[:]

    def interpolate_hoz(
        self, t: dt.datetime, name: str, tx: np.ndarray, ty: np.ndarray
    ):
        """Returns the variable `name` at time `t`, horizontally interpolated at points (tx, ty)"""

        assert tx.shape == ty.shape

        with nc.Dataset(self.times[t], "r") as ds:
            ds.set_auto_mask(False)

            var = ds[name]
            ti = cftime.date2index(t, ds["time"], calendar=ds["time"].calendar)

            latitude = ds["latitude"][:]
            latitude_sort = np.argsort(latitude)
            latitude = latitude[latitude_sort]

            longitude = ds["longitude"][:]
            longitude = ((longitude - 180) % 360) - 180
            longitude_sort = np.argsort(longitude)
            longitude = longitude[longitude_sort]

            out = np.empty((len(ds.dimensions["level"]), *tx.shape))

            arr = var[ti, ...]
            arr = arr[:, latitude_sort, :]
            arr = arr[:, :, longitude_sort]

            for l in range(len(ds.dimensions["level"])):
                interp = interpolate.RectBivariateSpline(
                    latitude,
                    longitude,
                    arr[l],
                )
                out[l, :, :] = np.clip(
                    interp(ty, tx, grid=False), a_min=0, a_max=None
                )  # Watch order of tx, ty!

        return out

    def __str__(self) -> str:
        return f"CAMS_EAC4({self.dir}), first date: {min(self.times)}, last date: {max(self.times)}"


class CAMS_EAC4_ML:
    """
    Handles interactions with a directory of EAC4-model level files. The given directory
    should contains a set of directories named `YYYY-MM-DD`, each containing a pair of
    `levtype_ml.nc` and `levtype_sfc.nc` files. This class handles model level data.
    """

    ml_vars = list[str]
    sfc_vars = list[str]

    def __init__(self, dir: Path):
        self.dir = dir
        self.times = {}

        self.read_time_info()
        self.read_var_names()

    def read_time_info(self):
        """
        Read the time information from the files in the directory, fills the self.times
        dictionary
        """

        self.times = {}
        for filepath in self.dir.iterdir():
            if not filepath.is_dir():
                continue

            sfc_lv = filepath / "levtype_sfc.nc"
            ml_lv = filepath / "levtype_ml.nc"

            if not sfc_lv.exists() or not ml_lv.exists():
                continue

            for date in self.get_dates_from_file(ml_lv):
                self.times[date] = {"sfc": sfc_lv, "ml": ml_lv}

        # Sort dictionary by key
        self.times = dict(sorted(self.times.items()))
        if len(self.times) == 0:
            raise RuntimeError(f"No files found in {self.dir}")

    def get_dates_from_file(self, filepath: Path):
        """Yields the dates from a single file, as datetime objects"""

        with nc.Dataset(filepath, "r") as ds:
            time_var = ds.variables["time"][:]
            units = ds.variables["time"].units
            calendar = ds.variables["time"].calendar
        return cftime.num2pydate(time_var, units, calendar)

    def read_var_names(self):
        """
        Reads the variable names from the first file in the directory, to know which
        variables are in the model level files and which are in the surface level files.
        """

        files = list(self.times.values())[0]
        with nc.Dataset(files["sfc"], "r") as ds:
            self.sfc_vars = list(ds.variables.keys())
        with nc.Dataset(files["ml"], "r") as ds:
            self.ml_vars = list(ds.variables.keys())

    def get_filepath(self, t: dt.datetime, name: str):
        """Returns the filepath that contains `name` for the given time `t`"""

        if name in self.sfc_vars:
            return self.times[t]["sfc"]
        elif name in self.ml_vars:
            return self.times[t]["ml"]
        else:
            raise RuntimeError(f"Variable {name} not found in CAMS_EAC4_ML")

    def get_var(self, t: dt.datetime, name: str):
        """Returns the variable `name` at time `t`"""

        with nc.Dataset(self.get_filepath(t, name), "r") as ds:
            ds.set_auto_mask(False)

            var = ds.variables[name]
            if var.ndim > 2:
                times = ds.variables["time"]
                time_index = cftime.date2index(
                    t, times, calendar=ds.variables["time"].calendar
                )
                return var[time_index, ...]
            else:
                return var[:]

    def get_pressure(self, t: dt.datetime):
        """Returns the pressure field at time `t`"""

        psfc = self.get_var(t, "sp")
        psfc = psfc.reshape(*psfc.shape, 1)

        levels = pd.read_csv(
            importlib.resources.files(res) / "cams_eac4_model_levels.csv"
        )
        a = levels["a [Pa]"].values.reshape(1, 1, 61)
        b = levels["b"].values.reshape(1, 1, 61)

        pres_hf = a + b * psfc

        # Full level pressure is the average of the adjacent half level pressures
        return (pres_hf[:, :, :-1] + pres_hf[:, :, 1:]) / 2

    def interpolate_hoz(
        self, t: dt.datetime, name: str, tx: np.ndarray, ty: np.ndarray
    ):
        """Returns the variable `name` at time `t`, horizontally interpolated at points (tx, ty)"""

        assert tx.shape == ty.shape

        with nc.Dataset(self.get_filepath(t, name), "r") as ds:
            ds.set_auto_mask(False)

            var = ds[name]
            ti = cftime.date2index(t, ds["time"], calendar=ds["time"].calendar)

            latitude = ds["latitude"][:]
            latitude_sort = np.argsort(latitude)
            latitude = latitude[latitude_sort]

            longitude = ds["longitude"][:]
            longitude = ((longitude - 180) % 360) - 180
            longitude_sort = np.argsort(longitude)
            longitude = longitude[longitude_sort]

            out = np.empty((len(ds.dimensions["level"]), *tx.shape))

            arr = var[ti, ...]
            arr = arr[:, latitude_sort, :]
            arr = arr[:, :, longitude_sort]

            for l in range(len(ds.dimensions["level"])):
                interp = interpolate.RectBivariateSpline(
                    latitude,
                    longitude,
                    arr[l],
                )
                out[l, :, :] = interp(ty, tx, grid=False)  # Watch order of tx, ty!

        return out

    def __str__(self) -> str:
        return f"CAMS_EAC4_ML({self.dir}), first date: {min(self.times)}, last date: {max(self.times)}"
