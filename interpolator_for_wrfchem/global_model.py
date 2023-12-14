import datetime as dt
from pathlib import Path
import netCDF4 as nc
import cftime
import numpy as np
from scipy import interpolate


class CAMS_EAC4:
    """
    Handles interactions with a directory of EAC4 files, named `cams_YYYY-MM-DD.nc`, each
    file containing one or more time steps.
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
                out[l, :, :] = interp(ty, tx, grid=False)  # Watch order of tx, ty!

        return out

    def __str__(self) -> str:
        return f"CAMS_EAC4({self.dir}), first date: {min(self.times)}, last date: {max(self.times)}"
