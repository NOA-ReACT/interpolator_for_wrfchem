import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np
import xarray as xr


class MetEm:
    """
    Handles reading data from a directory of met_em files
    """

    path: Path
    znu: np.ndarray
    p_top: int
    domain_size: tuple[int, int, int]

    def __init__(self, met_em_path: Path, wrf: xr.Dataset):
        """
        Read metadata from a met_em directory

        Args:
            met_em_path: Path to the met_em directory
            wrf: The wrfinput file, required to read the vertical levels. Must contain
                    the ZNU variable and P_TOP as an attribute. It must also have the
                    dimensions "bottom_top", "south_north" and "west_east".
        """

        self.path = met_em_path

        self.times = {}
        self._read_time_info()

        self.znu = wrf["ZNU"][:]
        self.p_top = wrf.attrs["P_TOP"]

        self.domain_size = (
            wrf.sizes["bottom_top"],
            wrf.sizes["south_north"],
            wrf.sizes["west_east"],
        )

    def _read_time_info(self):
        """
        Read the time information from the files in the directory, fills the self.times
        dictionary. The date of each met_em file is read from the filename.
        """

        self.times = {}
        for filepath in self.path.glob("met_em.d0*.nc"):
            name = filepath.stem
            datestr = name.split(".")[2]
            date = dt.datetime.strptime(datestr, "%Y-%m-%d_%H:%M:%S")

            self.times[date] = filepath

        # Sort dictionary by key
        self.times = dict(sorted(self.times.items()))

    def get_var(self, t: dt.datetime, name: str) -> np.ndarray:
        """
        Returns the variable `name` at time `t`.

        The time dimension is removed from the returned array.
        """

        with nc.Dataset(self.times[t], "r") as ds:
            ds.set_auto_mask(False)

            var = ds[name]
            if var.dimensions[0] == "Time":
                return ds[name][0, :]
            else:
                return ds[name][:]

    def get_pres(self, t: dt.datetime) -> xr.DataArray:
        """Returns the pressure field at time `t`"""

        psfc = self.get_var(t, "PSFC")
        pres = np.empty(self.domain_size)
        for l in reversed(range(len(self.znu))):
            pres[l, :, :] = psfc * self.znu[l].item() + self.p_top * (
                1 - self.znu[l].item()
            )

        return xr.DataArray(
            pres * 0.01, dims=("bottom_top", "south_north", "west_east")
        )

    def __str__(self) -> str:
        return f"MetEm(path={self.path})"
