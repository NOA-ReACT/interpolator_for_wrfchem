import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np
import xarray as xr


class WRFInput:
    path: Path
    nc_file: nc.Dataset
    time: dt.datetime

    size_south_north: int
    size_west_east: int
    size_bottom_top: int

    def __init__(self, wrfinput_path: Path) -> None:
        self.path = wrfinput_path

        self.nc_file = nc.Dataset(str(self.path), "r+")
        self.nc_file.set_auto_mask(False)
        self.time = self._get_time()

        self.size_south_north = self.nc_file.dimensions["south_north"].size
        self.size_west_east = self.nc_file.dimensions["west_east"].size
        self.size_bottom_top = self.nc_file.dimensions["bottom_top"].size

    def _get_time(self) -> dt.datetime:
        """Read the time of the wrfinput file"""
        t = self.nc_file.variables["Times"][0]
        t = nc.chartostring(t).item()
        return dt.datetime.strptime(t, "%Y-%m-%d_%H:%M:%S")

    def get_dataset(self):
        """Return a xarray.Dataset containing the basic coordinates and the pressure field."""

        xlong = self.nc_file.variables["XLONG"][0, :, :]
        xlat = self.nc_file.variables["XLAT"][0, :, :]
        pres = (
            self.nc_file.variables["P"][0, :, :, :]
            + self.nc_file.variables["PB"][0, :, :, :]
        ) * 0.01
        level = np.arange(1, self.size_bottom_top + 1)

        return xr.Dataset(
            {
                "pres": (("bottom_top", "south_north", "west_east"), pres),
            },
            coords={
                "XLONG": (("south_north", "west_east"), xlong),
                "XLAT": (("south_north", "west_east"), xlat),
                "ZNU": (("bottom_top",), self.nc_file.variables["ZNU"][0, :]),
                "level": (("bottom_top",), level),
            },
            attrs={
                "P_TOP": self.nc_file.variables["P_TOP"][0],
            },
        )

    def close(self) -> None:
        """Close the netCDF files"""
        self.nc_file.close()

    def __str__(self) -> str:
        return f"WRFInput(wrfinput={self.path} [{self.time:%Y-%m-%d:%H:%M:%S}])"


class WRFBoundary:
    path: Path
    nc_file: nc.Dataset

    times: list[dt.datetime]

    def __init__(self, wrfbdy_path: Path):
        self.path = wrfbdy_path
        self.nc_file = nc.Dataset(str(self.path), "r+")
        self.nc_file.set_auto_mask(False)

        self.times = self._get_times()

    def _get_times(self) -> list[dt.datetime]:
        """Read the lateral boundary update times from the wrfbdy file"""
        times = self.nc_file.variables["Times"]
        times = [
            dt.datetime.strptime(nc.chartostring(t).item(), "%Y-%m-%d_%H:%M:%S")
            for t in times
        ]
        return times

    def close(self) -> None:
        self.nc_file.close()

    def __str__(self) -> str:
        return f"WRFBoundary(wrfbdy={self.path})"
