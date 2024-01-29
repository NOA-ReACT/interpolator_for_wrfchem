import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np
import xarray as xr


class WRF:
    wrfinput_path: Path
    wrfbdy_path: Path
    wrfinput: nc.Dataset
    wrfbdy: nc.Dataset

    wrfinput_time: dt.datetime
    wrfbdy_time: list[dt.datetime]

    size_south_north: int
    size_west_east: int
    size_bottom_top: int

    def __init__(self, wrfinput_path: Path, wrfbdy_path: Path) -> None:
        self.wrfinput_path = wrfinput_path
        self.wrfbdy_path = wrfbdy_path
        self.wrfinput = nc.Dataset(str(self.wrfinput_path), "r+")
        self.wrfbdy = nc.Dataset(str(self.wrfbdy_path), "r+")

        self.wrfinput.set_auto_mask(False)
        self.wrfbdy.set_auto_mask(False)

        self.wrfinput_time = self.get_wrfinput_time()
        self.wrfbdy_time = self.get_wrfbdy_time()

        self.size_south_north = self.wrfinput.dimensions["south_north"].size
        self.size_west_east = self.wrfinput.dimensions["west_east"].size
        self.size_bottom_top = self.wrfinput.dimensions["bottom_top"].size

    def get_wrfinput_time(self) -> dt.datetime:
        """Read the time of the wrfinput file"""
        t = self.wrfinput.variables["Times"][0]
        t = nc.chartostring(t).item()
        return dt.datetime.strptime(t, "%Y-%m-%d_%H:%M:%S")

    def get_wrfbdy_time(self):
        """Read the lateral boundary update times from the wrfbdy file"""
        times = self.wrfbdy.variables["Times"]
        times = [
            dt.datetime.strptime(nc.chartostring(t).item(), "%Y-%m-%d_%H:%M:%S")
            for t in times
        ]
        return times

    def get_dataset(self):
        """Return a xarray.Dataset containing the basic coordinates and the pressure field."""

        xlong = self.wrfinput.variables["XLONG"][0, :, :]
        xlat = self.wrfinput.variables["XLAT"][0, :, :]
        pres = (
            self.wrfinput.variables["P"][0, :, :, :]
            + self.wrfinput.variables["PB"][0, :, :, :]
        ) * 0.01
        level = np.arange(1, self.size_bottom_top + 1)

        return xr.Dataset(
            {
                "pres": (("bottom_top", "south_north", "west_east"), pres),
            },
            coords={
                "XLONG": (("south_north", "west_east"), xlong),
                "XLAT": (("south_north", "west_east"), xlat),
                "ZNU": (("bottom_top",), self.wrfinput.variables["ZNU"][0, :]),
                "level": (("bottom_top",), level),
            },
            attrs={
                "P_TOP": self.wrfinput.variables["P_TOP"][0],
            },
        )

    def close(self) -> None:
        """Close the netCDF files"""
        self.wrfinput.close()
        self.wrfbdy.close()

    def __str__(self) -> str:
        return f"WRF(wrfinput={self.wrfinput_path} [{self.wrfinput_time:%Y-%m-%d:%H:%M:%S}], wrfbdy={self.wrfbdy_path})"
