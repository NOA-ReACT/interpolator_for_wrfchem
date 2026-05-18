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

    def __init__(self, wrfinput_path: Path, read_only=False) -> None:
        self.path = wrfinput_path

        self.nc_file = nc.Dataset(str(self.path), "r+" if not read_only else "r")
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
        """Return a xarray.Dataset containing the basic coordinates and the pressure field.

        Both full-level (mass-level) and half-level (interface) pressures are
        provided. The latter is reconstructed from the WRF dry-mass coordinate:
            p_hf[k] = ZNW[k] · (MU + MUB) + P_TOP
        and is needed by the mass-conservative vertical interpolation.
        """

        xlong = self.nc_file.variables["XLONG"][0, :, :]
        xlat = self.nc_file.variables["XLAT"][0, :, :]
        pres = (
            self.nc_file.variables["P"][0, :, :, :]
            + self.nc_file.variables["PB"][0, :, :, :]
        ) * 0.01
        level = np.arange(1, self.size_bottom_top + 1)
        level_hf = np.arange(0, self.size_bottom_top + 1)

        znw = self.nc_file.variables["ZNW"][0, :]
        mu = self.nc_file.variables["MU"][0, :, :]
        mub = self.nc_file.variables["MUB"][0, :, :]
        p_top = float(self.nc_file.variables["P_TOP"][0])
        mu_total = (mu + mub)[np.newaxis, :, :]
        pres_hf = (znw[:, np.newaxis, np.newaxis] * mu_total + p_top) * 0.01

        return xr.Dataset(
            {
                "pres": (("bottom_top", "south_north", "west_east"), pres),
                "pres_hf": (
                    ("bottom_top_stag", "south_north", "west_east"),
                    pres_hf,
                ),
            },
            coords={
                "XLONG": (("south_north", "west_east"), xlong),
                "XLAT": (("south_north", "west_east"), xlat),
                "ZNU": (("bottom_top",), self.nc_file.variables["ZNU"][0, :]),
                "ZNW": (("bottom_top_stag",), znw),
                "level": (("bottom_top",), level),
                "level_hf": (("bottom_top_stag",), level_hf),
            },
            attrs={
                "P_TOP": p_top,
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

        # The required timesteps are everything mentioned in `md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_`
        # and the last timestep mentioned in `md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_`.
        # The latter is needed only to compute the last tendency, it is not stored as a dimension of `Time`.
        times = np.concatenate(
            [
                self.nc_file.variables[
                    "md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_"
                ][:],
                self.nc_file.variables[
                    "md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_"
                ][-1:],
            ]
        )

        times = [
            dt.datetime.strptime(nc.chartostring(t).item(), "%Y-%m-%d_%H:%M:%S")
            for t in times
        ]
        return times

    def close(self) -> None:
        self.nc_file.close()

    def __str__(self) -> str:
        return f"WRFBoundary(wrfbdy={self.path})"
