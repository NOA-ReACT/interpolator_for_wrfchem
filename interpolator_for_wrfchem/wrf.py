import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np
import xarray as xr


def hybrid_pressure(c3f, c4f, mu_total, p_top: float) -> np.ndarray:
    """
    WRF's dry hydrostatic pressure on the interfaces (half levels), in hPa.

    Args:
        c3f, c4f: Hybrid coordinate coefficients on the interfaces, shape (bottom_top_stag,)
        mu_total: Total dry column mass (MU + MUB) in Pa, any horizontal shape
        p_top: Model top pressure in Pa
    """

    c3f = np.asarray(c3f)[(...,) + (np.newaxis,) * np.ndim(mu_total)]
    c4f = np.asarray(c4f)[(...,) + (np.newaxis,) * np.ndim(mu_total)]
    return (c3f * np.asarray(mu_total)[np.newaxis] + c4f + p_top) * 0.01


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
            p_hf[k] = C3F[k] · (MU + MUB) + C4F[k] + P_TOP
        and is needed by the mass-conservative vertical interpolation. C3F/C4F are the
        hybrid vertical coordinate coefficients; with the terrain-following coordinate
        (hybrid_opt = 0) they are ZNW and 0.

        The base state column mass (MUB) and the coefficients are included too, so the
        pressure can be recomputed at other times (see `boundary_pressure`).
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
        if "C3F" in self.nc_file.variables:
            c3f = self.nc_file.variables["C3F"][0, :]
            c4f = self.nc_file.variables["C4F"][0, :]
        else:
            # WRF before v3.9, terrain-following coordinate only
            c3f, c4f = znw, np.zeros_like(znw)
        pres_hf = hybrid_pressure(c3f, c4f, mu + mub, p_top)

        return xr.Dataset(
            {
                "pres": (("bottom_top", "south_north", "west_east"), pres),
                "pres_hf": (
                    ("bottom_top_stag", "south_north", "west_east"),
                    pres_hf,
                ),
                "MUB": (("south_north", "west_east"), mub),
            },
            coords={
                "XLONG": (("south_north", "west_east"), xlong),
                "XLAT": (("south_north", "west_east"), xlat),
                "ZNU": (("bottom_top",), self.nc_file.variables["ZNU"][0, :]),
                "ZNW": (("bottom_top_stag",), znw),
                "C3F": (("bottom_top_stag",), c3f),
                "C4F": (("bottom_top_stag",), c4f),
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

    def boundary_mu(self, t_idx: int, bdy: str) -> np.ndarray:
        """
        Perturbation dry column mass (MU) on the outermost row of boundary `bdy` (one of
        BXS, BXE, BYS, BYE) at `self.times[t_idx]`.

        The last time is not stored as a record, it's the end of the last record's
        interval. MU there is reconstructed from the last record's value and tendency,
        which is exactly how WRF gets it.
        """

        tend_name = "MU_BT" + bdy[1:]
        n_records = len(self.times) - 1
        if t_idx < n_records:
            return self.nc_file.variables[f"MU_{bdy}"][t_idx, 0, :]

        dt_last = (self.times[-1] - self.times[-2]).total_seconds()
        return (
            self.nc_file.variables[f"MU_{bdy}"][-1, 0, :]
            + self.nc_file.variables[tend_name][-1, 0, :] * dt_last
        )

    def close(self) -> None:
        self.nc_file.close()

    def __str__(self) -> str:
        return f"WRFBoundary(wrfbdy={self.path})"
