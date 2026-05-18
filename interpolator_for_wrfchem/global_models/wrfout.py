import datetime as dt
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pyproj
import xarray as xr

from interpolator_for_wrfchem.global_models.prototype import GlobalModel


class WRFOutput(GlobalModel):
    """
    Use WRF forecast output (wrfout files) as a global model source.

    The source WRF data is on a projected grid (e.g., Lambert Conformal Conic),
    not a regular lat-lon grid. To use the existing RectBivariateSpline-based
    interpolation, we project both source and target coordinates into the source
    WRF's map projection, where the source grid becomes regular (evenly spaced
    at DX/DY intervals).

    The given directory should contain wrfout files for a single domain
    (e.g., wrfout_d01_2021-08-01_00:00:00), one timestep per file.
    """

    dir: Path
    times: dict[dt.datetime, Path]
    required_vars: list[str]

    crs: pyproj.CRS
    transformer: pyproj.Transformer
    dx: float
    dy: float

    def __init__(self, dir: Path | str, required_vars: list[str] = []):
        if isinstance(dir, str):
            dir = Path(dir)
        self.dir = dir
        self.required_vars = required_vars

        self._explore_directory()
        self._init_projection()

    def _explore_directory(self):
        """Scan directory for wrfout files and build time→path mapping."""

        self.times = {}
        for filepath in self.dir.iterdir():
            if not filepath.is_file():
                continue
            if not filepath.name.startswith("wrfout_d"):
                continue

            with nc.Dataset(str(filepath), "r") as ds:
                ds.set_auto_mask(False)
                t = nc.chartostring(ds.variables["Times"][0]).item()
                t = dt.datetime.strptime(t, "%Y-%m-%d_%H:%M:%S")

            self.times[t] = filepath

        self.times = dict(sorted(self.times.items()))
        if len(self.times) == 0:
            raise RuntimeError(f"No wrfout files found in {self.dir}")

    def _init_projection(self):
        """Read projection parameters from the first wrfout file and create pyproj objects."""

        first_file = next(iter(self.times.values()))
        with nc.Dataset(str(first_file), "r") as ds:
            self.dx = float(ds.DX)
            self.dy = float(ds.DY)
            self.znu = ds.variables["ZNU"][0, :]
            self.znw = ds.variables["ZNW"][0, :]
            map_proj = int(ds.MAP_PROJ)

            if map_proj == 1:  # Lambert Conformal Conic
                self.crs = pyproj.CRS.from_proj4(
                    f"+proj=lcc "
                    f"+lat_1={ds.TRUELAT1} +lat_2={ds.TRUELAT2} "
                    f"+lat_0={ds.MOAD_CEN_LAT} +lon_0={ds.STAND_LON} "
                    f"+x_0=0 +y_0=0 +datum=WGS84 +units=m"
                )
            elif map_proj == 3:  # Mercator
                self.crs = pyproj.CRS.from_proj4(
                    f"+proj=merc "
                    f"+lat_ts={ds.TRUELAT1} "
                    f"+lon_0={ds.STAND_LON} "
                    f"+x_0=0 +y_0=0 +datum=WGS84 +units=m"
                )
            else:
                raise NotImplementedError(
                    f"MAP_PROJ={map_proj} is not supported. "
                    f"Supported projections: Lambert Conformal (1), Mercator (3)."
                )

        self.transformer = pyproj.Transformer.from_crs(
            "EPSG:4326", self.crs, always_xy=True
        )

    @property
    def available_times(self) -> list[dt.datetime]:
        return list(self.times.keys())

    def get_dataset(self, t: dt.datetime) -> xr.Dataset:
        """
        Returns the wrfout dataset at time `t`, projected to the source WRF's CRS.

        The returned dataset has 1D `x` and `y` coordinates (in projection space,
        regularly spaced at DX/DY) and sets attrs so the interpolation code uses
        these instead of lat/lon.
        """

        with nc.Dataset(str(self.times[t]), "r") as ds:
            ds.set_auto_mask(False)

            xlat = ds.variables["XLAT"][0, :, :]
            xlong = ds.variables["XLONG"][0, :, :]

            # Compute pressure in hPa
            pres = (
                ds.variables["P"][0, :, :, :] + ds.variables["PB"][0, :, :, :]
            ) * 0.01

            # Reconstruct half-level (interface) pressure from the dry-mass
            # coordinate: p_hf[k] = ZNW[k] · (MU + MUB) + P_TOP.
            mu = ds.variables["MU"][0, :, :]
            mub = ds.variables["MUB"][0, :, :]
            p_top = float(ds.variables["P_TOP"][0])
            mu_total = (mu + mub)[np.newaxis, :, :]
            pres_hf = (
                self.znw[:, np.newaxis, np.newaxis] * mu_total + p_top
            ) * 0.01

            # Read required chemistry variables
            data = {}
            for var in self.required_vars:
                data[var] = (("level", "y", "x"), ds.variables[var][0, :, :, :])

        # Project source coordinates to the source CRS
        x_2d, y_2d = self.transformer.transform(xlong, xlat)

        # Extract 1D axes — the grid is regular in projection space
        x_1d = x_2d[0, :]   # x varies along west_east (columns)
        y_1d = y_2d[:, 0]    # y varies along south_north (rows)

        n_levels = pres.shape[0]
        data["pres"] = (("level", "y", "x"), pres)
        data["pres_hf"] = (("level_hf", "y", "x"), pres_hf)
        data["x"] = x_1d
        data["y"] = y_1d
        data["level"] = np.arange(1, n_levels + 1)
        data["level_hf"] = np.arange(0, n_levels + 1)
        data["ZNU"] = ("level", self.znu)
        data["ZNW"] = ("level_hf", self.znw)

        ds = xr.Dataset(data)
        ds = ds.set_coords(["x", "y", "level", "level_hf", "ZNU", "ZNW"])

        # Sort: surface first (descending level = descending pressure)
        ds = ds.sortby("y")
        ds = ds.sortby("x")
        ds = ds.sortby("level", ascending=False)
        ds = ds.sortby("level_hf", ascending=False)

        # Set coordinate name attrs for the interpolation code
        ds.attrs["hoz_coord_x"] = "x"
        ds.attrs["hoz_coord_y"] = "y"

        return ds

    def transform_target_coords(self, wrf_ds: xr.Dataset) -> xr.Dataset:
        """
        Project the target WRF's XLAT/XLONG into this WRF's projection space.

        Returns a copy of wrf_ds with added x_target/y_target coordinates and
        attrs telling the interpolation code to use them.
        """

        x_target, y_target = self.transformer.transform(
            wrf_ds.XLONG.values, wrf_ds.XLAT.values
        )

        wrf_ds = wrf_ds.assign_coords(
            {
                "x_target": (("south_north", "west_east"), x_target),
                "y_target": (("south_north", "west_east"), y_target),
            }
        )
        wrf_ds.attrs["hoz_eval_x"] = "x_target"
        wrf_ds.attrs["hoz_eval_y"] = "y_target"

        return wrf_ds

    def __str__(self) -> str:
        return f"WRF Output({self.dir})"
