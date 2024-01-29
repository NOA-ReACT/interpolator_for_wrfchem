import numpy as np
import xarray as xr
from scipy import interpolate


def interpolate_to_wrf(wrf: xr.Dataset, gm: xr.Dataset):
    """
    Interpolate all variables of a global model `gm` to the coordinates of a WRF `wrf` dataset.

    Required variables in `gm`:
        - longitude
        - latitude
        - level
        - pres (pressure field)
        - Any other variable will be interpolated

    Required variables in `wrf`:
        - XLONG
        - XLAT
        - pres (pressure field, PB + P)
    """

    # Check required variables
    for var in ["longitude", "latitude", "level", "pres"]:
        if var not in gm.variables:
            raise RuntimeError(f"Variable {var} not found in global model dataset")

    for var in ["XLONG", "XLAT", "pres"]:
        if var not in wrf.variables:
            raise RuntimeError(f"Variable {var} not found in WRF dataset")

    # Horizontally interpolate pressure
    gm_pres = _interpolate_hoz(wrf, gm, "pres")

    # Interpolate each variable
    out = {}
    for var in gm.variables:
        if var in ["longitude", "latitude", "level", "pres"]:
            continue

        hoz_var = _interpolate_hoz(wrf, gm, var)
        out[f"{var}_hoz"] = hoz_var
        out[var] = _interpolate_ver(wrf.pres, gm_pres, hoz_var)

    return xr.Dataset(out)


def _interpolate_hoz(wrf: xr.Dataset, gm: xr.Dataset, var: str) -> xr.Dataset:
    """
    Horizontally interpolate variable `var` from global model `gm` to WRF `wrf`
    """

    out = np.empty((gm.sizes["level"], *wrf.XLONG.shape))
    for l in range(gm.sizes["level"]):
        interp = interpolate.RectBivariateSpline(
            gm.latitude,
            gm.longitude,
            gm[var].isel(level=l),
        )
        out[l, ...] = interp(wrf.XLAT, wrf.XLONG, grid=False)

    return xr.DataArray(out, dims=("level", "south_north", "west_east"))


def _interpolate_ver(
    wrf_pres: xr.DataArray, gm_pres: xr.DataArray, var: xr.DataArray
) -> xr.Dataset:
    """
    Vertically interpolate variable `var` from global model `gm` to WRF `wrf`

    Args:
        wrf_pres: WRF pressure field
        gm_pres: Global model pressure field, interpolated to WRF coordinates. Must have dimensions (level, south_north, west_east)
        var: Variable to interpolate, must have dimensions (level, south_north, west_east)
    """

    out = xr.zeros_like(wrf_pres)
    for x in range(wrf_pres.sizes["west_east"]):
        for y in range(wrf_pres.sizes["south_north"]):
            arr = interpolate.interp1d(
                gm_pres.isel(west_east=x, south_north=y),
                var.isel(west_east=x, south_north=y),
                kind="linear",
                fill_value="extrapolate",
                copy=False,
            )(wrf_pres.isel(west_east=x, south_north=y))
            out[dict(west_east=x, south_north=y)] = np.clip(arr, a_min=0, a_max=None)
    return out
