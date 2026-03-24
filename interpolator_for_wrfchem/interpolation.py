import numpy as np
import xarray as xr
from scipy import interpolate


def interpolate_to_wrf(wrf: xr.Dataset, gm: xr.Dataset):
    """
    Interpolate all variables of a global model `gm` to the coordinates of a WRF `wrf` dataset.

    Coordinate names are read from dataset attributes, with defaults for backward compatibility:
        - gm.attrs["hoz_coord_x"] (default: "longitude") — 1D regular x-axis in gm
        - gm.attrs["hoz_coord_y"] (default: "latitude")  — 1D regular y-axis in gm
        - wrf.attrs["hoz_eval_x"] (default: "XLONG") — 2D evaluation x-coords in wrf
        - wrf.attrs["hoz_eval_y"] (default: "XLAT")  — 2D evaluation y-coords in wrf

    Required variables in `gm`:
        - hoz_coord_x, hoz_coord_y (1D coordinates)
        - level
        - pres (pressure field)
        - Any other variable will be interpolated

    Required variables in `wrf`:
        - hoz_eval_x, hoz_eval_y (2D coordinates)
        - pres (pressure field, PB + P)
    """

    gm_x = gm.attrs.get("hoz_coord_x", "longitude")
    gm_y = gm.attrs.get("hoz_coord_y", "latitude")
    wrf_x = wrf.attrs.get("hoz_eval_x", "XLONG")
    wrf_y = wrf.attrs.get("hoz_eval_y", "XLAT")

    # Check required variables
    for var in [gm_x, gm_y, "level", "pres"]:
        if var not in gm.variables:
            raise RuntimeError(f"Variable {var} not found in global model dataset")

    for var in [wrf_x, wrf_y, "pres"]:
        if var not in wrf.variables:
            raise RuntimeError(f"Variable {var} not found in WRF dataset")

    skip_vertical = gm.attrs.get("skip_vertical", False)

    # Horizontally interpolate pressure (only needed for vertical interpolation)
    if not skip_vertical:
        gm_pres = _interpolate_hoz(wrf, gm, "pres")

    # Interpolate each variable
    coord_vars = {gm_x, gm_y, "level", "pres", "ZNU"}
    out = {}
    for var in gm.variables:
        if var in coord_vars:
            continue

        hoz_var = _interpolate_hoz(wrf, gm, var)

        if skip_vertical:
            vert_dim = wrf.pres.dims[0]
            # Sort surface-first (ZNU ≈ 1.0 at index 0) to match the target WRF's
            # bottom_top ordering, then rename the level dim.
            da = (
                hoz_var
                .assign_coords(ZNU=("level", gm["ZNU"].values))
                .sortby("ZNU", ascending=False)
                .drop_vars("ZNU")
                .rename({"level": vert_dim})
            )
            out[f"{var}_hoz"] = da
            out[var] = da
        else:
            out[f"{var}_hoz"] = hoz_var
            out[var] = _interpolate_ver(wrf.pres, gm_pres, hoz_var)

    return xr.Dataset(out)


def _interpolate_hoz(wrf: xr.Dataset, gm: xr.Dataset, var: str) -> xr.Dataset:
    """
    Horizontally interpolate variable `var` from global model `gm` to WRF `wrf`
    """

    gm_x = gm.attrs.get("hoz_coord_x", "longitude")
    gm_y = gm.attrs.get("hoz_coord_y", "latitude")
    wrf_x = wrf.attrs.get("hoz_eval_x", "XLONG")
    wrf_y = wrf.attrs.get("hoz_eval_y", "XLAT")

    assert (gm[gm_x].diff(gm_x) > 0).all(), f"{gm_x} values must be increasing"
    assert (gm[gm_y].diff(gm_y) > 0).all(), f"{gm_y} values must be increasing"

    out = np.empty((gm.sizes["level"], *wrf[wrf_x].shape))
    for l in range(gm.sizes["level"]):
        interp = interpolate.RectBivariateSpline(
            gm[gm_y],
            gm[gm_x],
            gm[var].isel(level=l),
        )
        out[l, ...] = interp(wrf[wrf_y], wrf[wrf_x], grid=False)

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
            var_at_xy = var.isel(west_east=x, south_north=y)

            arr = interpolate.interp1d(
                gm_pres.isel(west_east=x, south_north=y),
                var_at_xy,
                kind="linear",
                fill_value=(0, var_at_xy.isel(level=0).item()),
                copy=False,
                bounds_error=False,
            )(wrf_pres.isel(west_east=x, south_north=y))
            out[dict(west_east=x, south_north=y)] = np.clip(arr, a_min=0, a_max=None)
    return out
