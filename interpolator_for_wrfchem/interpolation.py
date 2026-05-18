from typing import Literal

import numpy as np
import xarray as xr
from scipy import interpolate


def interpolate_to_wrf(
    wrf: xr.Dataset,
    gm: xr.Dataset,
    below_surface: Literal["clamp", "extend"] = "clamp",
):
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
        - pres (full-level pressure field)
        - pres_hf (half-level pressure field on dim `level_hf`)
        - Any other variable will be interpolated

    Required variables in `wrf`:
        - hoz_eval_x, hoz_eval_y (2D coordinates)
        - pres (full-level pressure field)
        - pres_hf (half-level pressure field on dim `bottom_top_stag`)

    Args:
        below_surface: Policy for WRF half-levels below the deepest global-model
            half-level. ``"clamp"`` is strictly mass-conservative — the deep WRF
            layers receive zero additional mass. ``"extend"`` extrapolates the
            cumulative-mass curve at the slope of the deepest source layer
            (equivalent to holding the bottom mixing ratio constant down to the
            WRF surface). The latter is more physical for well-mixed near-surface
            fields but adds mass that was not in the source column.
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
        if "pres_hf" not in gm.variables or "pres_hf" not in wrf.variables:
            raise RuntimeError(
                "Mass-conservative vertical interpolation requires `pres_hf` "
                "(half-level pressures) in both the global-model and WRF datasets"
            )
        gm_pres = _interpolate_hoz(wrf, gm, "pres")
        gm_pres_hf = _interpolate_hoz(wrf, gm, "pres_hf", level_dim="level_hf")

    # Interpolate each variable
    coord_vars = {gm_x, gm_y, "level", "level_hf", "pres", "pres_hf", "ZNU", "ZNW"}
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
            out[var] = _interpolate_ver(
                wrf.pres,
                wrf.pres_hf,
                gm_pres,
                gm_pres_hf,
                hoz_var,
                below_surface=below_surface,
            )

    return xr.Dataset(out)


def _interpolate_hoz(
    wrf: xr.Dataset, gm: xr.Dataset, var: str, level_dim: str = "level"
) -> xr.DataArray:
    """
    Horizontally interpolate variable `var` from global model `gm` to WRF `wrf`.

    `level_dim` selects which level-like dimension of `var` to iterate over
    (defaults to "level"; use "level_hf" for half-level fields).
    """

    gm_x = gm.attrs.get("hoz_coord_x", "longitude")
    gm_y = gm.attrs.get("hoz_coord_y", "latitude")
    wrf_x = wrf.attrs.get("hoz_eval_x", "XLONG")
    wrf_y = wrf.attrs.get("hoz_eval_y", "XLAT")

    assert (gm[gm_x].diff(gm_x) > 0).all(), f"{gm_x} values must be increasing"
    assert (gm[gm_y].diff(gm_y) > 0).all(), f"{gm_y} values must be increasing"

    n_levels = gm.sizes[level_dim]
    out = np.empty((n_levels, *wrf[wrf_x].shape))
    for l in range(n_levels):
        interp = interpolate.RectBivariateSpline(
            gm[gm_y],
            gm[gm_x],
            gm[var].isel({level_dim: l}),
        )
        out[l, ...] = interp(wrf[wrf_y], wrf[wrf_x], grid=False)

    return xr.DataArray(out, dims=(level_dim, "south_north", "west_east"))


def _interpolate_ver(
    wrf_pres: xr.DataArray,
    wrf_pres_hf: xr.DataArray,
    gm_pres: xr.DataArray,
    gm_pres_hf: xr.DataArray,
    var: xr.DataArray,
    below_surface: Literal["clamp", "extend"] = "clamp",
) -> xr.DataArray:
    """
    Vertically interpolate `var` from the global model to the WRF grid using a
    mass-conservative cumulative-integral remap.

    Per column:
      1. Build per-layer mass on the source as ``m = q · dp`` (the ``1/g`` factor
         is constant and cancels when we divide back out).
      2. Build the cumulative source curve ``M`` at the source half-levels.
      3. Linearly interpolate ``M`` onto the target half-levels in pressure,
         with explicit out-of-range handling rather than ``interp1d`` fill
         values.
      4. Difference to get the per-target-layer mass and divide by the target
         layer thickness to recover the mixing ratio.

    Args:
        wrf_pres: WRF full-level pressure, dims (bottom_top, south_north, west_east).
        wrf_pres_hf: WRF half-level pressure, dims (bottom_top_stag, south_north, west_east).
        gm_pres: Global model full-level pressure on the WRF horizontal grid,
            dims (level, south_north, west_east).
        gm_pres_hf: Global model half-level pressure on the WRF horizontal grid,
            dims (level_hf, south_north, west_east).
        var: Variable to interpolate, dims (level, south_north, west_east).
        below_surface: See ``interpolate_to_wrf``.
    """

    out = xr.zeros_like(wrf_pres)
    var_np = var.values
    gm_pres_hf_np = gm_pres_hf.values
    wrf_pres_hf_np = wrf_pres_hf.values

    n_x = wrf_pres.sizes["west_east"]
    n_y = wrf_pres.sizes["south_north"]
    for x in range(n_x):
        for y in range(n_y):
            q_src = var_np[:, y, x]
            p_hf_src = gm_pres_hf_np[:, y, x]
            p_hf_tgt = wrf_pres_hf_np[:, y, x]

            q_tgt = _remap_column(p_hf_src, q_src, p_hf_tgt, below_surface)
            out[dict(west_east=x, south_north=y)] = np.clip(
                q_tgt, a_min=0, a_max=None
            )
    return out


def _remap_column(
    p_hf_src: np.ndarray,
    q_src: np.ndarray,
    p_hf_tgt: np.ndarray,
    below_surface: Literal["clamp", "extend"] = "clamp",
) -> np.ndarray:
    """
    Mass-conservative cumulative-integral remap of a single column.

    Inputs may be in any monotonic ordering (top-first or surface-first); both
    sides are independently sorted top-first internally and the result is
    written back in the original target ordering.

    Args:
        p_hf_src: Source half-level pressures, length n+1.
        q_src: Source full-level mixing ratios, length n. Element ``k`` is the
            value of the layer between ``p_hf_src[k]`` and ``p_hf_src[k+1]``
            in the input's native ordering.
        p_hf_tgt: Target half-level pressures, length m+1.
        below_surface: ``"clamp"`` (strictly conservative) or ``"extend"``
            (linear extension of the cumulative curve at the slope of the
            deepest source layer).

    Returns:
        Mixing ratios on the m target full levels, in the input's native
        ordering. Values may be slightly negative from floating-point noise;
        the caller is expected to clip.
    """

    # Sort source half-levels top-first (ascending pressure). Each full-level
    # value q_src[k] lives between half-levels k and k+1 in the input ordering;
    # in the sorted order, the interval between sorted positions i and i+1
    # came from original half-levels (src_order[i], src_order[i+1]) — the
    # full-level index that interval represents is ``min`` of the two.
    src_order = np.argsort(p_hf_src)
    p_hf_src_sorted = p_hf_src[src_order]
    full_src_idx = np.minimum(src_order[:-1], src_order[1:])
    q_src_sorted = q_src[full_src_idx]

    dp_src = np.diff(p_hf_src_sorted)
    m_src = q_src_sorted * dp_src
    M_src = np.concatenate(([0.0], np.cumsum(m_src)))
    M_total = M_src[-1]

    tgt_order = np.argsort(p_hf_tgt)
    p_hf_tgt_sorted = p_hf_tgt[tgt_order]

    # Linearly interpolate the cumulative-mass curve. Out-of-range above the
    # source top is clamped to 0 (no mass above the model). Below the source
    # surface, the default clamp gives M = M_total; the "extend" mode
    # overrides those entries with a linear extension at the deepest-layer
    # slope (constant near-surface mixing ratio).
    M_tgt = np.interp(
        p_hf_tgt_sorted, p_hf_src_sorted, M_src, left=0.0, right=M_total
    )
    if below_surface == "extend" and dp_src.size > 0 and dp_src[-1] > 0:
        slope = m_src[-1] / dp_src[-1]
        mask = p_hf_tgt_sorted > p_hf_src_sorted[-1]
        if mask.any():
            M_tgt[mask] = M_total + slope * (
                p_hf_tgt_sorted[mask] - p_hf_src_sorted[-1]
            )

    m_tgt_sorted = np.diff(M_tgt)
    dp_tgt_sorted = np.diff(p_hf_tgt_sorted)
    with np.errstate(divide="ignore", invalid="ignore"):
        q_tgt_sorted = np.where(
            dp_tgt_sorted > 0, m_tgt_sorted / dp_tgt_sorted, 0.0
        )

    full_tgt_idx = np.minimum(tgt_order[:-1], tgt_order[1:])
    q_tgt = np.empty_like(q_tgt_sorted)
    q_tgt[full_tgt_idx] = q_tgt_sorted
    return q_tgt
