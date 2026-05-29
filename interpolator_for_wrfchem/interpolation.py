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
                hoz_var.assign_coords(ZNU=("level", gm["ZNU"].values))
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
    gm_pres_hf: xr.DataArray,
    var: xr.DataArray,
    below_surface: Literal["clamp", "extend"] = "clamp",
) -> xr.DataArray:
    """
    Vertically interpolate `var` from the global model to the WRF grid using a
    mass-conservative cumulative-integral remap.

    Per column the algorithm is:
      1. Build per-layer mass on the source as ``m = q · dp`` (the ``1/g`` factor
         is constant and cancels when we divide back out).
      2. Build the cumulative source curve ``M`` at the source half-levels.
      3. Linearly interpolate ``M`` onto the target half-levels in pressure,
         with explicit out-of-range handling.
      4. Difference to get the per-target-layer mass and divide by the target
         layer thickness to recover the mixing ratio.

    All columns are processed at once by :func:`_remap_columns`, which relies on
    source and target half-level pressures being consistently ordered
    (surface-first, monotonically non-increasing) across the whole grid.

    Args:
        wrf_pres: WRF full-level pressure, dims (bottom_top, south_north, west_east).
        wrf_pres_hf: WRF half-level pressure, dims (bottom_top_stag, south_north, west_east).
        gm_pres_hf: Global model half-level pressure on the WRF horizontal grid,
            dims (level_hf, south_north, west_east).
        var: Variable to interpolate, dims (level, south_north, west_east).
        below_surface: See ``interpolate_to_wrf``.
    """

    q_tgt = _remap_columns(
        gm_pres_hf.values,
        var.values,
        wrf_pres_hf.values,
        below_surface=below_surface,
    )
    return wrf_pres.copy(data=np.clip(q_tgt, a_min=0, a_max=None))


def _remap_columns(
    p_hf_src: np.ndarray,
    q_src: np.ndarray,
    p_hf_tgt: np.ndarray,
    below_surface: Literal["clamp", "extend"] = "clamp",
) -> np.ndarray:
    """
    Mass-conservative cumulative-integral remap of every column at once.

    All arrays carry the level-like axis first: ``(level, south_north,
    west_east)``. Both source and target half-level pressures must be
    monotonically non-increasing (surface-first) along axis 0, *consistently*
    for every column — which lets the per-column ``argsort`` of the scalar
    implementation collapse into a single global flip. This holds for CAMS
    hybrid/pressure levels and WRF eta levels (see module callers); it is
    asserted below so a violation fails loudly rather than silently corrupting
    output.

    Args:
        p_hf_src: Source half-level pressures, shape (n+1, ny, nx).
        q_src: Source full-level mixing ratios, shape (n, ny, nx). Element ``k``
            is the value of the layer between ``p_hf_src[k]`` and
            ``p_hf_src[k+1]``.
        p_hf_tgt: Target half-level pressures, shape (m+1, ny, nx).
        below_surface: ``"clamp"`` (strictly conservative) or ``"extend"``
            (linear extension of the cumulative curve at the slope of the
            deepest source layer).

    Returns:
        Mixing ratios on the m target full levels, shape (m, ny, nx), in the
        input's (surface-first) ordering. Values may be slightly negative from
        floating-point noise; the caller is expected to clip.
    """

    if not (np.diff(p_hf_src, axis=0) <= 0).all():
        raise RuntimeError(
            "Source half-level pressures must be monotonically non-increasing "
            "(surface-first) along the level axis for every column"
        )
    if not (np.diff(p_hf_tgt, axis=0) <= 0).all():
        raise RuntimeError(
            "Target half-level pressures must be monotonically non-increasing "
            "(surface-first) along the level axis for every column"
        )

    # Flip to top-first (ascending pressure). After the flip each full-level
    # value q[k] lies between half-levels p_src[k] and p_src[k+1] in every column.
    p_src = p_hf_src[::-1]  # (L+1, ny, nx)
    q = q_src[::-1]  # (L, ny, nx)
    p_tgt = p_hf_tgt[::-1]  # (M+1, ny, nx)
    L = q.shape[0]

    # Cumulative source mass curve at the source half-levels.
    dp_src = np.diff(p_src, axis=0)
    m_src = q * dp_src
    M_src = np.concatenate(
        [np.zeros((1, *m_src.shape[1:])), np.cumsum(m_src, axis=0)], axis=0
    )  # (L+1, ny, nx)
    M_total = M_src[-1]  # (ny, nx)

    # Batched linear interpolation of M_src(p_src) at the target half-levels.
    # The bracket index is the number of source half-levels strictly below each
    # target value; the loop is over the small (vertical) source axis only.
    idx = np.zeros(p_tgt.shape, dtype=np.intp)
    for k in range(p_src.shape[0]):
        idx += p_src[k][None] < p_tgt
    left = np.clip(idx - 1, 0, L - 1)  # segment left index, (M+1, ny, nx)

    p0 = np.take_along_axis(p_src, left, axis=0)
    p1 = np.take_along_axis(p_src, left + 1, axis=0)
    M0 = np.take_along_axis(M_src, left, axis=0)
    M1 = np.take_along_axis(M_src, left + 1, axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (p_tgt - p0) / (p1 - p0)
    frac = np.where(p1 > p0, frac, 0.0)
    M_tgt = M0 + frac * (M1 - M0)

    # Out-of-range handling, matching np.interp(left=0, right=M_total): above the
    # source top there is no mass; below the source surface the clamp holds M at
    # the column total.
    above_top = p_tgt < p_src[0][None]
    below_surf = p_tgt > p_src[-1][None]
    M_tgt = np.where(above_top, 0.0, M_tgt)
    M_tgt = np.where(below_surf, M_total[None], M_tgt)

    if below_surface == "extend":
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = np.where(dp_src[-1] > 0, m_src[-1] / dp_src[-1], 0.0)  # (ny, nx)
        ext = M_total[None] + slope[None] * (p_tgt - p_src[-1][None])
        M_tgt = np.where(below_surf, ext, M_tgt)

    # Difference back to per-target-layer mass and recover the mixing ratio.
    m_tgt = np.diff(M_tgt, axis=0)
    dp_tgt = np.diff(p_tgt, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        q_tgt = np.where(dp_tgt > 0, m_tgt / dp_tgt, 0.0)

    # Flip back to surface-first to match the WRF ordering.
    return q_tgt[::-1]
