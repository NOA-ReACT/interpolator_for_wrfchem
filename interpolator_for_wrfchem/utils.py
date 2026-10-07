from typing import Literal

import numpy as np
import xarray as xr

from interpolator_for_wrfchem.wrf import hybrid_pressure


def get_boundary_profile(
    var: xr.DataArray | xr.Dataset, boundary: Literal["BXS", "BXE", "BYS", "BYE"]
) -> xr.DataArray | xr.Dataset:
    """
    Extract a boundary profile from a 3D variable.

    The boundaries are named:
        - BXS: Start of domain on the X axis, left edge (west)
        - BXE: End of domain on the X axis, right edge (east)
        - BYS: Start of domain on the Y axis, bottom edge (south)
        - BYE: End of domain on the Y axis, top edge (north)

    Args:
        var: The variable to extract the boundary from
        boundary: The boundary to extract.
    """

    dims = ["bottom_top", "south_north", "west_east"]
    for d in dims:
        if d not in var.dims:
            raise ValueError(f"Variable {var.name} does not have dimension {d}")

    if boundary == "BXS":
        return var.isel(west_east=slice(0, 1))
    elif boundary == "BXE":
        return var.isel(west_east=slice(-1, None))
    elif boundary == "BYS":
        return var.isel(south_north=slice(0, 1))
    elif boundary == "BYE":
        return var.isel(south_north=slice(-1, None))
    else:
        raise ValueError(f"Unknown boundary {boundary}")


def boundary_pressure_hf(wrf_bdy: xr.Dataset, mu: np.ndarray) -> xr.DataArray:
    """
    Interface pressure (`pres_hf`, hPa) of a boundary profile for a given perturbation
    dry column mass, e.g. the boundary's MU at some time from the wrfbdy file.

    Args:
        wrf_bdy: Boundary profile from `get_boundary_profile`, with MUB, C3F, C4F and the
                 P_TOP attribute (see `WRFInput.get_dataset`)
        mu: MU along the boundary, one value per point of the profile
    """

    mub = wrf_bdy["MUB"]
    mu_total = mub.to_numpy() + np.asarray(mu).reshape(mub.shape)
    pres_hf = hybrid_pressure(
        wrf_bdy["C3F"].to_numpy(),
        wrf_bdy["C4F"].to_numpy(),
        mu_total,
        wrf_bdy.attrs["P_TOP"],
    )
    return wrf_bdy["pres_hf"].copy(data=pres_hf)
