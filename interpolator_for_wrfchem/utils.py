from typing import Literal

import xarray as xr


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
