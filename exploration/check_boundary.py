"""
# Checking boundary conditions for cycle 20 of airsense_freerun experiment (2024-07-22 12:00 to 2024-07-22 15:00)
"""

# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import ipywidgets

# %%
wrfinput = xr.open_dataset(
    "/home/thgeorgiou/data/wrf_experiments/airsense_freerun/data/initial_boundary/wrfinput_d01_cycle_20"
).isel(Time=0)
wrfbdy = xr.open_dataset(
    "/home/thgeorgiou/data/wrf_experiments/airsense_freerun/data/initial_boundary/wrfbdy_d01_cycle_20"
).isel(Time=0)

wrfbdy
# %%
cams = (
    xr.open_dataset(
        "/home/thgeorgiou/data/Modelling/AIRSENSE/CAMS_FC_compressed/2024-07-22_12/data_mlev.nc",
        decode_timedelta=True,
    )
    .isel(forecast_reference_time=0)
    .set_index(forecast_period="valid_time")
    .rename({"forecast_period": "time"})
)
cams = cams.sel(time="2024-07-22T12:00:00")

# Longitude from (0-360) to (-180, 180)
cams["longitude"] = cams["longitude"].where(
    cams["longitude"] <= 180, cams["longitude"] - 360
)
cams = cams.sortby("longitude")

# Compute DUST_1 for CAMS
cams["DUST_1"] = (
    cams["aermr04"] * 0.963777 + cams["aermr05"] + cams["aermr06"] * 0.016303
)
cams = cams[["DUST_1"]]
cams

# %% Interpolate CAMS to WRF hoz. grid
cams_interp = cams.interp(
    longitude=wrfinput.XLONG, latitude=wrfinput.XLAT, method="linear"
)
cams_interp


# %% Pick a boundary and compare with CAMS
def boundary_to_name_and_sel(boundary):
    if boundary == "west":
        return "BXS", dict(west_east=0)
    elif boundary == "east":
        return "BXE", dict(west_east=-1)
    elif boundary == "south":
        return "BYS", dict(south_north=0)
    elif boundary == "north":
        return "BYE", dict(south_north=-1)
    raise ValueError(f"Unknown boundary: {boundary}")


name, sel = boundary_to_name_and_sel("south")
wrfbdy_var = wrfbdy[f"DUST_1_{name}"].isel(bdy_width=0)
cams_var = cams_interp["DUST_1"].sel(**sel) * 1e9

fig, axes = plt.subplots(2, 1, figsize=(10, 6))
wrfbdy_var.plot(ax=axes[0], cmap="viridis", vmax=2)
cams_var.plot(ax=axes[1], cmap="viridis", vmax=2)

# Flip Y axis for CAMS plot
axes[1].set_ylim(axes[1].get_ylim()[::-1])
# %%
