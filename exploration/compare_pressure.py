"""Compare pressure between met_em and wrfinput on the boundaries"""

# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import ipywidgets

# %%
diag = xr.open_dataset("/home/thgeorgiou/work/interpolator_for_wrfchem/diag_pres.nc")
diag

# %%
fig, axes = plt.subplots(3, 1, subplot_kw={"projection": ccrs.PlateCarree()})

level = 1

met_em_pres = diag["met_em_pres"].isel(bottom_top=level)
wrf_pres = diag["wrf_pres"].isel(bottom_top=level)
diff = met_em_pres - wrf_pres

met_em_pres.plot(x="XLONG", y="XLAT", ax=axes[0], cmap="viridis", add_colorbar=True)
wrf_pres.plot(x="XLONG", y="XLAT", ax=axes[1], cmap="viridis", add_colorbar=True)
diff.plot(x="XLONG", y="XLAT", ax=axes[2], cmap="coolwarm", add_colorbar=True)

# %% Compute mean difference in hPa for each model level
total_diff = diag["met_em_pres"] - diag["wrf_pres"]
mean_diff = total_diff.mean(dim=["west_east", "south_north"])
max_diff = total_diff.max(dim=["west_east", "south_north"])
min_diff = total_diff.min(dim=["west_east", "south_north"])
mean_diff.plot(y="bottom_top")
max_diff.plot(y="bottom_top")
min_diff.plot(y="bottom_top")
# Not insignificant, but they mostly relate to where there are mountains and the first level goes up or down.

# %% Check the boundaries only, where we are interested in the pressure difference

fig = plt.figure(figsize=(12, 8))
# Add a 3x3 gridspec, center for the main plot, and the rest for the boundaries
gs = fig.add_gridspec(3, 3)

top_ax = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
top_ax.set_title("Top Boundary")
bot_ax = fig.add_subplot(gs[2, 1], projection=ccrs.PlateCarree())
bot_ax.set_title("Bottom Boundary")
left_ax = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
left_ax.set_title("Left Boundary")
right_ax = fig.add_subplot(gs[1, 2], projection=ccrs.PlateCarree())
right_ax.set_title("Right Boundary")

top_diff = diag["met_em_pres"].isel(south_north=0) - diag["wrf_pres"].isel(
    south_north=0
)
bot_diff = diag["met_em_pres"].isel(south_north=-1) - diag["wrf_pres"].isel(
    south_north=-1
)
left_diff = diag["met_em_pres"].isel(west_east=0) - diag["wrf_pres"].isel(west_east=0)
right_diff = diag["met_em_pres"].isel(west_east=-1) - diag["wrf_pres"].isel(
    west_east=-1
)

top_diff.plot(x="XLONG", y="bottom_top", ax=top_ax, cmap="coolwarm", add_colorbar=True)
bot_diff.plot(x="XLONG", y="bottom_top", ax=bot_ax, cmap="coolwarm", add_colorbar=True)
left_diff.plot(x="XLAT", y="bottom_top", ax=left_ax, cmap="coolwarm", add_colorbar=True)
right_diff.plot(
    x="XLAT", y="bottom_top", ax=right_ax, cmap="coolwarm", add_colorbar=True
)

# Add surface difference to the center plot
center_ax = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())
center_ax.set_title("Surface Pressure Difference")
surface_diff = diag["met_em_pres"].isel(bottom_top=0) - diag["wrf_pres"].isel(
    bottom_top=0
)
surface_diff.plot(x="XLONG", y="XLAT", ax=center_ax, cmap="coolwarm", add_colorbar=True)

fig.tight_layout()
# %%
