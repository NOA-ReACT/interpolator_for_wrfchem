"""
Check interpolation of the O3 field on the surface level (for CAMS NCP)
"""

# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import ipywidgets

# %%
wrfinput = xr.open_dataset("/home/thgeorgiou/data/test_o3/wrf/wrfinput_d02").isel(
    Time=0
)
wrfinput_new = xr.open_dataset(
    "/home/thgeorgiou/data/test_o3/wrf/wrfinput_d02_new"
).isel(Time=0)
wrfinput

# %%
min_lat = wrfinput.XLAT.min().values.round()
max_lat = wrfinput.XLAT.max().values.round()
min_lon = wrfinput.XLONG.min().values.round()
max_lon = wrfinput.XLONG.max().values.round()

# %% cams ic
cams = (
    xr.open_dataset("/home/thgeorgiou/data/test_o3/CAMS/2025-05-06/data_mlev.nc")
    .stack(time=["forecast_period", "forecast_reference_time"])
    .set_index(time="valid_time")
    .sel(time="2025-05-06T00:00:00")
    .sel(longitude=slice(min_lon, max_lon), latitude=slice(max_lat, min_lat))
)
cams

# %%
cams_sfc = (
    xr.open_dataset("/home/thgeorgiou/data/test_o3/CAMS/2025-05-06/data_sfc.nc")
    .stack(time=["forecast_period", "forecast_reference_time"])
    .set_index(time="valid_time")
    .sel(time="2025-05-06T00:00:00")
    .sel(longitude=slice(min_lon, max_lon), latitude=slice(max_lat, min_lat))
)
cams_sfc

# %% Compute surface pressure from wrfinput
psfc = (wrfinput.MU + wrfinput.MUB + wrfinput.P_TOP) / 100  # Convert from Pa to hPa
psfc

# %%
fig, axes = plt.subplots(
    1, 3, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(15, 3)
)
psfc.plot(x="XLONG", y="XLAT", cmap="viridis", vmin=750, vmax=1020, ax=axes[0])
# (cams_sfc["sp"] / 100).plot(
#     x="longitude", y="latitude", cmap="viridis", vmin=750, vmax=1020, ax=axes[1]
# )
cams_sfc_sp_interp = (cams_sfc["sp"] / 100).interp(
    longitude=wrfinput.XLONG, latitude=wrfinput.XLAT
)
diff = psfc - cams_sfc_sp_interp

cams_sfc_sp_interp.plot(
    x="XLONG", y="XLAT", cmap="viridis", vmin=750, vmax=1020, ax=axes[1]
)
diff.plot(x="XLONG", y="XLAT", cmap="coolwarm", vmin=-10, vmax=10, ax=axes[2])

for ax in axes:
    ax.coastlines()
    ax.set_title("Surface Pressure (hPa)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

# %%
fig, axes = plt.subplots(
    1,
    3,
    subplot_kw={"projection": ccrs.PlateCarree()},
    figsize=(14, 4),
)

wrfinput["o3"].isel(bottom_top=0).plot(
    x="XLONG", y="XLAT", cmap="viridis", vmin=0, vmax=0.1, ax=axes[0]
)
wrfinput_new["o3"].isel(bottom_top=0).plot(
    x="XLONG", y="XLAT", cmap="viridis", vmin=0, vmax=0.1, ax=axes[1]
)
(cams["go3"] * 0.60345388 * 1e6).isel(model_level=-1).interp(
    latitude=wrfinput.XLAT, longitude=wrfinput.XLONG
).plot(
    x="XLONG",
    y="XLAT",
    cmap="viridis",
    vmin=0,
    vmax=0.1,
    ax=axes[2],
)

for ax in axes:
    ax.coastlines()


# %%
p = (
    (cams["go3"] * 0.60345388 * 1e6)
    .isel(model_level=-1)
    .plot(
        x="longitude",
        y="latitude",
        cmap="viridis",
        subplot_kws={"projection": ccrs.PlateCarree()},
        vmin=0,
        vmax=0.1,
    )
)


p.axes.coastlines()
# %%
