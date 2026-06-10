import datetime as dt
from abc import ABC
from collections import OrderedDict

import xarray as xr


class GlobalModel(ABC):
    """Prototype class for global models."""

    # Number of recently-read datasets to keep in memory. The boundary loop
    # processes times in increasing order, so the two bracketing source times
    # shift by at most one per WRF step — a 2-entry LRU keeps the current pair
    # hot and avoids reloading it for every dense WRF boundary step.
    _DS_CACHE_SIZE = 2

    # Perturbation: shift the fields by whole grid cells in (x, y), applied as
    # a periodic roll of the data while the coordinate axes stay fixed.
    # Positive values move the field toward increasing coordinate (east/north
    # on a lat-lon grid). Class-level default since subclass __init__ methods
    # don't chain to a base; set externally (e.g. from the CLI).
    hoz_shift: tuple[int, int] = (0, 0)

    def get_dataset(self, time: dt.datetime) -> xr.Dataset:
        """Return the required variables in a xarray dataset for the given time."""
        raise NotImplementedError("Subclasses must implement this method")

    def _get_dataset_cached(self, time: dt.datetime) -> xr.Dataset:
        """``get_dataset`` with a small bounded LRU cache, so repeated reads of
        the same source time (common when WRF boundary steps are denser than
        the global-model outputs) hit memory instead of disk."""
        # Lazily initialised: subclass __init__ methods don't chain to a base.
        cache = self.__dict__.setdefault("_ds_cache", OrderedDict())
        if time in cache:
            cache.move_to_end(time)
            return cache[time]
        ds = self.get_dataset(time)
        cache[time] = ds
        while len(cache) > self._DS_CACHE_SIZE:
            cache.popitem(last=False)
        return ds

    def _apply_hoz_shift(self, ds: xr.Dataset) -> xr.Dataset:
        """Roll all data variables by ``hoz_shift`` cells along the horizontal
        dimensions, leaving the coordinate axes untouched so they stay
        monotonic for the horizontal interpolation. Rolling by +n along an
        ascending axis moves the field toward increasing coordinate."""
        shift_x, shift_y = self.hoz_shift
        if shift_x == 0 and shift_y == 0:
            return ds
        dim_x = ds.attrs.get("hoz_coord_x", "longitude")
        dim_y = ds.attrs.get("hoz_coord_y", "latitude")
        return ds.roll({dim_x: shift_x, dim_y: shift_y}, roll_coords=False)

    @property
    def available_times(self) -> list[dt.datetime]:
        """Return the list of available times."""
        raise NotImplementedError("Subclasses must implement this method")

    def get_dataset_interpolated(self, time: dt.datetime) -> xr.Dataset:
        """Return the dataset at `time`, linearly interpolating in time between
        the two nearest available times when `time` is not itself available.

        Every ``get_dataset`` implementation returns the chemistry fields and
        the pressure fields (``pres``/``pres_hf``) as numeric data variables on
        a grid that is identical across times, while the grid axes are
        coordinates. A single linear blend ``ds0*(1-w) + ds1*w`` therefore
        interpolates every field (including the time-varying surface pressure)
        while leaving the coordinate axes untouched.

        Raises if `time` is outside the available range — temporal
        extrapolation is not supported.
        """

        times = sorted(self.available_times)
        if time in times:
            return self._apply_hoz_shift(self._get_dataset_cached(time))

        if time < times[0] or time > times[-1]:
            raise RuntimeError(
                f"Requested time {time} is outside the available global-model "
                f"range [{times[0]} .. {times[-1]}]; temporal extrapolation is "
                f"not supported"
            )

        t0 = max(t for t in times if t < time)
        t1 = min(t for t in times if t > time)
        w = (time - t0).total_seconds() / (t1 - t0).total_seconds()
        print(
            f"Temporal interpolation: {time} from {t0} and {t1} (weight {w:.3f})"
        )

        ds0 = self._get_dataset_cached(t0)
        ds1 = self._get_dataset_cached(t1)
        # keep_attrs preserves dataset/variable attrs from the left operand
        # (e.g. wrfout's hoz_coord_x/hoz_coord_y set in get_dataset).
        with xr.set_options(keep_attrs=True):
            return self._apply_hoz_shift(ds0 * (1.0 - w) + ds1 * w)
