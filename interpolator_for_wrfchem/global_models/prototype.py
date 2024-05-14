import datetime as dt
from abc import ABC

import xarray as xr


class GlobalModel(ABC):
    """Prototype class for global models."""

    def get_dataset(self, time: dt.datetime) -> xr.Dataset:
        """Return the required variables in a xarray dataset for the given time."""
        raise NotImplementedError("Subclasses must implement this method")

    @property
    def available_times(self) -> list[dt.datetime]:
        """Return the list of available times."""
        raise NotImplementedError("Subclasses must implement this method")
