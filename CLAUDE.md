# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Interpolator for WRF-CHEM is a preprocessing application that interpolates global chemistry model fields to WRF-CHEM grids and vertical levels. It functions as a chemistry-field equivalent to WPS (WRF Preprocessing System). The application reads global model data (currently CAMS EAC4 or CAMS Global Forecasts), interpolates horizontally and vertically to the WRF grid, applies species mapping via linear combinations, and updates `wrfinput` and `wrfbdy` files in-place.

## Commands

```bash
# Install dependencies
uv sync

# Run the CLI
uv run interpolator-for-wrfchem <global_model> <input_files> <mapping> <wrfinput>

# Example
uv run interpolator-for-wrfchem cams_eac4 /path/to/cams species_map.toml wrfinput_d01 --wrfbdy=wrfbdy_d01
```

No formal test suite exists. Testing is done via exploration scripts in `exploration/`.

## Architecture

### Core Processing Flow

1. **CLI Entry** (`__init__.py:main`) - Click-based CLI parses arguments and orchestrates workflow
2. **Species Mapping** (`species_map.py:SpeciesMap`) - Loads TOML config defining how global model species map to WRF-CHEM species via linear combinations with unit conversion
3. **WRF File Abstraction** (`wrf.py:WRFInput`, `WRFBoundary`) - Wraps netCDF4.Dataset for WRF files, provides coordinates as xarray Dataset via `get_dataset()`
4. **Global Model Abstraction** (`global_models/`) - `GlobalModel` ABC with implementations for CAMS. Returns data as xarray Dataset via `get_dataset(time)`
5. **Interpolation** (`interpolation.py`) - Two-step process: horizontal interpolation with RectBivariateSpline, then vertical linear interpolation to WRF pressure levels
6. **Output** - Updates wrfinput/wrfbdy NetCDF files directly

### Key Library Functions

- `do_initial_conditions(wrf, wrf_ds, global_model_ds, mapping, diagnostics)` - Interpolates and writes to wrfinput
- `do_boundary_conditions(wrfbdy_path, wrf_ds, global_model_ds, mapping)` - Interpolates lateral boundaries and computes time tendencies

### Adding New Global Models

1. Create subclass of `GlobalModel` in `global_models/`
2. Implement `get_dataset(time)` returning xarray Dataset and `available_times` property
3. Register in `GLOBAL_MODELS` dict in `global_models/__init__.py`

### Species Map Format

TOML files in `species_maps/` define mappings. Key sections:
- `[units]` - source/target units (kg, g, mg, ug)
- `[aliases_source]`/`[aliases_target]` - optional name remapping
- `[species_map]` - linear combination coefficients: `TARGET = { SRC1 = 1.0, SRC2 = 0.5 }`

See `species_maps/species_maps.md` for full specification.

## Important Assumptions

- Fields are mixing ratios
- Global model fields are on regular lat-lon grid
- WRF files use standard field names and dimensions: (Time, bottom_top, south_north, west_east)
- Global model data organized as `YYYY-MM-DD/data_sfc.nc`, `data_mlev.nc` subdirectories
