# Interpolator for WRF-CHEM

Interpolator for [WRF-CHEM](https://github.com/wrf-model/wrf) is a preprocessing application for using global chemistry model fields with WRF-CHEM. It interpolates the global chemistry model fields to the WRF-CHEM grid and vertical levels. Think of it like WPS but for chemistry fields.

Some features:

- Interpolates 3D fields
- Can do species/size bin mapping through linear combinations of source fields
- Potentially support many global models (currently supports CAMS EAC4)

The application is written in Python and is meant to be used a command-line tool. Some potential limitations:

- The application currently assumes that fields are mixing ratios.
- The global model fields should be on a regular lat-lon grid.

## Installation

The application is available on PyPI and can be installed using pip:

```bash
pip install interpolator-for-wrfchem
```

## Usage

The workflow for using the interpolator is as follows:

1. Use WPS and `real.exe` as usual to generate the `met_em`, `wrfinput`, and `wrfbdy` files.
2. Download global chemistry model fields (e.g. CAMS EAC4) for the same time period as the WRF simulation.
3. Run the interpolator to interpolate the global chemistry model fields to the WRF-CHEM grid and vertical levels.
4. Run WRF-CHEM.

The interpolator will update `wrfinput` and `wrfbdy` files to include the chemistry information.

![Workflow](./workflow.drawio.png)

The interpolator is a command-line tool and can be run as follows:

```bash
interpolator-for-wrfchem <global model name> <global model data path> <met_em path> <species map path> <wrfinput path>
```

The `wrfinput` and `wrfbdy` files **WILL BE MODIFIED**! The `global model name` can be one of the following:

- `cams_eac4`: CAMS EAC4 data (w/ 60 vertical levels)
- `cams_eac4_adsbeta``: CAMS EAC4 data from the [ADS-Beta](https://ads-beta.atmosphere.copernicus.eu/)
- `cams_global_forecasts`: CAMS global forecasts (w/ 137 vertical levels)

There are some optional flags:

- `--wrfbdy=`: Path to the `wrfbdy_d01` file, if not provided, the boundary is not updated.
- `--copy-icbc`: Make a backup of the `wrfinput` and `wrfbdy` files before updating them.
- `--no-ic`: Do not update the `wrfinput`. You must nonetheless provide the path to the `wrfinput` file as it is required to read some information.
- `--diagnostics`: Store some diagnostic information in the `diag_cams_interp.nc` file.

When you use nested domains, you can run the application multiple times, each time pointing to a different `wrfinput` file. You can omit the `wrfbdy` file when running the application for the nested domains' `wrfinput` files.
If you need to update `wrfbdy` files for the future without touching `wrfinput` (e.g. for a cycling run), point to a correct `wrfinput` file (correct means it's the same model grid and configuration) and use the `--no-ic` flag.

## Species mapping

In many cases, the available fields of the global model do not directly correspond to the ones used by the chemistry/dust scheme you want to use in WRF-CHEM. For example, you might have dust concentrations available in different size bins. The application supports "species mapping", through which the WRF-CHEM fields are created through a linear combination of global model fields, after interpolation.

Detailed description of the species file format is available in [species_maps.md](./species_maps/species_maps.md).

## Use in other projects

To automate use of the interpolator by using it as a library in your project, it is recommended to use the `do_initial_conditions()` and `do_boundary_conditions()` functions from the `interpolator_for_wrfchem` module. Provided with the correct arguments, these functions will mirror the command-line behavior of the interpolator. Exact details of how to set up all the magic objects required by these functions are available inside the main CLI entrypoint (file `interpolator_for_wrfchem/__init__.py`, `main()` function).

An example usage would be:

```python
# Contains information on how to map from global model species to WRF-CHEM species
mapping = SpeciesMap(mapping_path)

# Represents a wrfinput file. The original netCDF file can be accessed as a
# `netcdf4.Dataset` in the `WRFInput.nc_file` attribute.
wrf = WRFInput(wrfinput_path)

# You can get the main coordinates in dataset form by calling wrf.get_dataset()
wrf_ds = wrf.get_dataset()

# Represents a set of global model files. Available choices in `global_models/__init__.py`.
global_model = CAMS_ADSBeta_Base(input_files, mapping.required_source_species)
global_model_ds = global_model.get_dataset(wrf.time)

# Do initial condition interpolation and MODIFY wrfinput_d01
do_initial_conditions(wrf, wrf_ds, global_model_ds, mapping, diagnostics)
wrf.close()

# Represents a folder of met_em files. Required only for the boundary conditions.
met_em = MetEm(met_em_path, wrf_ds)
# Do boundary condition interpolation and MODIFY wrfbdy_d01
do_boundary_conditions(wrfbdy, met_em, wrf_ds, global_model_ds, mapping)
```

The interpolation routines are available inside the `interpolation.py` file and are applied to xarray Datasets, so they might be useful in other projects as well.


## License

The interpolator is licensed under the MIT License. See [LICENSE](./LICENSE) for more information.
Please cite the project if you use it for your research!
