reduce_ICESat2_ATL06_raster.py
==============================

- Create masks for reducing ICESat-2 ATL06 data using raster imagery

#### Calling Sequence
```bash
python3 reduce_ICESat2_ATL06_raster.py <path_to_ATL06_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/reduce_ICESat2_ATL06_raster.py)

#### Inputs
1. `ATL06_file`: full path to ATL06 file

#### Command Line Options
- `-R X`, `--raster X`: Input raster file
- `-F X`, `--format X`: Input raster file format
    * `'netCDF4'`
    * `'HDF5'`
    * `'geotiff'`
- `-v X`, `--variables X`: variable names of data in HDF5 or netCDF4 file
        `x`, `y` and `data` variable names
- `-P X`, `--projection X`: spatial projection as EPSG code or PROJ4 string
- `-O X`, `--output X`: Output mask file name
- `-V`, `--verbose`: output module information for process
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets
