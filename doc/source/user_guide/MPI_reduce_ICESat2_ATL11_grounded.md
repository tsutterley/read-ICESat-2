MPI_reduce_ICESat2_ATL11_ice_shelves.py
=======================================

- Create masks for reducing ICESat-2 data into grounded ice regions
- [IceBridge BedMachine Greenland, Version 4](https://doi.org/10.5067/VLJ5YXKCNGXO)
- [MEaSUREs Antarctic Boundaries for IPY 2007-2009](https://doi.org/10.5067/AXE4121732AD)

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_reduce_ICESat2_ATL11_grounded.py <path_to_ATL11_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL11_grounded.py)

#### Inputs
1. `ATL11_file`: full path to ATL11 file

#### Command Line Options
- `-D X`, `--directory X`: Working data directory for grounded ice shapefiles
- `-B X`, `--buffer X`: Distance in kilometers to buffer grounded ice mask
- `-V`, `--verbose`: output module information for process
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets
