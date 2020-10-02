MPI_reduce_ICESat2_ATL06_drainages.py
=====================================

- Create masks for reducing ICESat-2 data into [IMBIE-2 drainage regions](http://imbie.org/imbie-2016/drainage-basins/)  

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_reduce_ICESat2_ATL06_drainages.py <path_to_ATL06_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL06_drainages.py)  

#### Inputs
1. `ATL06_file`: full path to ATL06 file  

#### Command Line Options
- `-D X`, `--directory X`: Working data directory for the IMBIE-2 drainage shapefiles
- `-V`, `--verbose`: output module information for process  
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets  
