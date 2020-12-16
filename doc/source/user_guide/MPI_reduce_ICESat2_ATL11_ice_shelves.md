MPI_reduce_ICESat2_ATL11_ice_shelves.py
=======================================

- Create masks for reducing ICESat-2 data into [floating ice shelves](http://dx.doi.org/10.5067/AXE4121732AD)  

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_reduce_ICESat2_ATL11_ice_shelves.py <path_to_ATL11_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL11_ice_shelves.py)  

#### Inputs
1. `ATL11_file`: full path to ATL11 file  

#### Command Line Options
- `-D X`, `--directory X`: Working data directory for ice shelf shapefiles
- `-B X`, `--buffer X`: Distance in kilometers to buffer ice shelves mask
- `-V`, `--verbose`: output module information for process  
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets  
