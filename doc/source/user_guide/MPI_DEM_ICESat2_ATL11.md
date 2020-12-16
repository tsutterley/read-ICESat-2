MPI_DEM_ICESat2_ATL11.py
========================

- Determines which digital elevation model tiles to read for a given ATL11 file  
- Reads 3x3 array of tiles for points within bounding box of central mosaic tile  
- Interpolates digital elevation model to locations of ICESat-2 ATL11 segments  

- ArcticDEM 2m digital elevation model tiles  
    * http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/  
    * http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/  

- REMA 8m digital elevation model tiles  
    * http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/  
    * http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/  

- GIMP 30m digital elevation model tiles  
    * https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/  

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_DEM_ICESat2_ATL11.py <path_to_ATL11_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_DEM_ICESat2_ATL11.py)  

#### Inputs
1. `ATL11_file`: full path to ATL11 file  

#### Command Line Options
- `-D X`, `--directory X`: Working data directory for elevation models
- `--model X`: Set the digital elevation model (REMA, ArcticDEM, GIMP) to run
- `-V`, `--verbose`: output module information for process  
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets  
