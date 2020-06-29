MPI_ICESat2_ATL03_histogram.py
==============================

- Read ICESat-2 ATL03 and ATL09 data files to calculate average segment surfaces  
    * ATL03 datasets: Global Geolocated Photons  
    * ATL09 datasets: Atmospheric Characteristics  
- Alternative algorithm that uses gaussian/generalized gaussian decomposition to extract possibly multiple height surfaces from a histogram of photon events  

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_ICESat2_ATL03_histogram.py <path_to_ATL03_file> <path_to_ATL09_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/scripts/MPI_ICESat2_ATL03_histogram.py)  

#### Inputs
1. `ATL03_file`: full path to ATL03 file  
2. `ATL09_file`: full path to ATL09 file  

#### Command Line Options
- `-V`, `--verbose`: output module information for process  
- `-O X`, `--output=X`: full path to output file  
- `-M X`, `--mode=X`: permissions mode of output HDF5 datasets  
