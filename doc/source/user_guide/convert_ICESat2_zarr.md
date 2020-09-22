convert_ICESat2_zarr.py
=======================

 - Converts ICESat-2 HDF5 datafiles to zarr datafiles for a specified data product, release, granule and track.  
 - zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores

#### Calling Sequence
```bash
python convert_ICESat2_zarr.py --directory=<outgoing> \
	--release=001 --granule=10,11,12 --mode=0o775 ATL06
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/convert_ICESat2_zarr.py)  

#### Command Line Options
 - `-D X`, `--directory`: working data directory  
 - `--release=X`: ICESat-2 data release to run  
 - `--version=X`: ICESat-2 data version to run  
 - `--track=X`: ICESat-2 reference ground tracks to run  
 - `--granule=X`: ICESat-2 granule regions to run  
 - `-P X`, `--np=X`: Number of processes to use in file conversion  
 - `-M X`, `--mode=X`: Local permissions mode of the converted files  
 - `-C`, `--clobber`: Overwrite existing data in transfer  
