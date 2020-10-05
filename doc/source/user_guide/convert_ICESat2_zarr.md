convert_ICESat2_zarr.py
=======================

 - Converts ICESat-2 HDF5 datafiles to zarr datafiles for a specified data product, release, granule and track.  
 - zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores

#### Calling Sequence
```bash
python convert_ICESat2_zarr.py --directory <outgoing> \
	--release 003 --granule 10 11 12 --mode 0o775 ATL06
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/convert_ICESat2_zarr.py)  

#### Command Line Options
 - `-D X`, `--directory`: working data directory  
 - `-Y X`, `--year X`: years to run separated by commas
 - `-S X`, `--subdirectory X`: subdirectories to run separated by commas
 - `-r X`, `--release X`: ICESat-2 data release to run
 - `-v X`, `--version X:` ICESat-2 data version to run
 - `-t X`, `--track X`: ICESat-2 reference ground tracks to run
 - `-g X`, `--granule X`: ICESat-2 granule regions to run
 - `--chunks X`: Rechunk zarr files to size
 - `-P X`, `--np X`: Number of processes to use in file conversion  
 - `-C`, `--clobber`: Overwrite existing zarr files
 - `-V`, `--verbose`: Verbose output of processing run
 - `-M X`, `--mode X`: Local permissions mode of the converted files
