nsidc_icesat2_zarr.py
=====================

 - Syncs all available ICESat-2 data converting from HDF5 to zarr files for a specified data product, release, granule and track.  
 - zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores
 - The first time we run the script, it will copy the necessary dataset in the selected local directory.  
 - If we already have all the data, and we run the script again: only files added or modified on the remote server will downloaded.  

#### Calling Sequence
```bash
python nsidc_icesat2_zarr.py --user=<username> --directory=<outgoing> \
	--release=001 --granule=10,11,12 --mode=0o775 ATL06
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/nsidc_icesat2_zarr.py)  

#### Command Line Options
 - `-U X`, `--user=X`: username for NASA Earthdata Login  
 - `-N X`, `--netrc=X`: path to .netrc file for alternative authentication  
 - `-D X`, `--directory`: local working directory for receiving data  
 - `--release=X`: ICESat-2 data release to sync  
 - `--version=X`: ICESat-2 data version to sync  
 - `--track=X`: ICESat-2 reference ground tracks to sync  
 - `--granule=X`: ICESat-2 granule regions to sync  
 - `--auxiliary`: Sync ICESat-2 auxiliary files for each HDF5 file  
 - `-F`, `--flatten`: Do not create subdirectories  
 - `-P X`, `--np=X`: Number of processes to use in file downloads  
 - `-M X`, `--mode=X`: Local permissions mode of the directories and files synced  
 - `--log`: output log of files downloaded  
 - `--list`: print files to be transferred, but do not execute transfer  
 - `-C`, `--clobber`: Overwrite existing data in transfer  
