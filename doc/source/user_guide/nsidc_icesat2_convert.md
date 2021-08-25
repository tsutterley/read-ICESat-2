nsidc_icesat2_convert.py
========================

- Syncs all available ICESat-2 data converting from HDF5 to zarr or rechunked HDF5 files for a specified data product, release, granule and track.
- zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores
- rechunked HDF5 files can be more optimized for cloud-based object stores
- The first time we run the script, it will copy the necessary dataset in the selected local directory.
- If we already have all the data, and we run the script again: only files added or modified on the remote server will downloaded.

#### Calling Sequence
```bash
python nsidc_icesat2_convert.py --user <username> --directory <outgoing> \
	--release 001 --granule 10 11 12 --format zarr --mode 0o775 ATL06
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_convert.py)

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for alternative authentication
- `-D X`, `--directory`: local working directory for receiving data
- `-Y X`, `--year X`: years to sync
- `-S X`, `--subdirectory X`: specific subdirectories to sync
- `-r X`, `--release X`: ICESat-2 data release to sync
- `-v X`, `--version X:` ICESat-2 data version to sync
- `-t X`, `--track X`: ICESat-2 reference ground tracks to sync
- `-g X`, `--granule X`: ICESat-2 granule regions to sync
- `-f X`, `--format X`: output file format
	* `'zarr'`
	* `'HDF5'`
- `-c X`, `--chunks X`: Rechunk files to size
- `--auxiliary`: Sync ICESat-2 auxiliary files for each HDF5 file
- `-I X, --index X`: Input index of ICESat-2 files to sync
- `-F`, `--flatten`: Do not create subdirectories
- `-P X`, `--np X`: Number of processes to use in file downloads
- `-M X`, `--mode X`: Local permissions mode of the directories and files synced
- `-T X`, `--timeout X`: Timeout in seconds for blocking operations
- `-R X`, `--retry X`: Connection retry attempts
- `--log`: output log of files downloaded
- `--list`: print files to be transferred, but do not execute transfer
- `-C`, `--clobber`: Overwrite existing data in transfer
