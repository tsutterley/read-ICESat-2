scp_ICESat2_files.py
====================

 - Copies ICESat-2 HDF5 data from between a local host and a remote host

#### Calling Sequence
```bash
python scp_ICESat2_files.py --host=<host> --user=<username> \
    --product=ATL06 --release=205 --granule=10,11,12 --cycle=1,2 \
    --remote=<path_to_remote> --verbose --mode=0o775
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/scripts/scp_ICESat2_files.py)  

#### Command Line Options
 - `-h`, `--help`: list the command line options
 - `--host=X`: Remote server host
 - `--user=X`: Remote server username
 - `-D X`, `--directory=X`: Local working directory
 - `--remote=X`: Remote working directory
 - `--product=X`: ICESat-2 data product to copy
 - `--release=X`: ICESat-2 data release to copy
 - `--version=X`: ICESat-2 data version to copy
 - `--granule=X`: ICESat-2 granule regions to copy
 - `--cycle=X`: ICESat-2 cycle to copy
 - `--track=X`: ICESat-2 tracks to copy
 - `-C`, `--clobber`: overwrite existing data in transfer
 - `-V`, `--verbose`: output information about each synced file
 - `-M X`, `--mode=X`: permission mode of directories and files synced
 - `--push`: Transfer files from local computer to remote server
 - `-L`, `--list`: only list files to be transferred
