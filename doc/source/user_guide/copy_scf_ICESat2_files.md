copy_scf_ICESat2_files.py
=========================

 - Copies ICESat-2 HDF5 files from the SCF server  

#### Calling Sequence
```bash
python copy_scf_ICESat2_files.py --scf_host=<host> --scf_user=<username> \
	--product=ATL06 --release=001 --granule=10,11,12 --cycle=1,2 \
	--scf_outgoing=<path_to_outgoing> --verbose --mode=0o775
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/scripts/copy_scf_ICESat2_files.py)  

#### Command Line Options
 - `-h`, `--help`: list the command line options  
 - `--scf_host=X`: hostname of the SCF server  
 - `--scf_user=X`: SCF server username  
 - `-D X`, `--directory=X`: local working directory for receiving data  
 - `--product=X`: ICESat-2 data product to copy  
 - `--release=X`: ICESat-2 data release to copy  
 - `--version=X`: ICESat-2 data version to copy  
 - `--granule=X`: ICESat-2 granule regions to copy  
 - `--cycle=X`: ICESat-2 cycle to copy  
 - `--track=X`: ICESat-2 tracks to copy  
 - `--scf_incoming=X`: directory on the SCF where the rscf sends PANS  
 - `--scf_outgoing=X`: directory on the SCF where the data resides  
 - `-C`, `--clobber`: overwrite existing data in transfer  
 - `-V`, `--verbose`: output information about each synced file  
 - `-M X`, `--mode=X`: permission mode of directories and files synced  
 - `-L`, `--list`: only list files to be transferred  
