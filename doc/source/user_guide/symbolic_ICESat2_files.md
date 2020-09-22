symbolic_ICESat2_files.py
=========================

 - Creates symbolic links for ICESat-2 HDF5 files organized by date

#### Calling Sequence
```bash
python symbolic_ICESat2_files.py --product=ATL06 --release=001 \
    --granule=10,11,12 --cycle=1,2 --directory=<path_to_directory>
    --scf_outgoing=<path_to_outgoing> --verbose --mode=0o775
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/symbolic_ICESat2_files.py)  

#### Command Line Options
 - `-h`, `--help`: list the command line options
 - `-D X`, `--directory=X`: local working directory for creating symbolic links
 - `--product=X`: ICESat-2 data product to create symbolic links
 - `--release=X`: ICESat-2 data release to create symbolic links
 - `--version=X`: ICESat-2 data version to create symbolic links
 - `--granule=X`: ICESat-2 granule regions to create symbolic links
 - `--cycle=X`: ICESat-2 cycle to create symbolic links
 - `--track=X`: ICESat-2 tracks to create symbolic links
 - `--scf_incoming=X`: directory on the SCF where the rscf sends PANS
 - `--scf_outgoing=X`: directory on the SCF where the data resides
 - `-V`, `--verbose`: output information about each symbolic link
 - `-M X`, `--mode=X`: permission mode of directories
