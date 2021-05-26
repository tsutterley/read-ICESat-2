symbolic_ICESat2_files.py
=========================

- Creates symbolic links for ICESat-2 HDF5 files organized by date

#### Calling Sequence
```bash
python symbolic_ICESat2_files.py --product ATL06 --release 003 \
    --granule 10 11 12 --cycle 1 2 --directory <path_to_directory>
    --scf_outgoing <path_to_outgoing> --verbose --mode 0o775
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/symbolic_ICESat2_files.py)

#### Command Line Options
- `-h`, `--help`: list the command line options
- `-D X`, `--directory X`: local working directory for creating symbolic links
- `-p X`, `--product X`: ICESat-2 data product to create symbolic links
- `-r X`, `--release X`: ICESat-2 data release to create symbolic links
- `-v X`, `--version X:` ICESat-2 data version to create symbolic links
- `-c X`, `--cycle X`: ICESat-2 orbital cycle to create symbolic links
- `-t X`, `--track X`: ICESat-2 reference ground tracks to create symbolic links
- `-g X`, `--granule X`: ICESat-2 granule regions to create symbolic links
- `--scf_incoming X`: directory on the SCF where the rscf sends PANS
- `--scf_outgoing X`: directory on the SCF where the data resides
- `-V`, `--verbose`: output information about each symbolic link
- `-M X`, `--mode X`: permission mode of directories
