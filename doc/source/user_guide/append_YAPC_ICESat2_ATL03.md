append_YAPC_ICESat2_ATL03.py
============================

 - Reads ICESat-2 [ATL03 geolocated photon height product](https://nsidc.org/data/ATL03) and appends photon classification flags from [YAPC](https://github.com/tsutterley/yapc) (*Yet Another Photon Classifier*)

#### Calling Sequence
```bash
python append_YAPC_ICESat2_ATL03.py --verbose --mode 0o775 <path_to_ATL03_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/append_YAPC_ICESat2_ATL03.py)

#### Inputs
1. `ATL03_file`: full path to ATL03 file

#### Command Line Options
- `-V`, `--verbose`: output module information for process
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets
