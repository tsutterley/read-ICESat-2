read_ICESat2_ATL11.py
=====================

 - Reads ICESat-2 ATL11 Annual Land Ice Height data files

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11
IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_pairs = read_HDF5_ATL11(FILENAME)
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL11.py)

#### Inputs
 1. `FILENAME`: full path to ATL11 file

#### Options
 - `ATTRIBUTES`: read file, group and variable attributes
 - `REFERENCE`: read variables from ATL11 `ref_surf` group
 - `CROSSOVERS`: read variables from ATL11 `crossing_track_data` group
 - `VERBOSE`: output information about top-level HDF5 groups

#### Outputs
 - `IS2_atl06_mds`: dictionary with ATL11 variables
 - `IS2_atl06_attrs`: dictionary with ATL11 attributes
 - `IS2_atl06_pairs`: list with valid ICESat-2 beam pairs within ATL11 file
