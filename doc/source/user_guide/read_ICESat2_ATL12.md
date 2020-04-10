read_ICESat2_ATL12.py
=====================

 - Reads ICESat-2 ATL12 Ocean Surface Height data files  

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL12 import read_HDF5_ATL12
IS2_atl12_mds,IS2_atl12_attrs,IS2_atl12_beams = read_HDF5_ATL12(FILENAME)
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/master/icesat2_toolkit/read_ICESat2_ATL12.py)  

#### Inputs
 1. `FILENAME`: full path to ATL12 file

#### Options
 - `ATTRIBUTES`: read file, group and variable attributes
 - `VERBOSE`: output information about top-level HDF5 groups

#### Outputs
 - `IS2_atl12_mds`: dictionary with ATL12 variables
 - `IS2_atl12_attrs`: dictionary with ATL12 attributes
 - `IS2_atl12_beams`: list with valid ICESat-2 beams within ATL12 file
