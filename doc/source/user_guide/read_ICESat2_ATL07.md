read_ICESat2_ATL07.py
=====================

 - Reads ICESat-2 ATL07 Sea Ice Height data files  

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL07 import read_HDF5_ATL07
IS2_atl07_mds,IS2_atl07_attrs,IS2_atl07_beams = read_HDF5_ATL07(FILENAME)
```

#### Inputs
 1. `FILENAME`: full path to ATL07 file

#### Options
 - `ATTRIBUTES`: read file, group and variable attributes
 - `VERBOSE`: output information about top-level HDF5 groups

#### Outputs
 - `IS2_atl07_mds`: dictionary with ATL07 variables
 - `IS2_atl07_attrs`: dictionary with ATL07 attributes
 - `IS2_atl07_beams`: list with valid ICESat-2 beams within ATL07 file
