read_ICESat2_ATL10.py
=====================

- Reads ICESat-2 ATL10 Sea Ice Freeboad data files

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL10 import read_HDF5_ATL10
IS2_atl10_mds,IS2_atl10_attrs,IS2_atl10_beams = read_HDF5_ATL10(FILENAME)
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL10.py)

#### Arguments
1. `FILENAME`: full path to ATL10 file

#### Keyword arguments
- `ATTRIBUTES`: read file, group and variable attributes
- `VERBOSE`: output information about top-level HDF5 groups

#### Returns
- `IS2_atl10_mds`: dictionary with ATL10 variables
- `IS2_atl10_attrs`: dictionary with ATL10 attributes
- `IS2_atl10_beams`: list with valid ICESat-2 beams within ATL10 file
