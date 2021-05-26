read_ICESat2_ATL06.py
=====================

- Reads ICESat-2 ATL06 Land Ice Elevation data files

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL06 import read_HDF5_ATL06
IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_beams = read_HDF5_ATL06(FILENAME)
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL06.py)

#### Arguments
1. `FILENAME`: full path to ATL06 file

#### Keyword arguments
- `ATTRIBUTES`: read file, group and variable attributes
- `HISTOGRAM`: read variables from ATL06 `residual_histogram` group
- `QUALITY`: read variables from ATL06 `segment_quality` group
- `VERBOSE`: output information about top-level HDF5 groups

#### Returns
- `IS2_atl06_mds`: dictionary with ATL06 variables
- `IS2_atl06_attrs`: dictionary with ATL06 attributes
- `IS2_atl06_beams`: list with valid ICESat-2 beams within ATL06 file
