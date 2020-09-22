read_ICESat2_ATL03.py
=====================

- Reads ICESat-2 ATL03 Global Geolocated Photons and ATL09 Atmospheric Characteristics data files  

#### Calling Sequence
```python
from icesat2_toolkit.read_ICESat2_ATL03 import read_HDF5_ATL03
IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams = read_HDF5_ATL03(FILENAME)
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL03.py)  

#### Inputs
 1. `FILENAME`: full path to ATL03 file

#### Options
 - `ATTRIBUTES`: read file, group and variable attributes
 - `VERBOSE`: output information about top-level HDF5 groups

#### Outputs
 - `IS2_atl03_mds`: dictionary with ATL03 variables
 - `IS2_atl03_attrs`: dictionary with ATL03 attributes
 - `IS2_atl03_beams`: list with valid ICESat-2 beams within ATL03 file
