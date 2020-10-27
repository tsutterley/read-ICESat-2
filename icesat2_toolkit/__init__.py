"""
An ICESat-2 toolkit for Python
==============================

icesat2_toolkit contains Python tools for working with data
from the NASA Ice, Cloud and land Elevation Satellite-2 (ICESat-2)

The package works using Python packages (numpy, scipy, scikit-learn, shapely)
combined with data storage in HDF5 and zarr, and mapping with
matplotlib and cartopy

It aims to be a simple and efficient solution for using data from the
ICESat-2 mission and to support its science applications

Documentation is available at https://read-icesat-2.readthedocs.io
"""
import icesat2_toolkit.time
import icesat2_toolkit.utilities
from icesat2_toolkit.convert import convert
from icesat2_toolkit.convert_julian import convert_julian
from icesat2_toolkit.convert_calendar_decimal import convert_calendar_decimal
from icesat2_toolkit.convert_delta_time import convert_delta_time
from icesat2_toolkit.read_ICESat2_ATL03 import read_HDF5_ATL03, \
    find_HDF5_ATL03_beams, read_HDF5_ATL09, read_HDF5_ATL03_main, read_HDF5_ATL03_beam
from icesat2_toolkit.read_ICESat2_ATL06 import read_HDF5_ATL06, \
    find_HDF5_ATL06_beams, read_HDF5_ATL06_beam
from icesat2_toolkit.read_ICESat2_ATL07 import read_HDF5_ATL07, \
    find_HDF5_ATL07_beams
from icesat2_toolkit.read_ICESat2_ATL12 import read_HDF5_ATL12, \
    find_HDF5_ATL12_beams
