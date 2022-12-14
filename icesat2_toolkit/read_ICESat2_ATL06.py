#!/usr/bin/env python
u"""
read_ICESat2_ATL06.py (12/2022)
Read ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files

OPTIONS:
    ATTRIBUTES: read HDF5 attributes for groups and variables
    HISTOGRAM: read ATL06 residual_histogram variables
    QUALITY: read ATL06 segment_quality variables

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 12/2022: place some imports behind try/except statements
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: add check if input streaming from bytes
    Updated 10/2020: add small function to find valid beam groups
    Updated 07/2020: added function docstrings
    Updated 11/2019: create attribute dictionaries but don't fill if False
        add function for reading only beam level variables
    Updated 02/2019: continued writing read program with first ATL03 release
    Written 07/2017
"""
import warnings
import icesat2_toolkit.io

# PURPOSE: read ICESat-2 ATL06 HDF5 data files
def read_HDF5_ATL06(*args, **kwargs):
    """
    Reads ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL06 file
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables
    HISTOGRAM: bool, default False
        read ATL06 residual_histogram variables
    QUALITY: bool, default False
        read ATL06 segment_quality variables

    Returns
    -------
    IS2_atl06_mds: dict
        ATL06 variables
    IS2_atl06_attrs:
        ATL06 attributes
    IS2_atl06_beams: list
        valid ICESat-2 beams within ATL06 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL06.read_granule(*args, **kwargs)

# PURPOSE: find valid beam groups within ICESat-2 ATL06 HDF5 data files
def find_HDF5_ATL06_beams(*args, **kwargs):
    """
    Find valid beam groups within ICESat-2 ATL06 (Land Ice Along-Track
    Height Product) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL06 file

    Returns
    -------
    IS2_atl06_beams: list
        valid ICESat-2 beams within ATL06 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL06.find_beams(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL06 HDF5 data files for beam variables
def read_HDF5_ATL06_beam(*args, **kwargs):
    """
    Reads ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files
    for a specific beam

    Parameters
    ----------
    FILENAME: str
        full path to ATL06 file
    gtx: str
        beam name based on ground track and position

            - ``'gt1l'``
            - ``'gt1r'``
            - ``'gt2l'``
            - ``'gt2r'``
            - ``'gt3l'``
            - ``'gt3r'``
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables
    HISTOGRAM: bool, default False
        read ATL06 residual_histogram variables
    QUALITY: bool, default False
        read ATL06 segment_quality variables

    Returns
    -------
    IS2_atl06_mds: dict
        ATL06 variables
    IS2_atl06_attrs:
        ATL06 attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL06.read_beam(*args, **kwargs)