#!/usr/bin/env python
u"""
read_ICESat2_ATL07.py (12/2022)
Read ICESat-2 ATL07 (Sea Ice Height) data files

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
    Written 11/2019
"""
import warnings
import icesat2_toolkit.io

# PURPOSE: read ICESat-2 ATL07 HDF5 data files
def read_HDF5_ATL07(*args, **kwargs):
    """
    Reads ICESat-2 ATL07 (Sea Ice Height) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL07 file
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables

    Returns
    -------
    IS2_atl07_mds: dict
        ATL07 variables
    IS2_atl07_attrs:
        ATL07 attributes
    IS2_atl07_beams: list
        valid ICESat-2 beams within ATL07 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL07.read_granule(*args, **kwargs)

# PURPOSE: find valid beam groups within ICESat-2 ATL07 HDF5 data files
def find_HDF5_ATL07_beams(*args, **kwargs):
    """
    Find valid beam groups within ICESat-2 ATL07 (Sea Ice Height) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL07 file

    Returns
    -------
    IS2_atl07_beams: list
        valid ICESat-2 beams within ATL07 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL07.find_beams(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL07 HDF5 data files for beam variables
def read_HDF5_ATL07_beam(*args, **kwargs):
    """
    Reads ICESat-2 ATL07 (Sea Ice Height) data files for a specific beam

    Parameters
    ----------
    FILENAME: str
        full path to ATL07 file
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

    Returns
    -------
    IS2_atl07_mds: dict
        ATL07 variables
    IS2_atl07_attrs:
        ATL07 attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL07.read_beam(*args, **kwargs)
