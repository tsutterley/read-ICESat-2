#!/usr/bin/env python
u"""
read_ICESat2_ATL10.py (12/2022)
Read ICESat-2 ATL10 (Sea Ice Freeboard) data files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 12/2022: place some imports behind try/except statements
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 12/2021
"""
import warnings
import icesat2_toolkit.io

# PURPOSE: read ICESat-2 ATL10 HDF5 data files
def read_HDF5_ATL10(*args, **kwargs):
    """
    Reads ICESat-2 ATL10 (Sea Ice Freeboard) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL10 file
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables

    Returns
    -------
    IS2_atl10_mds: dict
        ATL10 variables
    IS2_atl10_attrs:
        ATL10 attributes
    IS2_atl10_beams: list
        valid ICESat-2 beams within ATL10 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL10.read_granule(*args, **kwargs)

# PURPOSE: find valid beam groups within ICESat-2 ATL10 HDF5 data files
def find_HDF5_ATL10_beams(*args, **kwargs):
    """
    Find valid beam groups within ICESat-2 ATL10 (Sea Ice Freeboard) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL10 file

    Returns
    -------
    IS2_atl10_beams: list
        valid ICESat-2 beams within ATL10 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL10.find_beams(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL10 HDF5 data files for beam variables
def read_HDF5_ATL10_beam(*args, **kwargs):
    """
    Reads ICESat-2 ATL10 (Sea Ice Freeboard) data files for a specific beam

    Parameters
    ----------
    FILENAME: str
        full path to ATL10 file
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
    IS2_atl10_mds: dict
        ATL10 variables
    IS2_atl10_attrs:
        ATL10 attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL10.read_beam(*args, **kwargs)
