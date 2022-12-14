#!/usr/bin/env python
u"""
read_ICESat2_ATL03.py (12/2022)
Read ICESat-2 ATL03 and ATL09 data files to calculate average segment surfaces
    ATL03 datasets: Global Geolocated Photons
    ATL09 datasets: Atmospheric Characteristics

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 12/2022: place some imports behind try/except statements
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: add check if input streaming from bytes
    Updated 10/2020: add small function to find valid beam groups
    Updated 09/2020: map ATL09 to ATL03 using delta times
    Updated 07/2020: added function docstrings
    Updated 06/2020: add additional beam check within heights groups
    Updated 11/2019: create attribute dictionaries but don't fill if False
    Updated 09/2019: add functions for reading main and beam level variables
    Updated 03/2019: extract a set of ATL09 parameters for each ATL03 segment_ID
    Updated 02/2019: continued writing read program with first ATL03 release
    Written 05/2017
"""
import warnings
import icesat2_toolkit.io

# PURPOSE: read ICESat-2 ATL03 HDF5 data files
def read_HDF5_ATL03(*args, **kwargs):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL03 file
    ATTRIBUTES: bool, default False
        read file, group and variable attributes

    Returns
    -------
    IS2_atl03_mds: dict
        ATL03 variables
    IS2_atl03_attrs: dict
        ATL03 attributes
    IS2_atl03_beams: list
        valid ICESat-2 beams within ATL03 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL03.read_granule(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL09 HDF5 data file for specific variables
def read_HDF5_ATL09(*args, **kwargs):
    """
    Reads ICESat-2 ATL09 Atmospheric Characteristics data files
    and interpolates a subset of variables to ATL03 segment lengths

    Parameters
    ----------
    FILENAME: str
        full path to ATL03 file
    pfl: str
        profile for a given beam
    dtime: float
        ATL03 reference photon delta_time
    ATTRIBUTES: bool, default False
        read file, group and variable attributes

    Returns
    -------
    IS2_atl09_mds: dict
        ATL09 variables remapped to ATL03 segments
    IS2_atl09_attrs: dict
        ATL09 attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL03.interpolate_ATL09(*args, **kwargs)

# PURPOSE: find valid beam groups within ICESat-2 ATL03 HDF5 data files
def find_HDF5_ATL03_beams(*args, **kwargs):
    """
    Find valid beam groups within ICESat-2 ATL03 Global Geolocated Photons
    data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL03 file

    Returns
    -------
    IS2_atl03_beams: list
        valid ICESat-2 beams within ATL03 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL03.find_beams(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL03 HDF5 data files for main level variables
def read_HDF5_ATL03_main(*args, **kwargs):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files
    for only the main-level variables and not the beam-level data

    Parameters
    ----------
    FILENAME: str
        full path to ATL03 file
    ATTRIBUTES: bool, default False
        read file, group and variable attributes

    Returns
    -------
    IS2_atl03_mds: dict
        ATL03 main-level variables
    IS2_atl03_attrs: dict
        ATL03 main-level attributes
    IS2_atl03_beams: list
        valid ICESat-2 beams within ATL03 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL03.read_main(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL03 HDF5 data files for beam variables
def read_HDF5_ATL03_beam(*args, **kwargs):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files
    for a specific beam

    Parameters
    ----------
    FILENAME: str
        full path to ATL03 file
    gtx: str
        beam name based on ground track and position

            - ``'gt1l'``
            - ``'gt1r'``
            - ``'gt2l'``
            - ``'gt2r'``
            - ``'gt3l'``
            - ``'gt3r'``
    ATTRIBUTES: bool, default False
        read file, group and variable attributes

    Returns
    -------
    IS2_atl03_mds: dict
        ATL03 beam-level variables
    IS2_atl03_attrs: dict
        ATL03 beam-level attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL03.read_beam(*args, **kwargs)
