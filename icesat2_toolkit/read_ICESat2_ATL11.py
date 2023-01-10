#!/usr/bin/env python
u"""
read_ICESat2_ATL11.py (12/2022)
Read ICESat-2 ATL11 (Annual Land Ice Height) data files

OPTIONS:
    GROUPS: HDF5 groups to read for each beam pair
    ATTRIBUTES: read HDF5 attributes for groups and variables
    REFERENCE: read ATL11 reference surface variables
    CROSSOVERS: read ATL11 crossover height variables
    SUBSETTING: read ATL11 subsetting variables

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
    Updated 03/2021: added function for reading only pair level variables
        simplified group reads to iterate within a try/except statement
        added keyword argument to explicitly define groups to read
    Updated 02/2021: add check if input streaming from bytes
        add small function to find valid beam pair groups
    Written 01/2021
"""
import warnings
import icesat2_toolkit.io

# PURPOSE: read ICESat-2 ATL11 HDF5 data files
def read_HDF5_ATL11(*args, **kwargs):
    """
    Reads ICESat-2 ATL11 (Annual Land Ice Height) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL11 file
    GROUPS: list, default ['cycle_stats']
        HDF5 groups to read for each beam pair
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables
    REFERENCE: bool, default False
        read ATL11 reference surface variables
    CROSSOVERS: bool, default False
        read ATL11 crossover height variables
    SUBSETTING: bool, default False
        read ATL11 subsetting variables

    Returns
    -------
    IS2_atl11_mds:
        ATL11 variables
    IS2_atl11_attrs:
        ATL11 attributes
    IS2_atl11_pairs: list
        valid ICESat-2 beam pairs within ATL11 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL11.read_granule(*args, **kwargs)

# PURPOSE: find valid beam pair groups within ICESat-2 ATL11 HDF5 data files
def find_HDF5_ATL11_pairs(*args, **kwargs):
    """
    Find valid beam pair groups within ICESat-2 ATL11 (Annual Land Ice Height)
    data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL11 file

    Returns
    -------
    IS2_atl11_pairs: list
        valid ICESat-2 beam pairs within ATL11 file
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL11.find_pairs(*args, **kwargs)

# PURPOSE: read ICESat-2 ATL11 HDF5 data files for a specific beam pair
def read_HDF5_ATL11_pair(*args, **kwargs):
    """
    Reads ICESat-2 ATL11 (Annual Land Ice Height) data files
    for a specific beam pair

    Parameters
    ----------
    FILENAME: str
        full path to ATL11 file
    ptx: str
        beam pair name

            - ``'pt1'``
            - ``'pt2'``
            - ``'pt3'``
    GROUPS: list, default ['cycle_stats']
        HDF5 groups to read for each beam pair
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables
    REFERENCE: bool, default False
        read ATL11 reference surface variables
    CROSSOVERS: bool, default False
        read ATL11 crossover height variables
    SUBSETTING: bool, default False
        read ATL11 subsetting variables

    Returns
    -------
    IS2_atl11_mds: dict
        ATL11 variables
    IS2_atl11_attrs: dict
        ATL11 attributes
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use icesat2_toolkit.io instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call renamed version to not break workflows
    return icesat2_toolkit.io.ATL11.read_pair(*args, **kwargs)
