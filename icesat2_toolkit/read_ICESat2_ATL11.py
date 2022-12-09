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
from __future__ import print_function

import os
import io
import re
import logging
import warnings
import numpy as np

# attempt imports
try:
    import h5py
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: read ICESat-2 ATL11 HDF5 data files
def read_HDF5_ATL11(FILENAME, GROUPS=['cycle_stats'], ATTRIBUTES=False,
    REFERENCE=False, CROSSOVERS=False, SUBSETTING=False, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL11 variables and attributes
    IS2_atl11_mds = {}
    IS2_atl11_attrs = {}

    # read each input beam pair within the file
    IS2_atl11_pairs = []
    for ptx in [k for k in fileID.keys() if bool(re.match(r'pt\d',k))]:
        # check if subsetted beam contains reference points
        try:
            fileID[ptx]['ref_pt']
        except KeyError:
            pass
        else:
            IS2_atl11_pairs.append(ptx)

    # groups to read from ATL11 file
    # ATL11 ref_surf group
    if REFERENCE:
        GROUPS.append('ref_surf')
    # ATL11 crossing_track_data group
    if CROSSOVERS:
        GROUPS.append('crossing_track_data')
    # ATL11 subsetting group
    if SUBSETTING:
        GROUPS.append('subsetting')

    # read each input pair track within the file
    for ptx in IS2_atl11_pairs:
        IS2_atl11_mds[ptx] = {}
        # get each main level HDF5 variable
        for key,val in fileID[ptx].items():
            if isinstance(val, h5py.Dataset):
                IS2_atl11_mds[ptx][key] = val[:]

        # get each HDF5 variable within named groups
        for group in GROUPS:
            try:
                IS2_atl11_mds[ptx][group] = {}
                for key,val in fileID[ptx][group].items():
                    IS2_atl11_mds[ptx][group][key] = val[:]
            except:
                pass

        # Getting attributes of included variables
        if ATTRIBUTES:
            # Getting attributes of ICESat-2 ATL11 main level variables
            IS2_atl11_attrs[ptx] = {}
            # Global Group Attributes for ATL11 beam pair
            for att_name,att_val in fileID[ptx].attrs.items():
                IS2_atl11_attrs[ptx][att_name] = att_val
            # getting attributes of main level ATL11 variables
            for key,val in fileID[ptx].items():
                IS2_atl11_attrs[ptx][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl11_attrs[ptx][key][att_name] = att_val
                # get fill value attributes if applicable
                if hasattr(val,'fillvalue'):
                    IS2_atl11_attrs[ptx][key]['_FillValue'] = \
                        getattr(val,'fillvalue')
            # getting attributes of variables within named groups
            for group in GROUPS:
                try:
                    IS2_atl11_attrs[ptx][group] = {}
                    for key,val in fileID[ptx][group].items():
                        IS2_atl11_attrs[ptx][group][key] = {}
                        for att_name,att_val in val.attrs.items():
                            IS2_atl11_attrs[ptx][group][key][att_name] = att_val
                        # get fill value attributes if applicable
                        if hasattr(val,'fillvalue'):
                            IS2_atl11_attrs[ptx][group][key]['_FillValue'] = \
                                getattr(val,'fillvalue')
                except:
                    pass

    # ICESat-2 orbit_info Group
    IS2_atl11_mds['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl11_mds['orbit_info'][key] = val[:]

    # ICESat-2 quality_assessment Group
    IS2_atl11_mds['quality_assessment'] = {}
    for key,val in fileID['quality_assessment'].items():
        IS2_atl11_mds['quality_assessment'][key] = val[:]

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    # could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    # and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl11_mds['ancillary_data'] = {}
    IS2_atl11_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl11_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        if ATTRIBUTES:
            # Variable Attributes
            IS2_atl11_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl11_attrs['ancillary_data'][key][att_name] = att_val

    # get each global attribute and the attributes for orbit and quality
    if ATTRIBUTES:
        # ICESat-2 HDF5 global attributes
        for att_name,att_val in fileID.attrs.items():
            IS2_atl11_attrs[att_name] = att_name
        # ICESat-2 orbit_info Group
        IS2_atl11_attrs['orbit_info'] = {}
        for key,val in fileID['orbit_info'].items():
            IS2_atl11_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl11_attrs['orbit_info'][key][att_name]= att_val
        # ICESat-2 quality_assessment Group
        IS2_atl11_attrs['quality_assessment'] = {}
        for key,val in fileID['quality_assessment'].items():
            IS2_atl11_attrs['quality_assessment'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl11_attrs['quality_assessment'][key][att_name]= att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs)

# PURPOSE: find valid beam pair groups within ICESat-2 ATL11 HDF5 data files
def find_HDF5_ATL11_pairs(FILENAME, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')
    # read each input beam pair within the file
    IS2_atl11_pairs = []
    for ptx in [k for k in fileID.keys() if bool(re.match(r'pt\d',k))]:
        # check if subsetted beam contains reference points
        try:
            fileID[ptx]['ref_pt']
        except KeyError:
            pass
        else:
            IS2_atl11_pairs.append(ptx)
    # Closing the HDF5 file
    fileID.close()
    # return the list of beam pairs
    return IS2_atl11_pairs

# PURPOSE: read ICESat-2 ATL11 HDF5 data files for a specific beam pair
def read_HDF5_ATL11_pair(FILENAME, ptx, GROUPS=['cycle_stats'],
    ATTRIBUTES=False, REFERENCE=False, CROSSOVERS=False,
    SUBSETTING=False, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL11 variables and attributes
    IS2_atl11_mds = {}
    IS2_atl11_attrs = {}

    # groups to read from ATL11 file
    # ATL11 ref_surf group
    if REFERENCE:
        GROUPS.append('ref_surf')
    # ATL11 crossing_track_data group
    if CROSSOVERS:
        GROUPS.append('crossing_track_data')
    # ATL11 subsetting group
    if SUBSETTING:
        GROUPS.append('subsetting')

    # read input pair track within the file
    IS2_atl11_mds[ptx] = {}
    # get each main level HDF5 variable
    for key,val in fileID[ptx].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl11_mds[ptx][key] = val[:]

    # get each cycle_stats HDF5 variable
    for group in GROUPS:
        try:
            IS2_atl11_mds[ptx][group] = {}
            for key,val in fileID[ptx][group].items():
                IS2_atl11_mds[ptx][group][key] = val[:]
        except:
            pass

    # Getting attributes of included variables
    if ATTRIBUTES:
        # Getting attributes of ICESat-2 ATL11 main level variables
        IS2_atl11_attrs[ptx] = {}
        # Global Group Attributes for ATL11 beam pair
        for att_name,att_val in fileID[ptx].attrs.items():
            IS2_atl11_attrs[ptx][att_name] = att_val
        # getting attributes of main level ATL11 variables
        for key,val in fileID[ptx].items():
            IS2_atl11_attrs[ptx][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl11_attrs[ptx][key][att_name] = att_val
            # get fill value attributes if applicable
            if hasattr(val,'fillvalue'):
                IS2_atl11_attrs[ptx][key]['_FillValue'] = \
                    getattr(val,'fillvalue')
        # getting attributes of variables within named groups
        for group in GROUPS:
            try:
                IS2_atl11_attrs[ptx][group] = {}
                for key,val in fileID[ptx][group].items():
                    IS2_atl11_attrs[ptx][group][key] = {}
                    for att_name,att_val in val.attrs.items():
                        IS2_atl11_attrs[ptx][group][key][att_name] = att_val
                        # get fill value attributes if applicable
                        if hasattr(val,'fillvalue'):
                            IS2_atl11_attrs[ptx][group][key]['_FillValue'] = \
                                getattr(val,'fillvalue')
            except:
                pass

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl11_mds,IS2_atl11_attrs)
