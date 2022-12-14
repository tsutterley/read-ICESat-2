#!/usr/bin/env python
u"""
ATL12.py (12/2022)
Read ICESat-2 ATL12 (Ocean Surface Height) data files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 12/2022: place some imports behind try/except statements
        refactor ICESat-2 data product read programs under io
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: using python logging for handling verbose output
    Updated 02/2021: add check if input streaming from bytes
    Updated 10/2020: add small function to find valid beam groups
    Updated 07/2020: added function docstrings
    Written 12/2019
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

# PURPOSE: read ICESat-2 ATL12 HDF5 data files
def read_granule(FILENAME, ATTRIBUTES=False, **kwargs):
    """
    Reads ICESat-2 ATL12 (Ocean Surface Height) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL12 file
    ATTRIBUTES: bool, default False
        read HDF5 attributes for groups and variables

    Returns
    -------
    IS2_atl12_mds: dict
        ATL12 variables
    IS2_atl12_attrs: dict
        TL12 attributes
    IS2_atl12_beams: list
        valid ICESat-2 beams within ATL12 file
    """
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL12 variables and attributes
    IS2_atl12_mds = {}
    IS2_atl12_attrs = {}

    # read each input beam within the file
    IS2_atl12_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains ocean surface height data
        try:
            fileID[gtx]['ssh_segments']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl12_beams.append(gtx)

    # read each input beam within the file
    for gtx in IS2_atl12_beams:
        IS2_atl12_mds[gtx] = {}
        IS2_atl12_mds[gtx]['ssh_segments'] = {}
        IS2_atl12_mds[gtx]['ssh_segments']['heights'] = {}
        IS2_atl12_mds[gtx]['ssh_segments']['stats'] = {}
        # get each HDF5 variable
        # ICESat-2 ssh_segments Group
        for key,val in fileID[gtx]['ssh_segments'].items():
            if isinstance(val, h5py.Dataset):
                IS2_atl12_mds[gtx]['ssh_segments'][key] = val[:]
            elif isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl12_mds[gtx]['ssh_segments'][key][k] = v[:]

        # Getting attributes of included variables
        if ATTRIBUTES:
            # Getting attributes of ICESat-2 ATL12 beam variables
            IS2_atl12_attrs[gtx] = {}
            IS2_atl12_attrs[gtx]['ssh_segments'] = {}
            IS2_atl12_attrs[gtx]['ssh_segments']['heights'] = {}
            IS2_atl12_attrs[gtx]['ssh_segments']['stats'] = {}
            # Global Group Attributes for ATL12 beam
            for att_name,att_val in fileID[gtx].attrs.items():
                IS2_atl12_attrs[gtx][att_name] = att_val
            for key,val in fileID[gtx]['ssh_segments'].items():
                IS2_atl12_attrs[gtx]['ssh_segments'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl12_attrs[gtx]['ssh_segments'][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl12_attrs[gtx]['ssh_segments'][key][k] = {}
                        for att_name,att_val in v.attrs.items():
                            IS2_atl12_attrs[gtx]['ssh_segments'][key][k][att_name] = att_val

    # ICESat-2 orbit_info Group
    IS2_atl12_mds['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl12_mds['orbit_info'][key] = val[:]
    # ICESat-2 quality_assessment Group
    IS2_atl12_mds['quality_assessment'] = {}
    for key,val in fileID['quality_assessment'].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl12_mds['quality_assessment'][key] = val[:]
        elif isinstance(val, h5py.Group):
            IS2_atl12_mds['quality_assessment'][key] = {}
            for k,v in val.items():
                IS2_atl12_mds['quality_assessment'][key][k] = v[:]

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    # could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    # and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl12_mds['ancillary_data'] = {}
    IS2_atl12_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl12_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        if ATTRIBUTES:
            # Variable Attributes
            IS2_atl12_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl12_attrs['ancillary_data'][key][att_name] = att_val

    # ocean height ancillary information (processing flags and parameters)
    IS2_atl12_mds['ancillary_data']['ocean'] = {}
    IS2_atl12_attrs['ancillary_data']['ocean'] = {}
    for key,val in fileID['ancillary_data']['ocean'].items():
        # get each HDF5 variable
        IS2_atl12_mds['ancillary_data']['ocean'][key] = val[:]
        # Getting attributes of group and included variables
        if ATTRIBUTES:
            # Variable Attributes
            IS2_atl12_attrs['ancillary_data']['ocean'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl12_attrs['ancillary_data']['ocean'][key][att_name] = att_val

    # get each global attribute and the attributes for orbit and quality
    if ATTRIBUTES:
        # ICESat-2 HDF5 global attributes
        for att_name,att_val in fileID.attrs.items():
            IS2_atl12_attrs[att_name] = att_name
        # ICESat-2 orbit_info Group
        IS2_atl12_attrs['orbit_info'] = {}
        for key,val in fileID['orbit_info'].items():
            IS2_atl12_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl12_attrs['orbit_info'][key][att_name]= att_val
        # ICESat-2 quality_assessment Group
        IS2_atl12_attrs['quality_assessment'] = {}
        for key,val in fileID['quality_assessment'].items():
            IS2_atl12_attrs['quality_assessment'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl12_attrs['quality_assessment'][key][att_name]= att_val
            if isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl12_attrs['quality_assessment'][key][k] = {}
                    for att_name,att_val in v.attrs.items():
                        IS2_atl12_attrs['quality_assessment'][key][k][att_name]= att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl12_mds,IS2_atl12_attrs,IS2_atl12_beams)

# PURPOSE: find valid beam groups within ICESat-2 ATL12 HDF5 data files
def find_beams(FILENAME, **kwargs):
    """
    Find valid beam groups within ICESat-2 ATL12 (Ocean Surface Height) data files

    Parameters
    ----------
    FILENAME: str
        full path to ATL12 file

    Returns
    -------
    IS2_atl12_beams: list
        valid ICESat-2 beams within ATL12 file
    """
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')
    # output list of beams
    IS2_atl12_beams = []
    # read each input beam within the file
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains ocean surface height data
        try:
            fileID[gtx]['ssh_segments']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl12_beams.append(gtx)
    # Closing the HDF5 file
    fileID.close()
    # return the list of beams
    return IS2_atl12_beams

# PURPOSE: read ICESat-2 ATL12 HDF5 data files for beam variables
def read_beam(FILENAME, gtx, ATTRIBUTES=False, **kwargs):
    """
    Reads ICESat-2 ATL12 (Ocean Surface Height) data files for a specific beam

    Parameters
    ----------
    FILENAME: str
        full path to ATL12 file
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
    IS2_atl12_mds: dict
        ATL12 variables
    IS2_atl12_attrs:
        ATL12 attributes
    """
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL12 variables and attributes
    IS2_atl12_mds = {}
    IS2_atl12_attrs = {}

    # read input beam within the file
    IS2_atl12_mds[gtx] = {}
    IS2_atl12_mds[gtx]['ssh_segments'] = {}
    IS2_atl12_mds[gtx]['ssh_segments']['heights'] = {}
    IS2_atl12_mds[gtx]['ssh_segments']['stats'] = {}
    # get each HDF5 variable
    # ICESat-2 ssh_segments Group
    for key,val in fileID[gtx]['ssh_segments'].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl12_mds[gtx]['ssh_segments'][key] = val[:]
        elif isinstance(val, h5py.Group):
            for k,v in val.items():
                IS2_atl12_mds[gtx]['ssh_segments'][key][k] = v[:]

    # Getting attributes of included variables
    if ATTRIBUTES:
        # Getting attributes of ICESat-2 ATL12 beam variables
        IS2_atl12_attrs[gtx] = {}
        IS2_atl12_attrs[gtx]['ssh_segments'] = {}
        IS2_atl12_attrs[gtx]['ssh_segments']['heights'] = {}
        IS2_atl12_attrs[gtx]['ssh_segments']['stats'] = {}
        # Global Group Attributes for ATL12 beam
        for att_name,att_val in fileID[gtx].attrs.items():
            IS2_atl12_attrs[gtx][att_name] = att_val
        for key,val in fileID[gtx]['ssh_segments'].items():
            IS2_atl12_attrs[gtx]['ssh_segments'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl12_attrs[gtx]['ssh_segments'][key][att_name] = att_val
            if isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl12_attrs[gtx]['ssh_segments'][key][k] = {}
                    for att_name,att_val in v.attrs.items():
                        IS2_atl12_attrs[gtx]['ssh_segments'][key][k][att_name] = att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl12_mds,IS2_atl12_attrs)
