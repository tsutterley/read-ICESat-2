#!/usr/bin/env python
u"""
ATL10.py (05/2024)
Read ICESat-2 ATL10 (Sea Ice Freeboard) data files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 05/2024: use wrapper to importlib for optional dependencies
        check if input filename is an open HDF5 file object
    Updated 03/2024: use pathlib to define and operate on paths
    Updated 11/2023: drop DIMENSION_LIST, CLASS and NAME attributes
    Updated 05/2023: extract more ancillary data from ATL10 files
    Updated 12/2022: place some imports behind try/except statements
        refactor ICESat-2 data product read programs under io
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 12/2021
"""
from __future__ import print_function

import io
import re
import logging
import pathlib
import warnings
import numpy as np
from icesat2_toolkit.utilities import import_dependency

# attempt imports
h5py = import_dependency('h5py')

# PURPOSE: read ICESat-2 ATL10 HDF5 data files
def read_granule(FILENAME, ATTRIBUTES=False, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    elif isinstance(FILENAME, h5py.File):
        fileID = FILENAME
    else:
        FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
        fileID = h5py.File(FILENAME, 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL10 variables and attributes
    IS2_atl10_mds = {}
    IS2_atl10_attrs = {}

    # read each input beam within the file
    IS2_atl10_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains sea ice freeboard data
        try:
            fileID[gtx]['freeboard_beam_segment']['delta_time']
            fileID[gtx]['leads']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl10_beams.append(gtx)

    # read each input beam within the file
    for gtx in IS2_atl10_beams:
        IS2_atl10_mds[gtx] = {}
        IS2_atl10_mds[gtx]['freeboard_beam_segment'] = {}
        IS2_atl10_mds[gtx]['freeboard_beam_segment']['beam_freeboard'] = {}
        IS2_atl10_mds[gtx]['freeboard_beam_segment']['geophysical'] = {}
        IS2_atl10_mds[gtx]['freeboard_beam_segment']['height_segments'] = {}
        IS2_atl10_mds[gtx]['leads'] = {}
        # get each HDF5 variable
        for group in ['freeboard_beam_segment','leads']:
            # ICESat-2 Group
            for key,val in fileID[gtx][group].items():
                if isinstance(val, h5py.Dataset):
                    IS2_atl10_mds[gtx][group][key] = val[:]
                elif isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl10_mds[gtx][group][key][k] = v[:]

        # Getting attributes of included variables
        if ATTRIBUTES:
            # Getting attributes of ICESat-2 ATL10 beam variables
            IS2_atl10_attrs[gtx] = {}
            IS2_atl10_attrs[gtx]['freeboard_beam_segment'] = {}
            IS2_atl10_attrs[gtx]['freeboard_beam_segment']['beam_freeboard'] = {}
            IS2_atl10_attrs[gtx]['freeboard_beam_segment']['geophysical'] = {}
            IS2_atl10_attrs[gtx]['freeboard_beam_segment']['height_segments'] = {}
            IS2_atl10_attrs[gtx]['leads'] = {}
            # Global Group Attributes for ATL10 beam
            for att_name,att_val in fileID[gtx].attrs.items():
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    IS2_atl10_attrs[gtx][att_name] = att_val
            for group in ['freeboard_beam_segment','leads']:
                for key,val in fileID[gtx][group].items():
                    IS2_atl10_attrs[gtx][group][key] = {}
                    for att_name,att_val in val.attrs.items():
                        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                            IS2_atl10_attrs[gtx][group][key][att_name] = att_val
                    if isinstance(val, h5py.Group):
                        for k,v in val.items():
                            IS2_atl10_attrs[gtx][group][key][k] = {}
                            for att_name,att_val in v.attrs.items():
                                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                                    IS2_atl10_attrs[gtx][group][key][k][att_name] = att_val

    # ICESat-2 orbit_info Group
    IS2_atl10_mds['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl10_mds['orbit_info'][key] = val[:]
    # ICESat-2 quality_assessment Group
    IS2_atl10_mds['quality_assessment'] = {}
    for key,val in fileID['quality_assessment'].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl10_mds['quality_assessment'][key] = val[:]
        elif isinstance(val, h5py.Group):
            IS2_atl10_mds['quality_assessment'][key] = {}
            for k,v in val.items():
                IS2_atl10_mds['quality_assessment'][key][k] = v[:]

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    # could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    # and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl10_mds['ancillary_data'] = {}
    IS2_atl10_attrs['ancillary_data'] = {}
    ancillary_keys = ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc',
        'end_cycle','end_geoseg','end_gpssow','end_gpsweek','end_orbit',
        'end_region','end_rgt','granule_end_utc','granule_start_utc','release',
        'start_cycle','start_geoseg','start_gpssow','start_gpsweek',
        'start_orbit','start_region','start_rgt','version']
    for key in ancillary_keys:
        # get each HDF5 variable
        IS2_atl10_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        if ATTRIBUTES:
            # Variable Attributes
            IS2_atl10_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    IS2_atl10_attrs['ancillary_data'][key][att_name] = att_val

    # sea ice ancillary information (processing flags and parameters)
    for cal in ('freeboard_estimation',):
        IS2_atl10_mds['ancillary_data'][cal] = {}
        IS2_atl10_attrs['ancillary_data'][cal] = {}
        for key,val in fileID['ancillary_data'][cal].items():
            # get each HDF5 variable
            IS2_atl10_mds['ancillary_data'][cal][key] = val[:]
            # Getting attributes of group and included variables
            if ATTRIBUTES:
                # Variable Attributes
                IS2_atl10_attrs['ancillary_data'][cal][key] = {}
                for att_name,att_val in val.attrs.items():
                    if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                        IS2_atl10_attrs['ancillary_data'][cal][key][att_name] = att_val

    # get each global attribute and the attributes for orbit and quality
    if ATTRIBUTES:
        # ICESat-2 HDF5 global attributes
        for att_name,att_val in fileID.attrs.items():
            if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                IS2_atl10_attrs[att_name] = att_name
        # ICESat-2 orbit_info Group
        IS2_atl10_attrs['orbit_info'] = {}
        for key,val in fileID['orbit_info'].items():
            IS2_atl10_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    IS2_atl10_attrs['orbit_info'][key][att_name]= att_val
        # ICESat-2 quality_assessment Group
        IS2_atl10_attrs['quality_assessment'] = {}
        for key,val in fileID['quality_assessment'].items():
            IS2_atl10_attrs['quality_assessment'][key] = {}
            for att_name,att_val in val.attrs.items():
                if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                    IS2_atl10_attrs['quality_assessment'][key][att_name]= att_val
            if isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl10_attrs['quality_assessment'][key][k] = {}
                    for att_name,att_val in v.attrs.items():
                        if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                            IS2_atl10_attrs['quality_assessment'][key][k][att_name]= att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl10_mds, IS2_atl10_attrs, IS2_atl10_beams)

# PURPOSE: find valid beam groups within ICESat-2 ATL10 HDF5 data files
def find_beams(FILENAME, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    elif isinstance(FILENAME, h5py.File):
        fileID = FILENAME
    else:
        FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
        fileID = h5py.File(FILENAME, 'r')
    # output list of beams
    IS2_atl10_beams = []
    # read each input beam within the file
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains sea ice data
        try:
            fileID[gtx]['freeboard_beam_segment']['height_segment_id']
            fileID[gtx]['freeboard_beam_segment']['leads']
        except KeyError:
            pass
        else:
            IS2_atl10_beams.append(gtx)
    # Closing the HDF5 file
    fileID.close()
    # return the list of beams
    return IS2_atl10_beams

# PURPOSE: read ICESat-2 ATL10 HDF5 data files for beam variables
def read_beam(FILENAME, gtx, ATTRIBUTES=False, **kwargs):
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
    # Open the HDF5 file for reading
    if isinstance(FILENAME, io.IOBase):
        fileID = h5py.File(FILENAME, 'r')
    elif isinstance(FILENAME, h5py.File):
        fileID = FILENAME
    else:
        FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
        fileID = h5py.File(FILENAME, 'r')

    # Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))

    # allocate python dictionaries for ICESat-2 ATL10 variables and attributes
    IS2_atl10_mds = {}
    IS2_atl10_attrs = {}

    # read input beam within the file
    IS2_atl10_mds[gtx] = {}
    IS2_atl10_mds[gtx]['freeboard_beam_segment'] = {}
    IS2_atl10_mds[gtx]['freeboard_beam_segment']['beam_freeboard'] = {}
    IS2_atl10_mds[gtx]['freeboard_beam_segment']['geophysical'] = {}
    IS2_atl10_mds[gtx]['freeboard_beam_segment']['height_segments'] = {}
    IS2_atl10_mds[gtx]['leads'] = {}
    # get each HDF5 variable
    for group in ['freeboard_beam_segment','leads']:
        # ICESat-2 Group
        for key,val in fileID[gtx][group].items():
            if isinstance(val, h5py.Dataset):
                IS2_atl10_mds[gtx][group][key] = val[:]
            elif isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl10_mds[gtx][group][key][k] = v[:]

    # Getting attributes of included variables
    if ATTRIBUTES:
        # Getting attributes of ICESat-2 ATL10 beam variables
        IS2_atl10_attrs[gtx] = {}
        IS2_atl10_attrs[gtx]['freeboard_beam_segment'] = {}
        IS2_atl10_attrs[gtx]['freeboard_beam_segment']['beam_freeboard'] = {}
        IS2_atl10_attrs[gtx]['freeboard_beam_segment']['geophysical'] = {}
        IS2_atl10_attrs[gtx]['freeboard_beam_segment']['height_segments'] = {}
        IS2_atl10_attrs[gtx]['leads'] = {}
        # Global Group Attributes for ATL10 beam
        for att_name,att_val in fileID[gtx].attrs.items():
            if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                IS2_atl10_attrs[gtx][att_name] = att_val
        for group in ['freeboard_beam_segment','leads']:
            for key,val in fileID[gtx][group].items():
                IS2_atl10_attrs[gtx][group][key] = {}
                for att_name,att_val in val.attrs.items():
                    if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                        IS2_atl10_attrs[gtx][group][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl10_attrs[gtx][group][key][k] = {}
                        for att_name,att_val in v.attrs.items():
                            if att_name not in ('DIMENSION_LIST','CLASS','NAME'):
                                IS2_atl10_attrs[gtx][group][key][k][att_name] = att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl10_mds, IS2_atl10_attrs)
