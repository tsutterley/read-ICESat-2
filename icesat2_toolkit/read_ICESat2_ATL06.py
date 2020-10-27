#!/usr/bin/env python
u"""
read_ICESat2_ATL06.py (10/2020)
Read ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files

OPTIONS:
    ATTRIBUTES: read HDF5 attributes for groups and variables
    HISTOGRAM: read ATL06 residual_histogram variables
    QUALITY: read ATL06 segment_quality variables
    VERBOSE: output information about input ATL06 file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Updated 10/2020: add small function to find valid beam groups
    Updated 07/2020: added function docstrings
    Updated 11/2019: create attribute dictionaries but don't fill if False
        add function for reading only beam level variables
    Updated 02/2019: continued writing read program with first ATL03 release
    Written 07/2017
"""
from __future__ import print_function

import os
import re
import h5py
import numpy as np

#-- PURPOSE: read ICESat-2 ATL06 HDF5 data files
def read_HDF5_ATL06(FILENAME, ATTRIBUTES=False, HISTOGRAM=False, QUALITY=False,
    VERBOSE=False):
    """
    Reads ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files

    Arguments
    ---------
    FILENAME: full path to ATL06 file

    Keyword arguments
    -----------------
    ATTRIBUTES: read HDF5 attributes for groups and variables
    HISTOGRAM: read ATL06 residual_histogram variables
    QUALITY: read ATL06 segment_quality variables
    VERBOSE: output information about input ATL06 file

    Returns
    -------
    IS2_atl06_mds: dictionary with ATL06 variables
    IS2_atl06_attrs: dictionary with ATL06 attributes
    IS2_atl06_beams: list with valid ICESat-2 beams within ATL06 file
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL06 variables and attributes
    IS2_atl06_mds = {}
    IS2_atl06_attrs = {}

    #-- read each input beam within the file
    IS2_atl06_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        #-- check if subsetted beam contains land ice data
        try:
            fileID[gtx]['land_ice_segments']['segment_id']
        except KeyError:
            pass
        else:
            IS2_atl06_beams.append(gtx)

    #-- read each input beam within the file
    for gtx in IS2_atl06_beams:
        IS2_atl06_mds[gtx] = {}
        IS2_atl06_mds[gtx]['land_ice_segments'] = {}
        IS2_atl06_mds[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl06_mds[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_mds[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl06_mds[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_mds[gtx]['land_ice_segments']['ground_track'] = {}
        #-- get each HDF5 variable
        #-- ICESat-2 land_ice_segments Group
        for key,val in fileID[gtx]['land_ice_segments'].items():
            if isinstance(val, h5py.Dataset):
                IS2_atl06_mds[gtx]['land_ice_segments'][key] = val[:]
            elif isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl06_mds[gtx]['land_ice_segments'][key][k] = v[:]

        #-- ICESat-2 residual_histogram Group
        if HISTOGRAM:
            IS2_atl06_mds[gtx]['residual_histogram'] = {}
            for key,val in fileID[gtx]['residual_histogram'].items():
                IS2_atl06_mds[gtx]['residual_histogram'][key] = val[:]

        #-- ICESat-2 segment_quality Group
        if QUALITY:
            IS2_atl06_mds[gtx]['segment_quality'] = {}
            IS2_atl06_mds[gtx]['segment_quality']['signal_selection_status'] = {}
            for key,val in fileID[gtx]['segment_quality'].items():
                if isinstance(val, h5py.Dataset):
                    IS2_atl06_mds[gtx]['segment_quality'][key] = val[:]
                elif isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl06_mds[gtx]['segment_quality'][key][k] = v[:]

        #-- Getting attributes of included variables
        if ATTRIBUTES:
            #-- Getting attributes of ICESat-2 ATL06 beam variables
            IS2_atl06_attrs[gtx] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments'] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments']['bias_correction'] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments']['dem'] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments']['geophysical'] = {}
            IS2_atl06_attrs[gtx]['land_ice_segments']['ground_track'] = {}
            #-- Global Group Attributes for ATL06 beam
            for att_name,att_val in fileID[gtx].attrs.items():
                IS2_atl06_attrs[gtx][att_name] = att_val
            for key,val in fileID[gtx]['land_ice_segments'].items():
                IS2_atl06_attrs[gtx]['land_ice_segments'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl06_attrs[gtx]['land_ice_segments'][key][k] = {}
                        for att_name,att_val in v.attrs.items():
                            IS2_atl06_attrs[gtx]['land_ice_segments'][key][k][att_name] = att_val
        #-- Getting attributes of histogram variables
        if ATTRIBUTES and HISTOGRAM:
            #-- ICESat-2 residual_histogram Group
            IS2_atl06_attrs[gtx]['residual_histogram'] = {}
            for key,val in fileID[gtx]['residual_histogram'].items():
                IS2_atl06_attrs[gtx]['residual_histogram'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl06_attrs[gtx]['residual_histogram'][key][att_name] = att_val
        #-- Getting attributes of quality variables
        if ATTRIBUTES and QUALITY:
            #-- ICESat-2 segment_quality Group
            IS2_atl06_attrs[gtx]['segment_quality'] = {}
            IS2_atl06_attrs[gtx]['segment_quality']['signal_selection_status'] = {}
            for key,val in fileID[gtx]['segment_quality'].items():
                IS2_atl06_attrs[gtx]['segment_quality'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl06_attrs[gtx]['segment_quality'][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl06_attrs[gtx]['segment_quality'][key][k] = {}
                        for att_name,att_val in v.attrs.items():
                            IS2_atl06_attrs[gtx]['segment_quality'][key][k][att_name] = att_val

    #-- ICESat-2 orbit_info Group
    IS2_atl06_mds['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl06_mds['orbit_info'][key] = val[:]
    #-- ICESat-2 quality_assessment Group
    IS2_atl06_mds['quality_assessment'] = {}
    for key,val in fileID['quality_assessment'].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl06_mds['quality_assessment'][key] = val[:]
        elif isinstance(val, h5py.Group):
            IS2_atl06_mds['quality_assessment'][key] = {}
            for k,v in val.items():
                IS2_atl06_mds['quality_assessment'][key][k] = v[:]

    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl06_mds['ancillary_data'] = {}
    IS2_atl06_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl06_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl06_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl06_attrs['ancillary_data'][key][att_name] = att_val

    #-- land ice ancillary information (first photon bias and statistics)
    cal1,cal2 = ('ancillary_data','land_ice')
    IS2_atl06_mds[cal1][cal2] = {}
    IS2_atl06_attrs[cal1][cal2] = {}
    for key,val in fileID[cal1][cal2].items():
        #-- get each HDF5 variable
        IS2_atl06_mds[cal1][cal2][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl06_attrs[cal1][cal2][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl06_attrs[cal1][cal2][key][att_name] = att_val

    #-- get each global attribute and the attributes for orbit and quality
    if ATTRIBUTES:
        #-- ICESat-2 HDF5 global attributes
        for att_name,att_val in fileID.attrs.items():
            IS2_atl06_attrs[att_name] = att_name
        #-- ICESat-2 orbit_info Group
        IS2_atl06_attrs['orbit_info'] = {}
        for key,val in fileID['orbit_info'].items():
            IS2_atl06_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl06_attrs['orbit_info'][key][att_name]= att_val
        #-- ICESat-2 quality_assessment Group
        IS2_atl06_attrs['quality_assessment'] = {}
        for key,val in fileID['quality_assessment'].items():
            IS2_atl06_attrs['quality_assessment'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl06_attrs['quality_assessment'][key][att_name]= att_val
            if isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl06_attrs['quality_assessment'][key][k] = {}
                    for att_name,att_val in v.attrs.items():
                        IS2_atl06_attrs['quality_assessment'][key][k][att_name]= att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_beams)

#-- PURPOSE: find valid beam groups within ICESat-2 ATL06 HDF5 data files
def find_HDF5_ATL06_beams(FILENAME):
    """
    Find valid beam groups within ICESat-2 ATL06 (Land Ice Along-Track
    Height Product) data files

    Arguments
    ---------
    FILENAME: full path to ATL06 file

    Returns
    -------
    IS2_atl06_beams: list with valid ICESat-2 beams within ATL06 file
    """
    #-- output list of beams
    IS2_atl06_beams = []
    #-- Open the HDF5 file for reading
    with h5py.File(os.path.expanduser(FILENAME), 'r') as fileID:
        #-- read each input beam within the file
        for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
            #-- check if subsetted beam contains land ice data
            try:
                fileID[gtx]['land_ice_segments']['segment_id']
            except KeyError:
                pass
            else:
                IS2_atl06_beams.append(gtx)
    #-- return the list of beams
    return IS2_atl06_beams

#-- PURPOSE: read ICESat-2 ATL06 HDF5 data files for beam variables
def read_HDF5_ATL06_beam(FILENAME, gtx, ATTRIBUTES=False, VERBOSE=False):
    """
    Reads ICESat-2 ATL06 (Land Ice Along-Track Height Product) data files
    for a specific beam

    Arguments
    ---------
    FILENAME: full path to ATL06 file
    gtx: beam name based on ground track and position
        GT1L
        GT1R
        GT2L
        GT2R
        GT3L
        GT3R

    Keyword arguments
    -----------------
    ATTRIBUTES: read HDF5 attributes for groups and variables
    HISTOGRAM: read ATL06 residual_histogram variables
    QUALITY: read ATL06 segment_quality variables
    VERBOSE: output information about input ATL06 file

    Returns
    -------
    IS2_atl06_mds: dictionary with ATL06 variables
    IS2_atl06_attrs: dictionary with ATL06 attributes
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL06 variables and attributes
    IS2_atl06_mds = {}
    IS2_atl06_attrs = {}

    #-- read input beam within the file
    IS2_atl06_mds[gtx] = {}
    IS2_atl06_mds[gtx]['land_ice_segments'] = {}
    IS2_atl06_mds[gtx]['land_ice_segments']['bias_correction'] = {}
    IS2_atl06_mds[gtx]['land_ice_segments']['dem'] = {}
    IS2_atl06_mds[gtx]['land_ice_segments']['fit_statistics'] = {}
    IS2_atl06_mds[gtx]['land_ice_segments']['geophysical'] = {}
    IS2_atl06_mds[gtx]['land_ice_segments']['ground_track'] = {}
    #-- get each HDF5 variable
    #-- ICESat-2 land_ice_segments Group
    for key,val in fileID[gtx]['land_ice_segments'].items():
        if isinstance(val, h5py.Dataset):
            IS2_atl06_mds[gtx]['land_ice_segments'][key] = val[:]
        elif isinstance(val, h5py.Group):
            for k,v in val.items():
                IS2_atl06_mds[gtx]['land_ice_segments'][key][k] = v[:]

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- Getting attributes of ICESat-2 ATL06 beam variables
        IS2_atl06_attrs[gtx] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments'] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_attrs[gtx]['land_ice_segments']['ground_track'] = {}
        #-- Global Group Attributes for ATL06 beam
        for att_name,att_val in fileID[gtx].attrs.items():
            IS2_atl06_attrs[gtx][att_name] = att_val
        for key,val in fileID[gtx]['land_ice_segments'].items():
            IS2_atl06_attrs[gtx]['land_ice_segments'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name] = att_val
            if isinstance(val, h5py.Group):
                for k,v in val.items():
                    IS2_atl06_attrs[gtx]['land_ice_segments'][key][k] = {}
                    for att_name,att_val in v.attrs.items():
                        IS2_atl06_attrs[gtx]['land_ice_segments'][key][k][att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl06_mds,IS2_atl06_attrs)
