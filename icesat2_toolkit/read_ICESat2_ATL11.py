#!/usr/bin/env python
u"""
read_ICESat2_ATL11.py (01/2021)
Read ICESat-2 ATL11 (Annual Land Ice Height) data files

OPTIONS:
    ATTRIBUTES: read HDF5 attributes for groups and variables
    REFERENCE: read ATL11 reference surface variables
    CROSSOVERS: read ATL11 crossover height variables
    VERBOSE: output information about input ATL11 file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/

UPDATE HISTORY:
    Written 01/2021
"""
from __future__ import print_function

import os
import re
import h5py
import numpy as np

#-- PURPOSE: read ICESat-2 ATL11 HDF5 data files
def read_HDF5_ATL11(FILENAME, ATTRIBUTES=False, REFERENCE=False,
    CROSSOVERS=False, VERBOSE=False):
    """
    Reads ICESat-2 ATL11 (Annual Land Ice Height) data files

    Arguments
    ---------
    FILENAME: full path to ATL11 file

    Keyword arguments
    -----------------
    ATTRIBUTES: read HDF5 attributes for groups and variables
    REFERENCE: read ATL11 reference surface variables
    CROSSOVERS: read ATL11 crossover height variables
    VERBOSE: output information about input ATL11 file

    Returns
    -------
    IS2_atl11_mds: dictionary with ATL11 variables
    IS2_atl11_attrs: dictionary with ATL11 attributes
    IS2_atl11_pairs: list with valid ICESat-2 beam pairs within ATL11 file
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL11 variables and attributes
    IS2_atl11_mds = {}
    IS2_atl11_attrs = {}

    #-- read each input beam pair within the file
    IS2_atl11_pairs = []
    for ptx in [k for k in fileID.keys() if bool(re.match(r'pt\d',k))]:
        #-- check if subsetted beam contains reference points
        try:
            fileID[ptx]['ref_pt']
        except KeyError:
            pass
        else:
            IS2_atl11_pairs.append(ptx)

    #-- read each input pair track within the file
    for ptx in IS2_atl11_pairs:
        IS2_atl11_mds[ptx] = {}
        #-- get each main level HDF5 variable
        for key,val in fileID[ptx].items():
            if isinstance(val, h5py.Dataset):
                IS2_atl11_mds[ptx][key] = val[:]

        #-- get each cycle_stats HDF5 variable
        IS2_atl11_mds[ptx]['cycle_stats'] = {}
        for key,val in fileID[ptx]['cycle_stats'].items():
            IS2_atl11_mds[ptx]['cycle_stats'][key] = val[:]

        #-- ICESat-2 ref_surf Group
        if REFERENCE:
            IS2_atl11_mds[ptx]['ref_surf'] = {}
            for key,val in fileID[ptx]['ref_surf'].items():
                IS2_atl11_mds[ptx]['ref_surf'][key] = val[:]

        #-- ICESat-2 crossing_track_data Group
        if CROSSOVERS:
            IS2_atl11_mds[ptx]['crossing_track_data'] = {}
            for key,val in fileID[ptx]['crossing_track_data'].items():
                IS2_atl11_mds[ptx]['crossing_track_data'][key] = val[:]

        #-- Getting attributes of included variables
        if ATTRIBUTES:
            #-- Getting attributes of ICESat-2 ATL11 main level variables
            IS2_atl11_attrs[ptx] = {}
            #-- Global Group Attributes for ATL11 beam pair
            for att_name,att_val in fileID[ptx].attrs.items():
                IS2_atl11_attrs[ptx][att_name] = att_val
            for key,val in fileID[ptx].items():
                IS2_atl11_attrs[ptx][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl11_attrs[ptx][key][att_name] = att_val
            #-- Getting attributes of ICESat-2 cycle_stats variables
            IS2_atl11_attrs[ptx]['cycle_stats'] = {}
            for key,val in fileID[ptx]['cycle_stats'].items():
                IS2_atl11_attrs[ptx]['cycle_stats'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl11_attrs[ptx]['cycle_stats'][key][att_name] = att_val

        #-- Getting attributes of ref_surf variables
        if ATTRIBUTES and REFERENCE:
            #-- ICESat-2 ref_surf Group
            IS2_atl11_attrs[ptx]['ref_surf'] = {}
            for key,val in fileID[ptx]['ref_surf'].items():
                IS2_atl11_attrs[ptx]['ref_surf'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl11_attrs[ptx]['ref_surf'][key][att_name] = att_val

        #-- Getting attributes of crossover variables
        if ATTRIBUTES and CROSSOVERS:
            #-- ICESat-2 crossing_track_data Group
            IS2_atl11_attrs[ptx]['crossing_track_data'] = {}
            for key,val in fileID[ptx]['crossing_track_data'].items():
                IS2_atl11_attrs[ptx]['crossing_track_data'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl11_attrs[ptx]['crossing_track_data'][key][att_name] = att_val


    #-- ICESat-2 orbit_info Group
    IS2_atl11_mds['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl11_mds['orbit_info'][key] = val[:]

    #-- ICESat-2 quality_assessment Group
    IS2_atl11_mds['quality_assessment'] = {}
    for key,val in fileID['quality_assessment'].items():
        IS2_atl11_mds['quality_assessment'][key] = val[:]

    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl11_mds['ancillary_data'] = {}
    IS2_atl11_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl11_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl11_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl11_attrs['ancillary_data'][key][att_name] = att_val

    #-- get each global attribute and the attributes for orbit and quality
    if ATTRIBUTES:
        #-- ICESat-2 HDF5 global attributes
        for att_name,att_val in fileID.attrs.items():
            IS2_atl11_attrs[att_name] = att_name
        #-- ICESat-2 orbit_info Group
        IS2_atl11_attrs['orbit_info'] = {}
        for key,val in fileID['orbit_info'].items():
            IS2_atl11_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl11_attrs['orbit_info'][key][att_name]= att_val
        #-- ICESat-2 quality_assessment Group
        IS2_atl11_attrs['quality_assessment'] = {}
        for key,val in fileID['quality_assessment'].items():
            IS2_atl11_attrs['quality_assessment'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl11_attrs['quality_assessment'][key][att_name]= att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs)
