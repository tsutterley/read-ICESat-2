#!/usr/bin/env python
u"""
read_ICESat2_ATL03.py (10/2020)
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
from __future__ import print_function, division

import os
import re
import h5py
import numpy as np
import scipy.interpolate

#-- PURPOSE: read ICESat-2 ATL03 HDF5 data files
def read_HDF5_ATL03(FILENAME, ATTRIBUTES=False, VERBOSE=False):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files

    Arguments
    ---------
    FILENAME: full path to ATL03 file

    Keyword arguments
    -----------------
    ATTRIBUTES: read file, group and variable attributes
    VERBOSE: output information about top-level HDF5 groups

    Returns
    -------
    IS2_atl03_mds: dictionary with ATL03 variables
    IS2_atl03_attrs: dictionary with ATL03 attributes
    IS2_atl03_beams: list with valid ICESat-2 beams within ATL03 file
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL03 variables and attributes
    IS2_atl03_mds = {}
    IS2_atl03_attrs = {}

    #-- read each input beam within the file
    IS2_atl03_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        #-- check if subsetted beam contains data
        #-- check in both the geolocation and heights groups
        try:
            fileID[gtx]['geolocation']['segment_id']
            fileID[gtx]['heights']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl03_beams.append(gtx)

    #-- for each included beam
    for gtx in IS2_atl03_beams:
        #-- get each HDF5 variable
        IS2_atl03_mds[gtx] = {}
        IS2_atl03_mds[gtx]['heights'] = {}
        IS2_atl03_mds[gtx]['geolocation'] = {}
        IS2_atl03_mds[gtx]['bckgrd_atlas'] = {}
        IS2_atl03_mds[gtx]['geophys_corr'] = {}
        #-- ICESat-2 Measurement Group
        for key,val in fileID[gtx]['heights'].items():
            IS2_atl03_mds[gtx]['heights'][key] = val[:]
        #-- ICESat-2 Geolocation Group
        for key,val in fileID[gtx]['geolocation'].items():
            IS2_atl03_mds[gtx]['geolocation'][key] = val[:]
        #-- ICESat-2 Background Photon Rate Group
        for key,val in fileID[gtx]['bckgrd_atlas'].items():
            IS2_atl03_mds[gtx]['bckgrd_atlas'][key] = val[:]
        #-- ICESat-2 Geophysical Corrections Group: Values for tides (ocean,
        #-- solid earth, pole, load, and equilibrium), inverted barometer (IB)
        #-- effects, and range corrections for tropospheric delays
        for key,val in fileID[gtx]['geophys_corr'].items():
            IS2_atl03_mds[gtx]['geophys_corr'][key] = val[:]

        #-- Getting attributes of included variables
        if ATTRIBUTES:
            #-- Getting attributes of IS2_atl03_mds beam variables
            IS2_atl03_attrs[gtx] = {}
            IS2_atl03_attrs[gtx]['heights'] = {}
            IS2_atl03_attrs[gtx]['geolocation'] = {}
            IS2_atl03_attrs[gtx]['bckgrd_atlas'] = {}
            IS2_atl03_attrs[gtx]['geophys_corr'] = {}
            #-- Global Group Attributes
            for att_name,att_val in fileID[gtx].attrs.items():
                IS2_atl03_attrs[gtx][att_name] = att_val
            #-- ICESat-2 Measurement Group
            for key,val in fileID[gtx]['heights'].items():
                IS2_atl03_attrs[gtx]['heights'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['heights'][key][att_name]=att_val
            #-- ICESat-2 Geolocation Group
            for key,val in fileID[gtx]['geolocation'].items():
                IS2_atl03_attrs[gtx]['geolocation'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['geolocation'][key][att_name]=att_val
            #-- ICESat-2 Background Photon Rate Group
            for key,val in fileID[gtx]['bckgrd_atlas'].items():
                IS2_atl03_attrs[gtx]['bckgrd_atlas'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['bckgrd_atlas'][key][att_name]=att_val
            #-- ICESat-2 Geophysical Corrections Group
            for key,val in fileID[gtx]['geophys_corr'].items():
                IS2_atl03_attrs[gtx]['geophys_corr'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['geophys_corr'][key][att_name]=att_val

    #-- ICESat-2 spacecraft orientation at time
    IS2_atl03_mds['orbit_info'] = {}
    IS2_atl03_attrs['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl03_mds['orbit_info'][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Global Group Attributes
            for att_name,att_val in fileID['orbit_info'].attrs.items():
                IS2_atl03_attrs['orbit_info'][att_name] = att_val
            #-- Variable Attributes
            IS2_atl03_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['orbit_info'][key][att_name] = att_val

    #-- information ancillary to the data product
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl03_mds['ancillary_data'] = {}
    IS2_atl03_attrs['ancillary_data'] = {}
    ancillary_keys = ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc',
        'end_cycle','end_geoseg','end_gpssow','end_gpsweek','end_orbit',
        'end_region','end_rgt','granule_end_utc','granule_start_utc','release',
        'start_cycle','start_geoseg','start_gpssow','start_gpsweek',
        'start_orbit','start_region','start_rgt','version']
    for key in ancillary_keys:
        #-- get each HDF5 variable
        IS2_atl03_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl03_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl03_attrs['ancillary_data'][key][att_name] = att_val

    #-- transmit-echo-path (tep) parameters
    IS2_atl03_mds['ancillary_data']['tep'] = {}
    IS2_atl03_attrs['ancillary_data']['tep'] = {}
    for key,val in fileID['ancillary_data']['tep'].items():
        #-- get each HDF5 variable
        IS2_atl03_mds['ancillary_data']['tep'][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl03_attrs['ancillary_data']['tep'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['ancillary_data']['tep'][key][att_name] = att_val

    #-- channel dead time and first photon bias derived from ATLAS calibration
    cal1,cal2 = ('ancillary_data','calibrations')
    for var in ['dead_time','first_photon_bias']:
        IS2_atl03_mds[cal1][var] = {}
        IS2_atl03_attrs[cal1][var] = {}
        for key,val in fileID[cal1][cal2][var].items():
            #-- get each HDF5 variable
            if isinstance(val, h5py.Dataset):
                IS2_atl03_mds[cal1][var][key] = val[:]
            elif isinstance(val, h5py.Group):
                IS2_atl03_mds[cal1][var][key] = {}
                for k,v in val.items():
                    IS2_atl03_mds[cal1][var][key][k] = v[:]
            #-- Getting attributes of group and included variables
            if ATTRIBUTES:
                #-- Variable Attributes
                IS2_atl03_attrs[cal1][var][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[cal1][var][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl03_attrs[cal1][var][key][k] = {}
                        for att_name,att_val in val.attrs.items():
                            IS2_atl03_attrs[cal1][var][key][k][att_name]=att_val

    #-- get ATLAS impulse response variables for the transmitter echo path (TEP)
    tep1,tep2 = ('atlas_impulse_response','tep_histogram')
    IS2_atl03_mds[tep1] = {}
    IS2_atl03_attrs[tep1] = {}
    for pce in ['pce1_spot1','pce2_spot3']:
        IS2_atl03_mds[tep1][pce] = {tep2:{}}
        IS2_atl03_attrs[tep1][pce] = {tep2:{}}
        #-- for each TEP variable
        for key,val in fileID[tep1][pce][tep2].items():
            IS2_atl03_mds[tep1][pce][tep2][key] = val[:]
            #-- Getting attributes of included variables
            if ATTRIBUTES:
                #-- Global Group Attributes
                for att_name,att_val in fileID[tep1][pce][tep2].attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][att_name] = att_val
                #-- Variable Attributes
                IS2_atl03_attrs[tep1][pce][tep2][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][key][att_name] = att_val

    #-- Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl03_attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams)

#-- PURPOSE: read ICESat-2 ATL09 HDF5 data file for specific variables
def read_HDF5_ATL09(FILENAME, pfl, dtime, ATTRIBUTES=True):
    """
    Reads ICESat-2 ATL09 Atmospheric Characteristics data files

    Arguments
    ---------
    FILENAME: full path to ATL03 file
    pfl: profile for a given beam
    dtime: ATL03 reference photon delta_time

    Keyword arguments
    -----------------
    ATTRIBUTES: read file, group and variable attributes
    VERBOSE: output information about top-level HDF5 groups

    Returns
    -------
    IS2_atl09_mds: dictionary with ATL09 variables
    IS2_atl09_attrs: dictionary with ATL09 attributes
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- allocate python dictionaries for ICESat-2 ATL09 variables and attributes
    IS2_atl09_mds = {}
    IS2_atl09_attrs = {}

    #-- read profile reported for the ATLAS strong beams within the file
    IS2_atl09_mds[pfl] = dict(high_rate={})
    #-- extract delta_time for mapping ATL09 atmospheric parameters to ATL03
    delta_time = fileID[pfl]['high_rate']['delta_time'][:]
    #-- Calibrated Attenuated Backscatter at 25 hz
    high_rate_keys = ['aclr_true','bsnow_con','bsnow_dens','bsnow_h',
        'bsnow_h_dens','bsnow_od','bsnow_psc','cloud_flag_asr','cloud_flag_atm',
        'cloud_fold_flag','column_od_asr','column_od_asr_qf','msw_flag',
        'snow_ice','solar_azimuth','solar_elevation','surf_refl_true']
    #-- extract variables of interest and map to ATL03 segments
    for key in high_rate_keys:
        val = np.copy(fileID[pfl]['high_rate'][key][:])
        fint = scipy.interpolate.interp1d(delta_time, val,
            kind='nearest', fill_value='extrapolate')
        IS2_atl09_mds[pfl]['high_rate'][key] = fint(dtime).astype(val.dtype)

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- Getting attributes of IS2_atl09_mds profile variables
        IS2_atl09_attrs[pfl] = dict(high_rate={})
        #-- Global Group Attributes
        for att_name,att_val in fileID[pfl].attrs.items():
            IS2_atl09_attrs[pfl][att_name] = att_val
        #-- Variable Attributes
        for key in high_rate_keys:
            IS2_atl09_attrs[pfl]['high_rate'][key] = {}
            for att_name,att_val in fileID[pfl]['high_rate'][key].attrs.items():
                IS2_atl09_attrs[pfl]['high_rate'][key][att_name] = att_val

    #-- Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl09_attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl09_mds,IS2_atl09_attrs)

#-- PURPOSE: find valid beam groups within ICESat-2 ATL03 HDF5 data files
def find_HDF5_ATL03_beams(FILENAME):
    """
    Find valid beam groups within ICESat-2 ATL03 Global Geolocated Photons
    data files

    Arguments
    ---------
    FILENAME: full path to ATL03 file

    Returns
    -------
    IS2_atl03_beams: list with valid ICESat-2 beams within ATL03 file
    """
    #-- output list of beams
    IS2_atl03_beams = []
    #-- Open the HDF5 file for reading
    with h5py.File(os.path.expanduser(FILENAME), 'r') as fileID:
        #-- read each input beam within the file
        for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
            #-- check if subsetted beam contains data
            #-- check in both the geolocation and heights groups
            try:
                fileID[gtx]['geolocation']['segment_id']
                fileID[gtx]['heights']['delta_time']
            except KeyError:
                pass
            else:
                IS2_atl03_beams.append(gtx)
    #-- return the list of beams
    return IS2_atl03_beams

#-- PURPOSE: read ICESat-2 ATL03 HDF5 data files for main level variables
def read_HDF5_ATL03_main(FILENAME, ATTRIBUTES=False, VERBOSE=False):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files
    for only the main-level variables and not the beam-level data

    Arguments
    ---------
    FILENAME: full path to ATL03 file

    Keyword arguments
    -----------------
    ATTRIBUTES: read file, group and variable attributes
    VERBOSE: output information about top-level HDF5 groups

    Returns
    -------
    IS2_atl03_mds: dictionary with ATL03 main-level variables
    IS2_atl03_attrs: dictionary with ATL03 main-level attributes
    IS2_atl03_beams: list with valid ICESat-2 beams within ATL03 file
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL03 variables and attributes
    IS2_atl03_mds = {}
    IS2_atl03_attrs = {}

    #-- read each input beam within the file
    IS2_atl03_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        #-- check if subsetted beam contains data
        #-- check in both the geolocation and heights groups
        try:
            fileID[gtx]['geolocation']['segment_id']
            fileID[gtx]['heights']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl03_beams.append(gtx)

    #-- ICESat-2 spacecraft orientation at time
    IS2_atl03_mds['orbit_info'] = {}
    IS2_atl03_attrs['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl03_mds['orbit_info'][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Global Group Attributes
            for att_name,att_val in fileID['orbit_info'].attrs.items():
                IS2_atl03_attrs['orbit_info'][att_name] = att_val
            #-- Variable Attributes
            IS2_atl03_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['orbit_info'][key][att_name] = att_val

    #-- information ancillary to the data product
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl03_mds['ancillary_data'] = {}
    IS2_atl03_attrs['ancillary_data'] = {}
    ancillary_keys = ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc',
        'end_cycle','end_geoseg','end_gpssow','end_gpsweek','end_orbit',
        'end_region','end_rgt','granule_end_utc','granule_start_utc','release',
        'start_cycle','start_geoseg','start_gpssow','start_gpsweek',
        'start_orbit','start_region','start_rgt','version']
    for key in ancillary_keys:
        #-- get each HDF5 variable
        IS2_atl03_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl03_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl03_attrs['ancillary_data'][key][att_name] = att_val

    #-- transmit-echo-path (tep) parameters
    IS2_atl03_mds['ancillary_data']['tep'] = {}
    IS2_atl03_attrs['ancillary_data']['tep'] = {}
    for key,val in fileID['ancillary_data']['tep'].items():
        #-- get each HDF5 variable
        IS2_atl03_mds['ancillary_data']['tep'][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl03_attrs['ancillary_data']['tep'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['ancillary_data']['tep'][key][att_name] = att_val

    #-- channel dead time and first photon bias derived from ATLAS calibration
    cal1,cal2 = ('ancillary_data','calibrations')
    for var in ['dead_time','first_photon_bias']:
        IS2_atl03_mds[cal1][var] = {}
        IS2_atl03_attrs[cal1][var] = {}
        for key,val in fileID[cal1][cal2][var].items():
            #-- get each HDF5 variable
            if isinstance(val, h5py.Dataset):
                IS2_atl03_mds[cal1][var][key] = val[:]
            elif isinstance(val, h5py.Group):
                IS2_atl03_mds[cal1][var][key] = {}
                for k,v in val.items():
                    IS2_atl03_mds[cal1][var][key][k] = v[:]
            #-- Getting attributes of group and included variables
            if ATTRIBUTES:
                #-- Variable Attributes
                IS2_atl03_attrs[cal1][var][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[cal1][var][key][att_name] = att_val
                if isinstance(val, h5py.Group):
                    for k,v in val.items():
                        IS2_atl03_attrs[cal1][var][key][k] = {}
                        for att_name,att_val in val.attrs.items():
                            IS2_atl03_attrs[cal1][var][key][k][att_name]=att_val

    #-- get ATLAS impulse response variables for the transmitter echo path (TEP)
    tep1,tep2 = ('atlas_impulse_response','tep_histogram')
    IS2_atl03_mds[tep1] = {}
    IS2_atl03_attrs[tep1] = {}
    for pce in ['pce1_spot1','pce2_spot3']:
        IS2_atl03_mds[tep1][pce] = {tep2:{}}
        IS2_atl03_attrs[tep1][pce] = {tep2:{}}
        #-- for each TEP variable
        for key,val in fileID[tep1][pce][tep2].items():
            IS2_atl03_mds[tep1][pce][tep2][key] = val[:]
            #-- Getting attributes of included variables
            if ATTRIBUTES:
                #-- Global Group Attributes
                for att_name,att_val in fileID[tep1][pce][tep2].attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][att_name] = att_val
                #-- Variable Attributes
                IS2_atl03_attrs[tep1][pce][tep2][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][key][att_name] = att_val

    #-- Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl03_attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams)

#-- PURPOSE: read ICESat-2 ATL03 HDF5 data files for beam variables
def read_HDF5_ATL03_beam(FILENAME, gtx, ATTRIBUTES=False, VERBOSE=False):
    """
    Reads ICESat-2 ATL03 Global Geolocated Photons data files
    for a specific beam

    Arguments
    ---------
    FILENAME: full path to ATL03 file
    gtx: beam name based on ground track and position
        GT1L
        GT1R
        GT2L
        GT2R
        GT3L
        GT3R

    Keyword arguments
    -----------------
    ATTRIBUTES: read file, group and variable attributes
    VERBOSE: output information about top-level HDF5 groups

    Returns
    -------
    IS2_atl03_mds: dictionary with ATL03 beam-level variables
    IS2_atl03_attrs: dictionary with ATL03 beam-level attributes
    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL03 variables and attributes
    IS2_atl03_mds = {}
    IS2_atl03_attrs = {}

    #-- get each HDF5 variable
    IS2_atl03_mds['heights'] = {}
    IS2_atl03_mds['geolocation'] = {}
    IS2_atl03_mds['bckgrd_atlas'] = {}
    IS2_atl03_mds['geophys_corr'] = {}
    #-- ICESat-2 Measurement Group
    for key,val in fileID[gtx]['heights'].items():
        IS2_atl03_mds['heights'][key] = val[:]
    #-- ICESat-2 Geolocation Group
    for key,val in fileID[gtx]['geolocation'].items():
        IS2_atl03_mds['geolocation'][key] = val[:]
    #-- ICESat-2 Background Photon Rate Group
    for key,val in fileID[gtx]['bckgrd_atlas'].items():
        IS2_atl03_mds['bckgrd_atlas'][key] = val[:]
    #-- ICESat-2 Geophysical Corrections Group: Values for tides (ocean,
    #-- solid earth, pole, load, and equilibrium), inverted barometer (IB)
    #-- effects, and range corrections for tropospheric delays
    for key,val in fileID[gtx]['geophys_corr'].items():
        IS2_atl03_mds['geophys_corr'][key] = val[:]

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- Getting attributes of IS2_atl03_mds beam variables
        IS2_atl03_attrs['heights'] = {}
        IS2_atl03_attrs['geolocation'] = {}
        IS2_atl03_attrs['bckgrd_atlas'] = {}
        IS2_atl03_attrs['geophys_corr'] = {}
        IS2_atl03_attrs['Atlas_impulse_response'] = {}
        #-- Global Group Attributes
        for att_name,att_val in fileID[gtx].attrs.items():
            IS2_atl03_attrs[att_name] = att_val
        #-- ICESat-2 Measurement Group
        for key,val in fileID[gtx]['heights'].items():
            IS2_atl03_attrs['heights'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['heights'][key][att_name]=att_val
        #-- ICESat-2 Geolocation Group
        for key,val in fileID[gtx]['geolocation'].items():
            IS2_atl03_attrs['geolocation'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['geolocation'][key][att_name]=att_val
        #-- ICESat-2 Background Photon Rate Group
        for key,val in fileID[gtx]['bckgrd_atlas'].items():
            IS2_atl03_attrs['bckgrd_atlas'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['bckgrd_atlas'][key][att_name]=att_val
        #-- ICESat-2 Geophysical Corrections Group
        for key,val in fileID[gtx]['geophys_corr'].items():
            IS2_atl03_attrs['geophys_corr'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['geophys_corr'][key][att_name]=att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl03_mds,IS2_atl03_attrs)
