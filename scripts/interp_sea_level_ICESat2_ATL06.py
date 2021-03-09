#!/usr/bin/env python
u"""
interp_sea_level_ICESat2_ATL06.py
Written by Tyler Sutterley (02/2021)
Interpolates sea level anomalies (sla), absolute dynamic topography (adt) and
    mean dynamic topography (mdt) to times and locations of ICESat-2 ATL06 data
    This data will be extrapolated onto land points
    (masking will be needed for accurate assessments)

https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/
    global/msla-h.html
ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/
    dataset-duacs-rep-global-merged-allsat-phy-l4-v3

Note that the AVISO sea level data are gzip compressed netCDF4 files

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -V, --verbose: Output information about each created file
    -M X, --mode X: Permission mode of directories and files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    scikit-learn: Machine Learning in Python
        https://scikit-learn.org/stable/index.html
        https://github.com/scikit-learn/scikit-learn
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    read_ICESat2_ATL06.py: reads ICESat-2 land ice along-track height data files
    time.py: utilities for calculating time operations
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 02/2021: using argparse to set command line options
        using scikit-learn neighbors to interpolate fields
        time utilities for converting times from JD and to decimal
    Updated 01/2021: using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 10/2019: changing Y/N flags to True/False
    Updated 05/2019: check if beam exists in a try except else clause
    Updated 04/2019: check if subsetted beam contains land ice data
    Updated 02/2019: find valid points for both sea level and dynamic topography
        interpolate and output mean dynamic topography (mdt)
        using lists for python3 compatibility
    Updated 01/2018: additionally interpolates absolute dynamic topography (adt)
    Updated 02/2018: can update existing files and run only specific groups
    Updated 09/2017: reduce grid domains for faster processing times
    Written 09/2017
"""
from __future__ import print_function

import os
import re
import gzip
import h5py
import pyproj
import netCDF4
import argparse
import datetime
import numpy as np
import sklearn.neighbors
import icesat2_toolkit.time
from icesat2_toolkit.read_ICESat2_ATL06 import read_HDF5_ATL06

#-- PURPOSE: set the hemisphere of interest based on the granule
def set_hemisphere(GRANULE):
    if GRANULE in ('10','11','12'):
        projection_flag = 'S'
    elif GRANULE in ('03','04','05'):
        projection_flag = 'N'
    return projection_flag

#-- PURPOSE: interpolates to coordinates with inverse distance weighting
def inverse_distance(x, y, z, xi, yi, SEARCH='BallTree', N=10, POWER=2.0):
    #-- number of output points
    npts = len(xi)
    #-- create neighbors object for coordinates
    if (SEARCH == 'BallTree'):
        tree = sklearn.neighbors.BallTree(np.c_[x,y])
    elif (SEARCH == 'KDTree'):
        tree = sklearn.neighbors.KDTree(np.c_[x,y])
    #-- query the search tree to find the N closest points
    dist,indices = tree.query(np.c_[xi,yi], k=N, return_distance=True)
    #-- normalized weights if POWER > 0 (typically between 1 and 3)
    #-- in the inverse distance weighting
    power_inverse_distance = dist**(-POWER)
    s = np.sum(power_inverse_distance, axis=1)
    w = power_inverse_distance/np.broadcast_to(s[:,None],(npts,N))
    #-- calculate interpolated fields by inverse distance weighting
    return np.sum(w*z[indices],axis=1)

#-- PURPOSE interpolate sea level anomalies to lat/lon and then to time
def interpolate_sea_level(base_dir, xi, yi, CJD, HEM):
    #-- EPSG projections for converting lat/lon to polar stereographic
    EPSG = dict(N=3413,S=3031)
    #-- pyproj transformer for converting to polar stereographic
    crs1 = pyproj.CRS.from_string('epsg:4326')
    crs2 = pyproj.CRS.from_string(EPSG[HEM])
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    #-- interpolate mean dynamic topography
    input_file = 'mdt_cnes_cls2013_global.nc.gz'
    #-- read bytes from compressed file
    fd = gzip.open(os.path.join(base_dir,input_file),'rb')
    #-- dictionary with input fields
    dinput = {}
    #-- read netCDF file for mean dynamic topography
    with netCDF4.Dataset('mdt', mode='r', memory=fd.read()) as fileID:
        dinput['lon'] = fileID['lon'][:].copy()
        dinput['lat'] = fileID['lat'][:].copy()
        dinput['mdt'] = np.ma.array(fileID['mdt'][0,:,:].copy(),
            fill_value=fileID['mdt']._FillValue)
        dinput['mdt'].mask = (dinput['mdt'].data == dinput['mdt'].fill_value)
    #-- close the compressed file objects
    fd.close()
    #-- create 2-D grid coordinates from longitude and latitude vectors
    gridlon,gridlat = np.meshgrid(dinput['lon'],dinput['lat'])
    #-- convert from latitude/longitude into polar stereographic
    xg,yg = transformer.transform(gridlon,gridlat)

    #-- reduce to local coordinates to improve computational time
    gridmask = np.logical_not(dinput['mdt'].mask)
    if (HEM.upper() == 'N'):
        gridmask &= (gridlat >= 50.0)
    elif (HEM.upper() == 'S'):
        gridmask &= (gridlat <= -50.0)
    indy,indx = np.nonzero(gridmask)
    #-- calculate mean dynamic topography by inverse distance weighting
    MDT = inverse_distance(xg[indy,indx], yg[indy,indx],
        dinput['mdt'].data[indy,indx], xi, yi)

    #-- CNES Julian Days before and after measurement
    CJD1 = np.floor(CJD)
    #-- scale for linearly interpolating to date
    dt = (CJD - CJD1[0])
    #-- output sea level anomaly and absolute dynamic topography
    SLA = np.zeros_like(CJD)
    ADT = np.zeros_like(CJD)
    #-- for the range of dates
    for day in range(2):
        #-- convert from CNES Julians Days to calendar dates for time
        JD1 = CJD1 + day + 2433282.5
        YY,MM,DD,HH,MN,SS = icesat2_toolkit.time.convert_julian(JD1[0],
            FORMAT='tuple', ASTYPE=np.int)
        #-- sea level directory
        ddir = os.path.join(base_dir, '{0:0.0f}'.format(YY))
        #-- input file for day before the measurement
        regex = re.compile(('dt_global_allsat_phy_l4_{0:4d}{1:02d}{2:02d}_'
            '(\d{{4}})(\d{{2}})(\d{{2}}).nc.gz').format(YY,MM,DD))
        input_file, = [fi for fi in os.listdir(ddir) if regex.match(fi)]
        #-- dictionary with input fields
        dinput = {}
        #-- read bytes from compressed file
        fd = gzip.open(os.path.join(ddir,input_file),'rb')
        #-- read netCDF file for time
        with netCDF4.Dataset('sla', mode='r', memory=fd.read()) as fileID:
            dinput['lon'] = fileID['lon'][:].copy()
            dinput['lat'] = fileID['lat'][:].copy()
            dinput['sla'] = np.ma.array(fileID['sla'][0,:,:].copy(),
                fill_value=fileID['sla']._FillValue)
            dinput['adt'] = np.ma.array(fileID['adt'][0,:,:].copy(),
                fill_value=fileID['adt']._FillValue)
        #-- close the compressed file objects
        fd.close()
        #-- for each variable to interpolate
        out = {}
        for var in ['sla','adt']:
            #-- reduce to local coordinates to improve computational time
            gridmask = np.logical_not(dinput[var].mask)
            if (HEM.upper() == 'N'):
                gridmask &= (gridlat >= 50.0)
            elif (HEM.upper() == 'S'):
                gridmask &= (gridlat <= -50.0)
            indy,indx = np.nonzero(gridmask)
            #-- calculate variable by inverse distance weighting
            out[var] = inverse_distance(xg[indy,indx], yg[indy,indx],
                dinput[var].data[indy,indx], xi, yi)
        #-- linearly interpolate to date for iteration
        SLA += out['sla']*(2.0*dt*day - dt - day + 1.0)
        ADT += out['adt']*(2.0*dt*day - dt - day + 1.0)
    #-- return interpolated values
    return (MDT,SLA,ADT)

#-- PURPOSE: read ICESat-2 land ice data (ATL06) from NSIDC
#-- interpolate AVISO sea level at points and times
def interp_sea_level_ICESat2(base_dir, FILE, VERBOSE=False, MODE=0o775):

    #-- read data from input_file
    print('{0} -->'.format(os.path.basename(FILE))) if VERBOSE else None
    IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_beams = read_HDF5_ATL06(FILE,
        ATTRIBUTES=True)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = rx.findall(FILE).pop()
    #-- set the hemisphere flag based on ICESat-2 granule
    HEM = set_hemisphere(GRAN)

    #-- HDF5 file attributes
    attrib = {}
    #-- mean dynamic topography
    attrib['h_mdt'] = {}
    attrib['h_mdt']['long_name'] = 'Mean Dynamic Topography'
    attrib['h_mdt']['description'] = 'Sea surface height above geoid'
    attrib['h_mdt']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')
    #-- sea level anomalies
    attrib['h_sla'] = {}
    attrib['h_sla']['long_name'] = 'Sea Level Anomaly'
    attrib['h_sla']['description'] = 'Sea surface anomalies'
    attrib['h_sla']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')
    #-- absolute dynamic topography
    attrib['h_adt'] = {}
    attrib['h_adt']['long_name'] = 'Absolute Dynamic Topography'
    attrib['h_adt']['description'] = ('Sea surface height above geoid calculated '
        'by adding the mean dynamic topography to the sea level anomalies')
    attrib['h_adt']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')

    #-- EPSG projections for converting lat/lon to polar stereographic
    EPSG = dict(N=3413,S=3031)
    #-- pyproj transformer for converting to polar stereographic
    crs1 = pyproj.CRS.from_string('epsg:4326')
    crs2 = pyproj.CRS.from_string(EPSG[HEM])
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl06_mds['ancillary_data']['atlas_sdp_gps_epoch']

    #-- copy variables for outputting to HDF5 file
    IS2_atl06_corr = {}
    IS2_atl06_fill = {}
    IS2_atl06_dims = {}
    IS2_atl06_corr_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl06_corr['ancillary_data'] = {}
    IS2_atl06_corr_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl06_corr['ancillary_data'][key] = IS2_atl06_mds['ancillary_data'][key]
        #-- Getting attributes of group and included variables
        IS2_atl06_corr_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][key].items():
            IS2_atl06_corr_attrs['ancillary_data'][key][att_name] = att_val
    #-- for each input beam within the file
    for gtx in sorted(IS2_atl06_beams):
        #-- output data dictionaries for beam
        IS2_atl06_corr[gtx] = dict(land_ice_segments={})
        IS2_atl06_fill[gtx] = dict(land_ice_segments={})
        IS2_atl06_dims[gtx] = dict(land_ice_segments={})
        IS2_atl06_corr_attrs[gtx] = dict(land_ice_segments={})

        #-- number of segments
        val = IS2_atl06_mds[gtx]['land_ice_segments']
        n_seg = len(val['segment_id'])
        #-- find valid segments for beam
        fv = IS2_atl06_attrs[gtx]['land_ice_segments']['h_li']['_FillValue']
        mask = (val['h_li'] == fv)
        valid, = np.nonzero(np.logical_not(mask))

        #-- convert time from ATLAS SDP to CNES JD
        #-- days relative to 1950-01-01T00:00:00
        gps_seconds = atlas_sdp_gps_epoch + val['delta_time']
        leap_seconds = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
        cnes_time = icesat2_toolkit.time.convert_delta_time(gps_seconds-leap_seconds,
            epoch1=(1980,1,6,0,0,0), epoch2=(1950,1,1,0,0,0), scale=1.0/86400.0)

        #-- extract lat/lon and convert to polar stereographic
        X,Y = transformer.transform(val['longitude'],val['latitude'])

        #-- interpolate sea level anomalies and dynamic topographies
        MDT,SLA,ADT = interpolate_sea_level(base_dir,X,Y,cnes_time,HEM)
        #-- map to original segments
        sea_level = {}
        sea_level['h_mdt'] = np.ma.zeros((n_seg),fill_value=fv)
        sea_level['h_mdt'].data[valid] = np.copy(MDT)
        sea_level['h_sla'] = np.ma.zeros((n_seg),fill_value=fv)
        sea_level['h_sla'].data[valid] = np.copy(SLA)
        sea_level['h_adt'] = np.ma.zeros((n_seg),fill_value=fv)
        sea_level['h_adt'].data[valid] = np.copy(ADT)
        #-- replace nan values with fill value
        for key in ['h_mdt','h_sla','h_adt']:
            sea_level[key].mask = np.copy(mask)
            sea_level[key].mask |= np.isnan(sea_level[key].data)
            sea_level[key].data[sea_level[key].mask] = fv

        #-- group attributes for beam
        IS2_atl06_corr_attrs[gtx]['Description'] = IS2_atl06_attrs[gtx]['Description']
        IS2_atl06_corr_attrs[gtx]['atlas_pce'] = IS2_atl06_attrs[gtx]['atlas_pce']
        IS2_atl06_corr_attrs[gtx]['atlas_beam_type'] = IS2_atl06_attrs[gtx]['atlas_beam_type']
        IS2_atl06_corr_attrs[gtx]['groundtrack_id'] = IS2_atl06_attrs[gtx]['groundtrack_id']
        IS2_atl06_corr_attrs[gtx]['atmosphere_profile'] = IS2_atl06_attrs[gtx]['atmosphere_profile']
        IS2_atl06_corr_attrs[gtx]['atlas_spot_number'] = IS2_atl06_attrs[gtx]['atlas_spot_number']
        IS2_atl06_corr_attrs[gtx]['sc_orientation'] = IS2_atl06_attrs[gtx]['sc_orientation']
        #-- group attributes for land_ice_segments
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        #-- geolocation, time and segment ID
        #-- delta time
        delta_time = np.ma.array(val['delta_time'], fill_value=fv,
            mask=(val['delta_time']==fv))
        IS2_atl06_corr[gtx]['land_ice_segments']['delta_time'] = delta_time
        IS2_atl06_fill[gtx]['land_ice_segments']['delta_time'] = delta_time.fill_value
        IS2_atl06_dims[gtx]['land_ice_segments']['delta_time'] = None
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['standard_name'] = "time"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['calendar'] = "standard"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['delta_time']['coordinates'] = \
            "segment_id latitude longitude"
        #-- latitude
        latitude = np.ma.array(val['latitude'], fill_value=fv,
            mask=(val['latitude']==fv))
        IS2_atl06_corr[gtx]['land_ice_segments']['latitude'] = latitude
        IS2_atl06_fill[gtx]['land_ice_segments']['latitude'] = latitude.fill_value
        IS2_atl06_dims[gtx]['land_ice_segments']['latitude'] = ['delta_time']
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['units'] = "degrees_north"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['valid_min'] = -90.0
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['valid_max'] = 90.0
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['latitude']['coordinates'] = \
            "segment_id delta_time longitude"
        #-- longitude
        longitude = np.ma.array(val['longitude'], fill_value=fv,
            mask=(val['longitude']==fv))
        IS2_atl06_corr[gtx]['land_ice_segments']['longitude'] = longitude
        IS2_atl06_fill[gtx]['land_ice_segments']['longitude'] = longitude.fill_value
        IS2_atl06_dims[gtx]['land_ice_segments']['longitude'] = ['delta_time']
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['units'] = "degrees_east"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['valid_min'] = -180.0
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['valid_max'] = 180.0
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['longitude']['coordinates'] = \
            "segment_id delta_time latitude"
        #-- segment ID
        IS2_atl06_corr[gtx]['land_ice_segments']['segment_id'] = val['segment_id']
        IS2_atl06_fill[gtx]['land_ice_segments']['segment_id'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['segment_id'] = ['delta_time']
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id']['units'] = "1"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id']['contentType'] = "referenceInformation"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id']['long_name'] = "Along-track segment ID number"
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id']['description'] = ("A 7 digit number "
            "identifying the along-track geolocation segment number.  These are sequential, starting with "
            "1 for the first segment after an ascending equatorial crossing node. Equal to the segment_id for "
            "the second of the two 20m ATL03 segments included in the 40m ATL06 segment")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['segment_id']['coordinates'] = \
            "delta_time latitude longitude"

        #-- geophysical variables
        IS2_atl06_corr[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_dims[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical']['Description'] = ("The geophysical group "
            "contains parameters used to correct segment heights for geophysical effects, parameters "
            "related to solar background and parameters indicative of the presence or absence of clouds.")
        IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        #-- interpolated sea level products
        for key,val in sea_level.items():
            #-- add to output
            IS2_atl06_corr[gtx]['land_ice_segments']['geophysical'][key] = val
            IS2_atl06_fill[gtx]['land_ice_segments']['geophysical'][key] = fv
            IS2_atl06_dims[gtx]['land_ice_segments']['geophysical'][key] = ['delta_time']
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key] = {}
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['units'] = "meters"
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['contentType'] = "referenceInformation"
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['long_name'] = attrib[key]['long_name']
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['description'] = attrib[key]['description']
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['source'] = 'AVISO/Copernicus'
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['reference'] = attrib[key]['reference']
            IS2_atl06_corr_attrs[gtx]['land_ice_segments']['geophysical'][key]['coordinates'] = \
                "../segment_id ../delta_time ../latitude ../longitude"

    #-- output HDF5 files with interpolated sea level data
    args = (PRD,'AVISO_SEA_LEVEL',YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX)
    file_format='{0}_{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
    #-- print file information
    print('\t{0}'.format(file_format.format(*args))) if VERBOSE else None
    HDF5_ATL06_corr_write(IS2_atl06_corr, IS2_atl06_corr_attrs,
        CLOBBER=True, INPUT=os.path.basename(FILE),
        FILL_VALUE=IS2_atl06_fill, DIMENSIONS=IS2_atl06_dims,
        FILENAME=os.path.join(DIRECTORY,file_format.format(*args)))
    #-- change the permissions mode
    os.chmod(os.path.join(DIRECTORY,file_format.format(*args)), MODE)

#-- PURPOSE: outputting the correction values for ICESat-2 data to HDF5
def HDF5_ATL06_corr_write(IS2_atl06_corr, IS2_atl06_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CLOBBER=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)

    #-- create HDF5 records
    h5 = {}

    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k,v in IS2_atl06_corr['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl06_corr.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        #-- add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl06_attrs[gtx][att_name]
        #-- create land_ice_segments group
        fileID[gtx].create_group('land_ice_segments')
        h5[gtx] = dict(land_ice_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl06_attrs[gtx]['land_ice_segments'][att_name]
            fileID[gtx]['land_ice_segments'].attrs[att_name] = att_val

        #-- delta_time, geolocation and segment_id variables
        for k in ['delta_time','latitude','longitude','segment_id']:
            #-- values and attributes
            v = IS2_atl06_corr[gtx]['land_ice_segments'][k]
            attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments',k)
            if fillvalue:
                h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                    compression='gzip')
            else:
                h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, compression='gzip')
            #-- create or attach dimensions for HDF5 variable
            if DIMENSIONS[gtx]['land_ice_segments'][k]:
                #-- attach dimensions
                for i,dim in enumerate(DIMENSIONS[gtx]['land_ice_segments'][k]):
                    h5[gtx]['land_ice_segments'][k].dims[i].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
            else:
                #-- make dimension
                h5[gtx]['land_ice_segments'][k].make_scale(k)
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][k].attrs[att_name] = att_val

        #-- add to geophysical corrections
        key = 'geophysical'
        fileID[gtx]['land_ice_segments'].create_group(key)
        h5[gtx]['land_ice_segments'][key] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name]
            fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
        for k,v in IS2_atl06_corr[gtx]['land_ice_segments'][key].items():
            #-- attributes
            attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][key][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][key][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}/{3}'.format(gtx,'land_ice_segments',key,k)
            if fillvalue:
                h5[gtx]['land_ice_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
            else:
                h5[gtx]['land_ice_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, compression='gzip')
            #-- attach dimensions
            for i,dim in enumerate(DIMENSIONS[gtx]['land_ice_segments'][key][k]):
                h5[gtx]['land_ice_segments'][key][k].dims[i].attach_scale(
                    h5[gtx]['land_ice_segments'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Land Ice Height'
    fileID.attrs['summary'] = ('Estimates of the ice-sheet correction parameters '
        'needed to interpret and assess the quality of land height estimates.')
    fileID.attrs['description'] = ('Land ice parameters for each beam.  All '
        'parameters are calculated for the same along-track increments for '
        'each beam and repeat.')
    date_created = datetime.datetime.today()
    fileID.attrs['date_created'] = date_created.isoformat()
    project = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = project
    platform = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = platform
    #-- add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    #-- add attributes for input ATL06 file
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl06_corr[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl06_corr[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl06_corr[gtx]['land_ice_segments']['delta_time']
        #-- setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time.min() if (delta_time.min() < tmn) else tmn
        tmx = delta_time.max() if (delta_time.max() > tmx) else tmx
    #-- add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    #-- convert start and end time from ATLAS SDP seconds into GPS seconds
    atlas_sdp_gps_epoch=IS2_atl06_corr['ancillary_data']['atlas_sdp_gps_epoch']
    gps_seconds = atlas_sdp_gps_epoch + np.array([tmn,tmx])
    #-- calculate leap seconds
    leaps = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
    #-- convert from seconds since 1980-01-06T00:00:00 to Modified Julian days
    MJD = icesat2_toolkit.time.convert_delta_time(gps_seconds - leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date
    YY,MM,DD,HH,MN,SS = icesat2_toolkit.time.convert_julian(MJD + 2400000.5,
        FORMAT='tuple')
    #-- add attributes with measurement date start, end and duration
    tcs = datetime.datetime(np.int(YY[0]), np.int(MM[0]), np.int(DD[0]),
        np.int(HH[0]), np.int(MN[0]), np.int(SS[0]), np.int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(np.int(YY[1]), np.int(MM[1]), np.int(DD[1]),
        np.int(HH[1]), np.int(MN[1]), np.int(SS[1]), np.int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    #-- Closing the HDF5 file
    fileID.close()

#-- Main program that calls interp_sea_level_ICESat2()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Interpolates AVISO sea level anomalies, absolute
            dynamic topography and mean dynamic topography to ICESat-2
            ATL06 land ice elevation data
            """
    )
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat-2 ATL06 file to run')
    #-- directory with sea level data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information about each created file')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files created')
    args = parser.parse_args()

    #-- run for each input ATL06 file
    for FILE in args.infile:
        interp_sea_level_ICESat2(args.directory, FILE,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()