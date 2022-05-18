#!/usr/bin/env python
u"""
interp_sea_level_ICESat2_ATL11.py
Written by Tyler Sutterley (05/2022)
Interpolates sea level anomalies (sla), absolute dynamic topography (adt) and
    mean dynamic topography (mdt) to times and locations of ICESat-2 ATL11 data
    This data will be extrapolated onto land points
    (masking will be needed for accurate assessments)

https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/
    global/msla-h.html
ftp://ftp.sltac.cls.fr/Core/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/
    dataset-duacs-rep-global-merged-allsat-phy-l4-v3

Note that the AVISO sea level data are gzip compressed netCDF4 files

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -C, --crossovers: Run ATL11 Crossovers
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
    read_ICESat2_ATL11.py: reads ICESat-2 annual land ice height data files
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Written 02/2021
"""
from __future__ import print_function

import os
import re
import gzip
import h5py
import pyproj
import logging
import netCDF4
import argparse
import datetime
import numpy as np
import collections
import sklearn.neighbors
import icesat2_toolkit.time
import icesat2_toolkit.utilities
from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11

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
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG[HEM]))
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
            FORMAT='tuple', ASTYPE=int)
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

#-- PURPOSE: read ICESat-2 annual land ice height data (ATL11) from NSIDC
#-- interpolate AVISO sea level at points and times
def interp_sea_level_ICESat2(base_dir, FILE, CROSSOVERS=False, VERBOSE=False,
    MODE=0o775):

    #-- create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- read data from input file
    logging.info('{0} -->'.format(os.path.basename(FILE)))
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = read_HDF5_ATL11(FILE,
        ATTRIBUTES=True, CROSSOVERS=CROSSOVERS)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_'
        r'(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()
    #-- set the hemisphere flag based on ICESat-2 granule
    HEM = set_hemisphere(GRAN)

    #-- HDF5 file attributes
    attrib = {}
    #-- mean dynamic topography
    attrib['mdt'] = {}
    attrib['mdt']['long_name'] = 'Mean Dynamic Topography'
    attrib['mdt']['description'] = 'Sea surface height above geoid'
    attrib['mdt']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')
    #-- sea level anomalies
    attrib['sla'] = {}
    attrib['sla']['long_name'] = 'Sea Level Anomaly'
    attrib['sla']['description'] = 'Sea surface anomalies'
    attrib['sla']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')
    #-- absolute dynamic topography
    attrib['adt'] = {}
    attrib['adt']['long_name'] = 'Absolute Dynamic Topography'
    attrib['adt']['description'] = ('Sea surface height above geoid calculated '
        'by adding the mean dynamic topography to the sea level anomalies')
    attrib['adt']['reference'] = ('https://www.aviso.altimetry.fr/en/data/'
        'products/sea-surface-height-products/global/msla-h.html')

    #-- EPSG projections for converting lat/lon to polar stereographic
    EPSG = dict(N=3413,S=3031)
    #-- pyproj transformer for converting to polar stereographic
    crs1 = pyproj.CRS.from_string('epsg:4326')
    crs2 = pyproj.CRS.from_string(EPSG[HEM])
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl11_mds['ancillary_data']['atlas_sdp_gps_epoch']

    #-- copy variables for outputting to HDF5 file
    IS2_atl11_corr = {}
    IS2_atl11_fill = {}
    IS2_atl11_dims = {}
    IS2_atl11_corr_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl11_corr['ancillary_data'] = {}
    IS2_atl11_corr_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl11_corr['ancillary_data'][key] = IS2_atl11_mds['ancillary_data'][key]
        #-- Getting attributes of group and included variables
        IS2_atl11_corr_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][key].items():
            IS2_atl11_corr_attrs['ancillary_data'][key][att_name] = att_val
    #-- HDF5 group name for across-track data
    XT = 'crossing_track_data'

    #-- for each input beam pair within the file
    for ptx in sorted(IS2_atl11_pairs):
        #-- output data dictionaries for beam pair
        IS2_atl11_corr[ptx] = dict(cycle_stats=collections.OrderedDict(),
            crossing_track_data=collections.OrderedDict())
        IS2_atl11_fill[ptx] = dict(cycle_stats={},crossing_track_data={})
        IS2_atl11_dims[ptx] = dict(cycle_stats={},crossing_track_data={})
        IS2_atl11_corr_attrs[ptx] = dict(cycle_stats={},crossing_track_data={})

        #-- extract along-track and across-track variables
        ref_pt = {}
        latitude = {}
        longitude = {}
        delta_time = {}
        groups = ['AT']
        #-- dictionary with output sea level variables
        MDT,SLA,ADT = ({},{},{})
        #-- number of average segments and number of included cycles
        #-- fill_value for invalid heights and corrections
        fv = IS2_atl11_attrs[ptx]['h_corr']['_FillValue']
        #-- shape of along-track data
        n_points,n_cycles = IS2_atl11_mds[ptx]['delta_time'].shape
        #-- along-track (AT) reference point, latitude, longitude and time
        ref_pt['AT'] = IS2_atl11_mds[ptx]['ref_pt'].copy()
        latitude['AT'] = np.ma.array(IS2_atl11_mds[ptx]['latitude'],
            fill_value=IS2_atl11_attrs[ptx]['latitude']['_FillValue'])
        latitude['AT'].mask = (latitude['AT'] == latitude['AT'].fill_value)
        longitude['AT'] = np.ma.array(IS2_atl11_mds[ptx]['longitude'],
            fill_value=IS2_atl11_attrs[ptx]['longitude']['_FillValue'])
        longitude['AT'].mask = (longitude['AT'] == longitude['AT'].fill_value)
        delta_time['AT'] = np.ma.array(IS2_atl11_mds[ptx]['delta_time'],
            fill_value=IS2_atl11_attrs[ptx]['delta_time']['_FillValue'])
        delta_time['AT'].mask = (delta_time['AT'] == delta_time['AT'].fill_value)
        #-- along-track (AT) sea level corrections
        MDT['AT'] = np.ma.empty((n_points,n_cycles),fill_value=fv)
        MDT['AT'].mask = (delta_time['AT'] == delta_time['AT'].fill_value)
        SLA['AT'] = np.ma.empty((n_points,n_cycles),fill_value=fv)
        SLA['AT'].mask = (delta_time['AT'] == delta_time['AT'].fill_value)
        ADT['AT'] = np.ma.empty((n_points,n_cycles),fill_value=fv)
        ADT['AT'].mask = (delta_time['AT'] == delta_time['AT'].fill_value)
        #-- if running ATL11 crossovers
        if CROSSOVERS:
            #-- add to group
            groups.append('XT')
            #-- shape of across-track data
            n_cross, = IS2_atl11_mds[ptx][XT]['delta_time'].shape
            #-- across-track (XT) reference point, latitude, longitude and time
            ref_pt['XT'] = IS2_atl11_mds[ptx][XT]['ref_pt'].copy()
            latitude['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['latitude'],
                fill_value=IS2_atl11_attrs[ptx][XT]['latitude']['_FillValue'])
            latitude['XT'].mask = (latitude['XT'] == latitude['XT'].fill_value)
            longitude['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['longitude'],
                fill_value=IS2_atl11_attrs[ptx][XT]['longitude']['_FillValue'])
            latitude['XT'].mask = (latitude['XT'] == longitude['XT'].fill_value)
            delta_time['XT'] = np.ma.array(IS2_atl11_mds[ptx][XT]['delta_time'],
                fill_value=IS2_atl11_attrs[ptx][XT]['delta_time']['_FillValue'])
            delta_time['XT'].mask = (delta_time['XT'] == delta_time['XT'].fill_value)
            #-- across-track (XT) sea level corrections
            MDT['XT'] = np.ma.empty((n_cross),fill_value=fv)
            MDT['XT'].mask = (delta_time['XT'] == delta_time['XT'].fill_value)
            SLA['XT'] = np.ma.empty((n_cross),fill_value=fv)
            SLA['XT'].mask = (delta_time['XT'] == delta_time['XT'].fill_value)
            ADT['XT'] = np.ma.empty((n_cross),fill_value=fv)
            ADT['XT'].mask = (delta_time['XT'] == delta_time['XT'].fill_value)

        #-- calculate corrections for along-track and across-track data
        for track in groups:
            #-- convert time from ATLAS SDP to CNES Julian Days
            #-- days relative to 1950-01-01T00:00:00
            gps_seconds = atlas_sdp_gps_epoch + delta_time[track]
            leap_seconds = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
            cnes_time = icesat2_toolkit.time.convert_delta_time(gps_seconds-leap_seconds,
                epoch1=(1980,1,6,0,0,0), epoch2=(1950,1,1,0,0,0), scale=1.0/86400.0)

            #-- extract lat/lon and convert to polar stereographic
            X,Y = transformer.transform(longitude[track],longitude[track])

            #-- calculate sea level corrections for track type
            if (track == 'AT'):
                #-- calculate for each cycle if along-track
                for cycle in range(n_cycles):
                    #-- interpolate sea level anomalies and dynamic topographies
                    MDT[track][:,cycle],SLA[track][:,cycle],ADT[track][:,cycle] = \
                        interpolate_sea_level(base_dir,X,Y,cnes_time[:,cycle],HEM)
            elif (track == 'XT'):
                #-- for each unique CNES day to interpolate in the crossovers
                CJD,inverse = np.unique(np.floor(cnes_time),return_inverse=True)
                for indice,_ in enumerate(CJD):
                    #-- indices in original arrays for the CNES day
                    i, = np.nonzero(inverse == indice)
                    #-- interpolate sea level anomalies and dynamic topographies
                    MDT[track][i],SLA[track][i],ADT[track][i] = \
                        interpolate_sea_level(base_dir,X[i],Y[i],cnes_time[i],HEM)

        #-- group attributes for beam
        IS2_atl11_corr_attrs[ptx]['description'] = ('Contains the primary science parameters '
            'for this data set')
        IS2_atl11_corr_attrs[ptx]['beam_pair'] = IS2_atl11_attrs[ptx]['beam_pair']
        IS2_atl11_corr_attrs[ptx]['ReferenceGroundTrack'] = IS2_atl11_attrs[ptx]['ReferenceGroundTrack']
        IS2_atl11_corr_attrs[ptx]['first_cycle'] = IS2_atl11_attrs[ptx]['first_cycle']
        IS2_atl11_corr_attrs[ptx]['last_cycle'] = IS2_atl11_attrs[ptx]['last_cycle']
        IS2_atl11_corr_attrs[ptx]['equatorial_radius'] = IS2_atl11_attrs[ptx]['equatorial_radius']
        IS2_atl11_corr_attrs[ptx]['polar_radius'] = IS2_atl11_attrs[ptx]['polar_radius']

        #-- geolocation, time and reference point
        #-- reference point
        IS2_atl11_corr[ptx]['ref_pt'] = ref_pt['AT'].copy()
        IS2_atl11_fill[ptx]['ref_pt'] = None
        IS2_atl11_dims[ptx]['ref_pt'] = None
        IS2_atl11_corr_attrs[ptx]['ref_pt'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['ref_pt']['units'] = "1"
        IS2_atl11_corr_attrs[ptx]['ref_pt']['contentType'] = "referenceInformation"
        IS2_atl11_corr_attrs[ptx]['ref_pt']['long_name'] = "Reference point number"
        IS2_atl11_corr_attrs[ptx]['ref_pt']['source'] = "ATL06"
        IS2_atl11_corr_attrs[ptx]['ref_pt']['description'] = ("The reference point is the "
            "7 digit segment_id number corresponding to the center of the ATL06 data used "
            "for each ATL11 point.  These are sequential, starting with 1 for the first "
            "segment after an ascending equatorial crossing node.")
        IS2_atl11_corr_attrs[ptx]['ref_pt']['coordinates'] = \
            "delta_time latitude longitude"
        #-- cycle_number
        IS2_atl11_corr[ptx]['cycle_number'] = IS2_atl11_mds[ptx]['cycle_number'].copy()
        IS2_atl11_fill[ptx]['cycle_number'] = None
        IS2_atl11_dims[ptx]['cycle_number'] = None
        IS2_atl11_corr_attrs[ptx]['cycle_number'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['cycle_number']['units'] = "1"
        IS2_atl11_corr_attrs[ptx]['cycle_number']['long_name'] = "Orbital cycle number"
        IS2_atl11_corr_attrs[ptx]['cycle_number']['source'] = "ATL06"
        IS2_atl11_corr_attrs[ptx]['cycle_number']['description'] = ("Number of 91-day periods "
            "that have elapsed since ICESat-2 entered the science orbit. Each of the 1,387 "
            "reference ground track (RGTs) is targeted in the polar regions once "
            "every 91 days.")
        #-- delta time
        IS2_atl11_corr[ptx]['delta_time'] = delta_time['AT'].copy()
        IS2_atl11_fill[ptx]['delta_time'] = delta_time['AT'].fill_value
        IS2_atl11_dims[ptx]['delta_time'] = ['ref_pt','cycle_number']
        IS2_atl11_corr_attrs[ptx]['delta_time'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl11_corr_attrs[ptx]['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl11_corr_attrs[ptx]['delta_time']['standard_name'] = "time"
        IS2_atl11_corr_attrs[ptx]['delta_time']['calendar'] = "standard"
        IS2_atl11_corr_attrs[ptx]['delta_time']['source'] = "ATL06"
        IS2_atl11_corr_attrs[ptx]['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl11_corr_attrs[ptx]['delta_time']['coordinates'] = \
            "ref_pt cycle_number latitude longitude"
        #-- latitude
        IS2_atl11_corr[ptx]['latitude'] = latitude['AT'].copy()
        IS2_atl11_fill[ptx]['latitude'] = latitude['AT'].fill_value
        IS2_atl11_dims[ptx]['latitude'] = ['ref_pt']
        IS2_atl11_corr_attrs[ptx]['latitude'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['latitude']['units'] = "degrees_north"
        IS2_atl11_corr_attrs[ptx]['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_corr_attrs[ptx]['latitude']['long_name'] = "Latitude"
        IS2_atl11_corr_attrs[ptx]['latitude']['standard_name'] = "latitude"
        IS2_atl11_corr_attrs[ptx]['latitude']['source'] = "ATL06"
        IS2_atl11_corr_attrs[ptx]['latitude']['description'] = ("Center latitude of "
            "selected segments")
        IS2_atl11_corr_attrs[ptx]['latitude']['valid_min'] = -90.0
        IS2_atl11_corr_attrs[ptx]['latitude']['valid_max'] = 90.0
        IS2_atl11_corr_attrs[ptx]['latitude']['coordinates'] = \
            "ref_pt delta_time longitude"
        #-- longitude
        IS2_atl11_corr[ptx]['longitude'] = longitude['AT'].copy()
        IS2_atl11_fill[ptx]['longitude'] = longitude['AT'].fill_value
        IS2_atl11_dims[ptx]['longitude'] = ['ref_pt']
        IS2_atl11_corr_attrs[ptx]['longitude'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['longitude']['units'] = "degrees_east"
        IS2_atl11_corr_attrs[ptx]['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_corr_attrs[ptx]['longitude']['long_name'] = "Longitude"
        IS2_atl11_corr_attrs[ptx]['longitude']['standard_name'] = "longitude"
        IS2_atl11_corr_attrs[ptx]['longitude']['source'] = "ATL06"
        IS2_atl11_corr_attrs[ptx]['longitude']['description'] = ("Center longitude of "
            "selected segments")
        IS2_atl11_corr_attrs[ptx]['longitude']['valid_min'] = -180.0
        IS2_atl11_corr_attrs[ptx]['longitude']['valid_max'] = 180.0
        IS2_atl11_corr_attrs[ptx]['longitude']['coordinates'] = \
            "ref_pt delta_time latitude"

        #-- cycle statistics variables
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['Description'] = ("The cycle_stats subgroup "
            "contains summary information about segments for each reference point, including "
            "the uncorrected mean heights for reference surfaces, blowing snow and cloud "
            "indicators, and geolocation and height misfit statistics.")
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['data_rate'] = ("Data within this group "
            "are stored at the average segment rate.")

        #-- interpolated sea level products
        sea_level = dict(mdt=MDT['AT'],sla=SLA['AT'],adt=ADT['AT'])
        for key,val in sea_level.items():
            #-- add to output
            IS2_atl11_corr[ptx]['cycle_stats'][key] = val.copy()
            IS2_atl11_fill[ptx]['cycle_stats'][key] = val.fill_value
            IS2_atl11_dims[ptx]['cycle_stats'][key] = ['ref_pt','cycle_number']
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['units'] = "meters"
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['contentType'] = "referenceInformation"
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['long_name'] = attrib[key]['long_name']
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['description'] = attrib[key]['description']
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['source'] = 'AVISO/Copernicus'
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['reference'] = attrib[key]['reference']
            IS2_atl11_corr_attrs[ptx]['cycle_stats'][key]['coordinates'] = \
                "../ref_pt ../cycle_number ../delta_time ../latitude ../longitude"

        #-- if crossover measurements were calculated
        if CROSSOVERS:
            #-- crossing track variables
            IS2_atl11_corr_attrs[ptx][XT]['Description'] = ("The crossing_track_data "
                "subgroup contains elevation data at crossover locations. These are "
                "locations where two ICESat-2 pair tracks cross, so data are available "
                "from both the datum track, for which the granule was generated, and "
                "from the crossing track.")
            IS2_atl11_corr_attrs[ptx][XT]['data_rate'] = ("Data within this group are "
                "stored at the average segment rate.")

            #-- reference point
            IS2_atl11_corr[ptx][XT]['ref_pt'] = ref_pt['XT'].copy()
            IS2_atl11_fill[ptx][XT]['ref_pt'] = None
            IS2_atl11_dims[ptx][XT]['ref_pt'] = None
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['units'] = "1"
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['contentType'] = "referenceInformation"
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['long_name'] = ("fit center reference point number, "
                "segment_id")
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['source'] = "derived, ATL11 algorithm"
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['description'] = ("The reference-point number of the "
                "fit center for the datum track. The reference point is the 7 digit segment_id number "
                "corresponding to the center of the ATL06 data used for each ATL11 point.  These are "
                "sequential, starting with 1 for the first segment after an ascending equatorial "
                "crossing node.")
            IS2_atl11_corr_attrs[ptx][XT]['ref_pt']['coordinates'] = \
                "delta_time latitude longitude"

            #-- reference ground track of the crossing track
            IS2_atl11_corr[ptx][XT]['rgt'] = IS2_atl11_mds[ptx][XT]['rgt'].copy()
            IS2_atl11_fill[ptx][XT]['rgt'] = IS2_atl11_attrs[ptx][XT]['rgt']['_FillValue']
            IS2_atl11_dims[ptx][XT]['rgt'] = None
            IS2_atl11_corr_attrs[ptx][XT]['rgt'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['units'] = "1"
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['contentType'] = "referenceInformation"
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['long_name'] = "crossover reference ground track"
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['source'] = "ATL06"
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['description'] = "The RGT number for the crossing data."
            IS2_atl11_corr_attrs[ptx][XT]['rgt']['coordinates'] = \
                "ref_pt delta_time latitude longitude"
            #-- cycle_number of the crossing track
            IS2_atl11_corr[ptx][XT]['cycle_number'] = IS2_atl11_mds[ptx][XT]['cycle_number'].copy()
            IS2_atl11_fill[ptx][XT]['cycle_number'] = IS2_atl11_attrs[ptx][XT]['cycle_number']['_FillValue']
            IS2_atl11_dims[ptx][XT]['cycle_number'] = None
            IS2_atl11_corr_attrs[ptx][XT]['cycle_number'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['cycle_number']['units'] = "1"
            IS2_atl11_corr_attrs[ptx][XT]['cycle_number']['long_name'] = "crossover cycle number"
            IS2_atl11_corr_attrs[ptx][XT]['cycle_number']['source'] = "ATL06"
            IS2_atl11_corr_attrs[ptx][XT]['cycle_number']['description'] = ("Cycle number for the "
                "crossing data. Number of 91-day periods that have elapsed since ICESat-2 entered "
                "the science orbit. Each of the 1,387 reference ground track (RGTs) is targeted "
                "in the polar regions once every 91 days.")
            #-- delta time of the crossing track
            IS2_atl11_corr[ptx][XT]['delta_time'] = delta_time['XT'].copy()
            IS2_atl11_fill[ptx][XT]['delta_time'] = delta_time['XT'].fill_value
            IS2_atl11_dims[ptx][XT]['delta_time'] = ['ref_pt']
            IS2_atl11_corr_attrs[ptx][XT]['delta_time'] = {}
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['units'] = "seconds since 2018-01-01"
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['long_name'] = "Elapsed GPS seconds"
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['standard_name'] = "time"
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['calendar'] = "standard"
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['source'] = "ATL06"
            IS2_atl11_corr_attrs[ptx][XT]['delta_time']['description'] = ("Number of GPS "
                "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
                "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
                "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
                "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
                "time in gps_seconds relative to the GPS epoch can be computed.")
            IS2_atl11_corr_attrs[ptx]['delta_time']['coordinates'] = \
                "ref_pt latitude longitude"
            #-- latitude of the crossover measurement
            IS2_atl11_corr[ptx][XT]['latitude'] = latitude['XT'].copy()
            IS2_atl11_fill[ptx][XT]['latitude'] = latitude['XT'].fill_value
            IS2_atl11_dims[ptx][XT]['latitude'] = ['ref_pt']
            IS2_atl11_corr_attrs[ptx][XT]['latitude'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['units'] = "degrees_north"
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['contentType'] = "physicalMeasurement"
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['long_name'] = "crossover latitude"
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['standard_name'] = "latitude"
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['source'] = "ATL06"
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['description'] = ("Center latitude of "
                "selected segments")
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['valid_min'] = -90.0
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['valid_max'] = 90.0
            IS2_atl11_corr_attrs[ptx][XT]['latitude']['coordinates'] = \
                "ref_pt delta_time longitude"
            #-- longitude of the crossover measurement
            IS2_atl11_corr[ptx][XT]['longitude'] = longitude['XT'].copy()
            IS2_atl11_fill[ptx][XT]['longitude'] = longitude['XT'].fill_value
            IS2_atl11_dims[ptx][XT]['longitude'] = ['ref_pt']
            IS2_atl11_corr_attrs[ptx][XT]['longitude'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['units'] = "degrees_east"
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['contentType'] = "physicalMeasurement"
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['long_name'] = "crossover longitude"
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['standard_name'] = "longitude"
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['source'] = "ATL06"
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['description'] = ("Center longitude of "
                "selected segments")
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['valid_min'] = -180.0
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['valid_max'] = 180.0
            IS2_atl11_corr_attrs[ptx][XT]['longitude']['coordinates'] = \
                "ref_pt delta_time latitude"

            #-- interpolated sea level at the crossover measurement
            sea_level = dict(mdt=MDT['XT'],sla=SLA['XT'],adt=ADT['XT'])
            for key,val in sea_level.items():
                #-- add to output
                IS2_atl11_corr[ptx][XT][key] = val.copy()
                IS2_atl11_fill[ptx][XT][key] = val.fill_value
                IS2_atl11_dims[ptx][XT][key] = ['ref_pt']
                IS2_atl11_corr_attrs[ptx][XT][key] = collections.OrderedDict()
                IS2_atl11_corr_attrs[ptx][XT][key]['units'] = "meters"
                IS2_atl11_corr_attrs[ptx][XT][key]['contentType'] = "referenceInformation"
                IS2_atl11_corr_attrs[ptx][XT][key]['long_name'] = attrib[key]['long_name']
                IS2_atl11_corr_attrs[ptx][XT][key]['description'] = attrib[key]['description']
                IS2_atl11_corr_attrs[ptx][XT][key]['source'] = 'AVISO/Copernicus'
                IS2_atl11_corr_attrs[ptx][XT][key]['reference'] = attrib[key]['reference']
                IS2_atl11_corr_attrs[ptx][XT][key]['coordinates'] = \
                    "ref_pt delta_time latitude longitude"

    #-- output HDF5 files with interpolated sea level data
    fargs = (PRD,'AVISO_SEA_LEVEL',TRK,GRAN,SCYC,ECYC,RL,VERS,AUX)
    file_format = '{0}_{1}_{2}{3}_{4}{5}_{6}_{7}{8}.h5'
    output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
    #-- print file information
    logging.info('\t{0}'.format(output_file))
    HDF5_ATL11_corr_write(IS2_atl11_corr, IS2_atl11_corr_attrs,
        CLOBBER=True, INPUT=os.path.basename(FILE), CROSSOVERS=CROSSOVERS,
        FILL_VALUE=IS2_atl11_fill, DIMENSIONS=IS2_atl11_dims,
        FILENAME=output_file)
    #-- change the permissions mode
    os.chmod(output_file, MODE)

#-- PURPOSE: outputting the correction values for ICESat-2 data to HDF5
def HDF5_ATL11_corr_write(IS2_atl11_corr, IS2_atl11_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CROSSOVERS=False,
    CLOBBER=False):
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
    for k,v in IS2_atl11_corr['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam pair
    pairs = [k for k in IS2_atl11_corr.keys() if bool(re.match(r'pt\d',k))]
    for ptx in pairs:
        fileID.create_group(ptx)
        h5[ptx] = {}
        #-- add HDF5 group attributes for beam
        for att_name in ['description','beam_pair','ReferenceGroundTrack',
            'first_cycle','last_cycle','equatorial_radius','polar_radius']:
            fileID[ptx].attrs[att_name] = IS2_atl11_attrs[ptx][att_name]

        #-- ref_pt, cycle number, geolocation and delta_time variables
        for k in ['ref_pt','cycle_number','delta_time','latitude','longitude']:
            #-- values and attributes
            v = IS2_atl11_corr[ptx][k]
            attrs = IS2_atl11_attrs[ptx][k]
            fillvalue = FILL_VALUE[ptx][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}'.format(ptx,k)
            if fillvalue:
                h5[ptx][k] = fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
            else:
                h5[ptx][k] = fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, compression='gzip')
            #-- create or attach dimensions for HDF5 variable
            if DIMENSIONS[ptx][k]:
                #-- attach dimensions
                for i,dim in enumerate(DIMENSIONS[ptx][k]):
                    h5[ptx][k].dims[i].attach_scale(h5[ptx][dim])
            else:
                #-- make dimension
                h5[ptx][k].make_scale(k)
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[ptx][k].attrs[att_name] = att_val

        #-- add to cycle_stats variables
        groups = ['cycle_stats']
        #-- if running crossovers: add to crossing_track_data variables
        if CROSSOVERS:
            groups.append('crossing_track_data')
        for key in groups:
            fileID[ptx].create_group(key)
            h5[ptx][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl11_attrs[ptx][key][att_name]
                fileID[ptx][key].attrs[att_name] = att_val
            for k,v in IS2_atl11_corr[ptx][key].items():
                #-- attributes
                attrs = IS2_atl11_attrs[ptx][key][k]
                fillvalue = FILL_VALUE[ptx][key][k]
                #-- Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}'.format(ptx,key,k)
                if fillvalue:
                    h5[ptx][key][k] = fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
                else:
                    h5[ptx][key][k] = fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, compression='gzip')
                #-- create or attach dimensions for HDF5 variable
                if DIMENSIONS[ptx][key][k]:
                    #-- attach dimensions
                    for i,dim in enumerate(DIMENSIONS[ptx][key][k]):
                        if (key == 'cycle_stats'):
                            h5[ptx][key][k].dims[i].attach_scale(h5[ptx][dim])
                        else:
                            h5[ptx][key][k].dims[i].attach_scale(h5[ptx][key][dim])
                else:
                    #-- make dimension
                    h5[ptx][key][k].make_scale(k)
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[ptx][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Annual Land Ice Height'
    fileID.attrs['summary'] = ('The purpose of ATL11 is to provide an ICESat-2 '
        'satellite cycle summary of heights and height changes of land-based '
        'ice and will be provided as input to ATL15 and ATL16, gridded '
        'estimates of heights and height-changes.')
    fileID.attrs['description'] = ('Land ice parameters for each beam pair. '
        'All parameters are calculated for the same along-track increments '
        'for each beam pair and repeat.')
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
    #-- add attributes for input ATL11 files
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for ptx in pairs:
        lon = IS2_atl11_corr[ptx]['longitude']
        lat = IS2_atl11_corr[ptx]['latitude']
        delta_time = IS2_atl11_corr[ptx]['delta_time']
        valid = np.nonzero(delta_time != FILL_VALUE[ptx]['delta_time'])
        #-- setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time[valid].min() if (delta_time[valid].min() < tmn) else tmn
        tmx = delta_time[valid].max() if (delta_time[valid].max() > tmx) else tmx
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
    atlas_sdp_gps_epoch=IS2_atl11_corr['ancillary_data']['atlas_sdp_gps_epoch']
    gps_seconds = atlas_sdp_gps_epoch + np.array([tmn,tmx])
    #-- calculate leap seconds
    leaps = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
    #-- convert from seconds since 1980-01-06T00:00:00 to Julian days
    MJD = icesat2_toolkit.time.convert_delta_time(gps_seconds - leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date
    YY,MM,DD,HH,MN,SS = icesat2_toolkit.time.convert_julian(MJD + 2400000.5,
        FORMAT='tuple')
    #-- add attributes with measurement date start, end and duration
    tcs = datetime.datetime(int(YY[0]), int(MM[0]), int(DD[0]),
        int(HH[0]), int(MN[0]), int(SS[0]), int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(int(YY[1]), int(MM[1]), int(DD[1]),
        int(HH[1]), int(MN[1]), int(SS[1]), int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    #-- Closing the HDF5 file
    fileID.close()

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Interpolates AVISO sea level anomalies, absolute
            dynamic topography and mean dynamic topography to ICESat-2
            ATL11 annual land ice height data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        icesat2_toolkit.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat-2 ATL11 file to run')
    #-- directory with sea level data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- run with ATL11 crossovers
    parser.add_argument('--crossovers','-C',
        default=False, action='store_true',
        help='Run ATL11 Crossovers')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information about each created file')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files created')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- run for each input ATL11 file
    for FILE in args.infile:
        interp_sea_level_ICESat2(args.directory, FILE,
            CROSSOVERS=args.crossovers, VERBOSE=args.verbose,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()