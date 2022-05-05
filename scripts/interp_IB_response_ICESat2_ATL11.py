#!/usr/bin/env python
u"""
interp_IB_response_ICESat2_ATL11.py
Written by Tyler Sutterley (05/2022)
Calculates and interpolates inverse-barometer responses to times and
    locations of ICESat-2 ATL11 annual land ice height data
    This data will be interpolated for all valid points
    (masking land values will be needed for accurate assessments)

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -R X, --reanalysis X: Reanalysis model to run
        ERA-Interim: http://apps.ecmwf.int/datasets/data/interim-full-moda
        ERA5: http://apps.ecmwf.int/data-catalogues/era5/?class=ea
        MERRA-2: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
    -m X, --mean X: Start and end year range for mean
    -d X, --density X: Density of seawater in kg/m^3
    -C, --crossovers: Run ATL11 Crossovers
    -V, --verbose: Output information about each created file
    -M X, --mode X: Permission mode of directories and files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    read_ICESat2_ATL11.py: reads ICESat-2 annual land ice height data files
    time.py: utilities for calculating time operations
    utilities.py: download and management utilities for syncing files
    calc_delta_time.py: calculates difference between universal and dynamic time

REFERENCES:
    C Wunsch and D Stammer, Atmospheric loading and the oceanic "inverted
        barometer" effect, Reviews of Geophysics, 35(1), 79--107, (1997).
        https://doi.org/10.1029/96RG03037
    P S Callahan, TOPEX/POSEIDON Project GDR Users Handbook, JPL Doc. D-8944,
        Rev. A, 84 pp., (1994)

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 03/2021: spatially subset sea level pressure maps to conserve memory
        additionally calculate conventional IB response using an average MSLP
        replaced numpy bool/int to prevent deprecation warnings
    Written 02/2021
"""
from __future__ import print_function

import os
import re
import h5py
import pyproj
import logging
import netCDF4
import argparse
import datetime
import operator
import itertools
import numpy as np
import collections
import scipy.interpolate
import icesat2_toolkit.time
import icesat2_toolkit.utilities
from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11

#-- PURPOSE: compress complete list values into a set of ranges
def compress_list(i,n):
    for _,b in itertools.groupby(enumerate(i), lambda v: ((v[1]-v[0])//n)*n):
        group = list(map(operator.itemgetter(1),b))
        yield (group[0], group[-1])

#-- PURPOSE: read land sea mask to get indices of oceanic values
def ncdf_landmask(FILENAME,MASKNAME,OCEAN):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        landsea = np.squeeze(fileID.variables[MASKNAME][:].copy())
    return (landsea == OCEAN)

#-- PURPOSE: read reanalysis mean sea level pressure
def ncdf_mean_pressure(FILENAME,VARNAME,LONNAME,LATNAME):
    with netCDF4.Dataset(FILENAME,'r') as fileID:
        #-- extract pressure and remove singleton dimensions
        mean_pressure = np.array(fileID.variables[VARNAME][:].squeeze())
        longitude = fileID.variables[LONNAME][:].squeeze()
        latitude = fileID.variables[LATNAME][:].squeeze()
    return (mean_pressure,longitude,latitude)

#-- PURPOSE: find pressure files in a directory
def find_pressure_files(ddir, MODEL, MJD):
    #-- regular expression pattern for finding files
    if (MODEL == 'ERA-Interim'):
        regex_pattern = r'ERA\-Interim\-Hourly\-MSL\-({0})\.nc$'
        joiner = r'\-'
    elif (MODEL == 'ERA5'):
        regex_pattern = r'ERA5\-Hourly\-MSL\-({0})\.nc$'
        joiner = r'\-'
    elif (MODEL == 'MERRA-2'):
        regex_pattern = r'MERRA2_\d{{3}}.tavg1_2d_slv_Nx.({0}).(.*?).nc$'
        joiner = r''
    #-- list of dates to read
    dates = []
    #-- for each unique Modified Julian Day (MJD)
    for mjd in np.unique(np.floor(MJD)):
        #-- append day prior, day of and day after
        JD = mjd + np.arange(-1,2) + 2400000.5
        #-- convert from Julian Days to calendar dates
        Y,M,D,_,_,_ = icesat2_toolkit.time.convert_julian(JD,
            ASTYPE=int, FORMAT='tuple')
        #-- append day as formatted strings
        for y,m,d in zip(Y,M,D):
            dates.append(joiner.join([str(y),str(m).zfill(2),str(d).zfill(2)]))
    #-- compile regular expression pattern for finding dates
    rx = re.compile(regex_pattern.format('|'.join(dates)))
    flist = [os.path.join(ddir,f) for f in os.listdir(ddir) if rx.match(f)]
    #-- return the sorted list of unique files
    return sorted(set(flist))

#-- PURPOSE: read sea level pressure fields and calculate anomalies
#-- spatially subset pressure maps to latitudinal range of ATL11 points
def ncdf_pressure(FILENAMES,VARNAME,TIMENAME,LATNAME,MEAN,OCEAN,INDICES,AREA):
    #-- shape of subsetted pressure field
    ny,nx = np.shape(MEAN[INDICES,:])
    nfiles = len(FILENAMES)
    #-- allocate for pressure fields
    SLP = np.ma.zeros((24*nfiles,ny,nx))
    TPX = np.ma.zeros((24*nfiles,ny,nx))
    MJD = np.zeros((24*nfiles))
    #-- calculate total area of reanalysis ocean
    #-- ocean pressure points will be based on reanalysis mask
    ii,jj = np.nonzero(OCEAN)
    ocean_area = np.sum(AREA[ii,jj])
    #-- parameters for conventional TOPEX/POSEIDON IB correction
    rho0 = 1025.0
    g0 = -9.80665
    p0 = 101325.0
    #-- counter for filling arrays
    c = 0
    #-- for each file
    for FILENAME in FILENAMES:
        with netCDF4.Dataset(FILENAME,'r') as fileID:
            #-- extract coordinates
            latitude = fileID.variables[LATNAME][:].squeeze()
            #-- convert time to Modified Julian Days
            delta_time = np.copy(fileID.variables[TIMENAME][:])
            units = fileID.variables[TIMENAME].units
            epoch,to_secs = icesat2_toolkit.time.parse_date_string(units)
            for t,dt in enumerate(delta_time):
                MJD[c] = icesat2_toolkit.time.convert_delta_time(dt*to_secs,
                    epoch1=epoch, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
                #-- check dimensions for expver slice
                if (fileID.variables[VARNAME].ndim == 4):
                    _,nexp,_,_ = fileID.variables[VARNAME].shape
                    #-- sea level pressure for time
                    pressure = fileID.variables[VARNAME][t,:,:,:].copy()
                    #-- iterate over expver slices to find valid outputs
                    for j in range(nexp):
                        #-- check if any are valid for expver
                        if np.any(pressure[j,:,:]):
                            #-- remove average with respect to time
                            AveRmvd = pressure[j,:,:] - MEAN
                            #-- conventional TOPEX/POSEIDON IB correction
                            TPX[c,:,:] = (pressure[j,INDICES,:]-p0)/(rho0*g0)
                            break
                else:
                    #-- sea level pressure for time
                    pressure = fileID.variables[VARNAME][t,:,:].copy()
                    #-- remove average with respect to time
                    AveRmvd = pressure - MEAN
                    #-- conventional TOPEX/POSEIDON IB correction
                    TPX[c,:,:] = (pressure[INDICES,:]-p0)/(rho0*g0)
                #-- calculate average oceanic pressure values
                AVERAGE = np.sum(AveRmvd[ii,jj]*AREA[ii,jj])/ocean_area
                #-- calculate sea level pressure anomalies and
                #-- reduce to latitudinal range of ATL11 points
                SLP[c,:,:] = AveRmvd[INDICES,:] - AVERAGE
                #-- clear temp variables for iteration to free up memory
                pressure,AveRmvd = (None,None)
                #-- add to counter
                c += 1
    #-- verify latitudes are sorted in ascending order
    LAT = latitude[INDICES]
    ilat = np.argsort(LAT)
    SLP = SLP[:,ilat,:]
    TPX = TPX[:,ilat,:]
    LAT = LAT[ilat]
    #-- verify time is sorted in ascending order
    itime = np.argsort(MJD)
    SLP = SLP[itime,:,:]
    TPX = TPX[itime,:,:]
    MJD = MJD[itime]
    #-- return the sea level pressure anomalies, latitudes and times
    return (SLP,TPX,LAT,MJD)

#-- PURPOSE: read ICESat-2 annual land ice height data (ATL11) from NSIDC
#-- calculate and interpolate the instantaneous inverse barometer response
def interp_IB_response_ICESat2(base_dir, FILE, MODEL, RANGE=None,
    DENSITY=None, CROSSOVERS=False, VERBOSE=False, MODE=0o775):

    #-- create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- directory setup for reanalysis model
    ddir = os.path.join(base_dir,MODEL)
    #-- set model specific parameters
    if (MODEL == 'ERA-Interim'):
        #-- mean sea level pressure file
        input_mean_file = 'ERA-Interim-Mean-MSL-{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA-Interim-Invariant-Parameters.nc'
        VARNAME = 'msl'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        #-- projection string
        proj4_params = ('+proj=longlat +ellps=WGS84 +datum=WGS84 '
            '+no_defs lon_wrap=180')
    elif (MODEL == 'ERA5'):
        #-- mean sea level pressure file
        input_mean_file = 'ERA5-Mean-MSL-{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'ERA5-Invariant-Parameters.nc'
        VARNAME = 'msl'
        LONNAME = 'longitude'
        LATNAME = 'latitude'
        TIMENAME = 'time'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'lsm'
        OCEAN = 0
        #-- projection string
        proj4_params = ('+proj=longlat +ellps=WGS84 +datum=WGS84 '
            '+no_defs lon_wrap=180')
    elif (MODEL == 'MERRA-2'):
        #-- mean sea level pressure file
        input_mean_file = 'MERRA2.Mean_SLP.{0:4d}-{1:4d}.nc'
        #-- input land-sea mask for ocean redistribution
        input_mask_file = 'MERRA2_101.const_2d_asm_Nx.00000000.nc4'
        VARNAME = 'SLP'
        LONNAME = 'lon'
        LATNAME = 'lat'
        TIMENAME = 'time'
        #-- land-sea mask variable name and value of oceanic points
        MASKNAME = 'FROCEAN'
        OCEAN = 1
        #-- projection string
        proj4_params = 'epsg:4326'

    #-- read data from FILE
    logging.info('{0} -->'.format(os.path.basename(FILE)))
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = read_HDF5_ATL11(FILE,
        ATTRIBUTES=True, CROSSOVERS=CROSSOVERS)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_'
        r'(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = IS2_atl11_mds['ancillary_data']['atlas_sdp_gps_epoch']
    #-- HDF5 group name for across-track data
    XT = 'crossing_track_data'

    #-- read mean pressure field
    mean_file = os.path.join(ddir,input_mean_file.format(RANGE[0],RANGE[1]))
    mean_pressure,lon,lat=ncdf_mean_pressure(mean_file,VARNAME,LONNAME,LATNAME)

    #-- pyproj transformer for converting from input coordinates (EPSG)
    #-- to model coordinates
    crs1 = pyproj.CRS.from_string('epsg:{0:d}'.format(4326))
    crs2 = pyproj.CRS.from_string(proj4_params)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    #-- grid step size in radians
    dphi = np.pi*np.abs(lon[1] - lon[0])/180.0
    dth = np.pi*np.abs(lat[1] - lat[0])/180.0
    #-- calculate meshgrid from latitude and longitude
    gridlon,gridlat = np.meshgrid(lon,lat)
    gridphi = gridlon*np.pi/180.0
    #-- calculate colatitude
    gridtheta = (90.0 - gridlat)*np.pi/180.0

    #-- ellipsoidal parameters of WGS84 ellipsoid
    #-- semimajor axis of the ellipsoid [m]
    a_axis = 6378137.0
    #-- flattening of the ellipsoid
    flat = 1.0/298.257223563
    #-- semiminor axis of the ellipsoid [m]
    b_axis = (1.0 -flat)*a_axis
    #-- calculate grid areas globally
    AREA = dphi*dth*np.sin(gridtheta)*np.sqrt((a_axis**2)*(b_axis**2) *
        ((np.sin(gridtheta)**2)*(np.cos(gridphi)**2) +
        (np.sin(gridtheta)**2)*(np.sin(gridphi)**2)) +
        (a_axis**4)*(np.cos(gridtheta)**2))
    #-- read land-sea mask to find ocean values
    #-- ocean pressure points will be based on reanalysis mask
    MASK = ncdf_landmask(os.path.join(ddir,input_mask_file),MASKNAME,OCEAN)

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
        #-- dictionary with output inverse barometer variables
        IB,TPX = ({},{})
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
        #-- along-track (AT) inverse barometer corrections
        IB['AT'] = np.ma.zeros((n_points,n_cycles),fill_value=fv)
        IB['AT'].mask = np.copy(delta_time['AT'].mask)
        #-- along-track (AT) conventional TOPEX/POSEIDON IB corrections
        TPX['AT'] = np.ma.zeros((n_points,n_cycles),fill_value=fv)
        TPX['AT'].mask = np.copy(delta_time['AT'].mask)
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
            #-- across-track (XT) inverse barometer corrections
            IB['XT'] = np.ma.zeros((n_cross),fill_value=fv)
            IB['XT'].mask = np.copy(delta_time['XT'].mask)
            #-- across-track (AT) conventional TOPEX/POSEIDON IB corrections
            TPX['XT'] = np.ma.zeros((n_cross),fill_value=fv)
            TPX['XT'].mask = np.copy(delta_time['XT'].mask)

        #-- calculate corrections for along-track and across-track data
        for track in groups:
            #-- convert time from ATLAS SDP to Modified Julian Days
            gps_seconds = atlas_sdp_gps_epoch + delta_time[track]
            leap_seconds = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
            MJD = icesat2_toolkit.time.convert_delta_time(gps_seconds-leap_seconds,
                epoch1=(1980,1,6,0,0,0),epoch2=(1858,11,17,0,0,0),scale=1.0/86400.0)

            #-- calculate projected coordinates of input coordinates
            ix,iy = transformer.transform(longitude[track], latitude[track])
            #-- reduce range of SLP values to be read for a given file
            #-- buffer by grid spacing to prevent edge effects in interpolation
            latrange = np.zeros((2))
            latrange[0] = np.min(iy) - 4.0*(dth*180.0/np.pi)
            latrange[1] = np.max(iy) + 4.0*(dth*180.0/np.pi)
            #-- truncate to valid latitudinal values of input grid
            np.clip(latrange, np.min(lat), np.max(lat), out=latrange)
            #-- find grid points within bounds of ATL11 beam pair
            indices, = np.nonzero((lat >= latrange[0]) & (lat <= latrange[1]))

            #-- colatitudes of the ATL11 measurements
            th = (90.0 - latitude[track])*np.pi/180.0
            #-- gravitational acceleration at mean sea level at the equator
            ge = 9.780356
            #-- gravitational acceleration at mean sea level over colatitudes
            #-- from Heiskanen and Moritz, Physical Geodesy, (1967)
            gs = ge*(1.0 + 5.2885e-3*np.cos(th)**2 - 5.9e-6*np.cos(2.0*th)**2)

            #-- calculate sea level corrections for track type
            if (track == 'AT'):
                #-- calculate for each cycle if along-track
                for cycle in range(n_cycles):
                    #-- find valid points points
                    if np.all(delta_time[track].mask[:,cycle]):
                        IB[track].mask[:,cycle] = True
                        continue
                    valid, = np.nonzero(np.logical_not(delta_time[track].mask[:,cycle]))
                    #-- find each reanalysis pressure field for cycle
                    FILENAMES = find_pressure_files(ddir,MODEL,MJD[valid,cycle])
                    n_files = len(FILENAMES)
                    if (n_files == 0):
                        IB[track].mask[valid] = True
                        continue
                    #-- read sea level pressure and calculate anomalies
                    islp,itpx,ilat,imjd = ncdf_pressure(FILENAMES,VARNAME,TIMENAME,
                        LATNAME,mean_pressure,MASK,indices,AREA)
                    #-- create an interpolator for mean sea level pressure anomalies
                    R1 = scipy.interpolate.RegularGridInterpolator((imjd,ilat,lon),
                        islp, bounds_error=False)
                    R2 = scipy.interpolate.RegularGridInterpolator((imjd,ilat,lon),
                        itpx, bounds_error=False)
                    SLP = R1.__call__(np.c_[MJD[valid,cycle],iy[valid],ix[valid]])
                    #-- calculate inverse barometer response
                    IB[track].data[valid,cycle] = -SLP*(DENSITY*gs[valid])**-1
                    TPX[track].data[valid,cycle] = R2.__call__(np.c_[MJD[valid,cycle],
                        iy[valid],ix[valid]])
                    #-- clear variables for iteration to free up memory
                    islp,itpx,R1,R2 = (None,None,None,None)
            elif (track == 'XT'):
                #-- find valid points to points
                if np.all(delta_time[track].mask):
                    IB[track].mask[:] = True
                    continue
                #-- compress list of Modified Julian Days into ranges
                #-- limits the number of files to be read for a given date
                #-- and the number of times a file needs to be read
                MJD_list = compress_list(np.unique(np.floor(MJD)),10)
                #-- for each date range
                for imin,imax in MJD_list:
                    #-- find where valid and within the range of MJD
                    valid, = np.nonzero(np.logical_not(delta_time[track].mask) &
                        (MJD >= imin) & (MJD <= (imax+1)))
                    #-- find each reanalysis pressure field
                    FILENAMES = find_pressure_files(ddir,MODEL,np.arange(imin,imax+1))
                    n_files = len(FILENAMES)
                    if (n_files == 0):
                        IB[track].mask[valid] = True
                        continue
                    #-- read sea level pressure and calculate anomalies
                    islp,itpx,ilat,imjd = ncdf_pressure(FILENAMES,VARNAME,TIMENAME,
                        LATNAME,mean_pressure,MASK,indices,AREA)
                    #-- create an interpolator for mean sea level pressure anomalies
                    R1 = scipy.interpolate.RegularGridInterpolator((imjd,ilat,lon),
                        islp, bounds_error=False)
                    R2 = scipy.interpolate.RegularGridInterpolator((imjd,ilat,lon),
                        itpx, bounds_error=False)
                    SLP = R1.__call__(np.c_[MJD[valid],iy[valid],ix[valid]])
                    #-- calculate inverse barometer response
                    IB[track][valid] = -SLP*(DENSITY*gs[valid])**-1
                    TPX[track].data[valid] = R2.__call__(np.c_[MJD[valid],
                        iy[valid],ix[valid]])
                    #-- clear variables for iteration to free up memory
                    islp,itpx,R1,R2 = (None,None,None,None)
            #-- replace any nan values with fill value
            IB[track].mask |= np.isnan(IB[track].data)
            TPX[track].mask |= np.isnan(TPX[track].data)
            IB[track].data[IB[track].mask] = IB[track].fill_value
            TPX[track].data[IB[track].mask] = TPX[track].fill_value

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

        #-- inverse barometer response
        IS2_atl11_corr[ptx]['cycle_stats']['ib'] = IB['AT'].copy()
        IS2_atl11_fill[ptx]['cycle_stats']['ib'] = IB['AT'].fill_value
        IS2_atl11_dims[ptx]['cycle_stats']['ib'] = ['ref_pt','cycle_number']
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['units'] = "meters"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['contentType'] = "referenceInformation"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['long_name'] = "inverse barometer"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['description'] = ("Instantaneous inverse "
            "barometer effect due to atmospheric loading")
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['source'] = MODEL
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['reference'] = \
            'https://doi.org/10.1029/96RG03037'
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['ib']['coordinates'] = \
            "../ref_pt ../cycle_number ../delta_time ../latitude ../longitude"

        #-- conventional (TOPEX/POSEIDON) inverse barometer response
        IS2_atl11_corr[ptx]['cycle_stats']['tpx'] = TPX['AT'].copy()
        IS2_atl11_fill[ptx]['cycle_stats']['tpx'] = TPX['AT'].fill_value
        IS2_atl11_dims[ptx]['cycle_stats']['tpx'] = ['ref_pt','cycle_number']
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx'] = collections.OrderedDict()
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['units'] = "meters"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['contentType'] = "referenceInformation"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['long_name'] = "inverse barometer"
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['description'] = ("Conventional "
            "(TOPEX/POSEIDON) instantaneous inverse barometer effect due to "
            "atmospheric loading")
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['source'] = MODEL
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['reference'] = \
            ' TOPEX/POSEIDON Project GDR Users Handbook'
        IS2_atl11_corr_attrs[ptx]['cycle_stats']['tpx']['coordinates'] = \
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

            #-- inverse barometer response at the crossover measurement
            IS2_atl11_corr[ptx][XT]['ib'] = IB['XT']
            IS2_atl11_fill[ptx][XT]['ib'] = IB['XT'].fill_value
            IS2_atl11_dims[ptx][XT]['ib'] = ['ref_pt']
            IS2_atl11_corr_attrs[ptx][XT]['ib'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['ib']['units'] = "meters"
            IS2_atl11_corr_attrs[ptx][XT]['ib']['contentType'] = "referenceInformation"
            IS2_atl11_corr_attrs[ptx][XT]['ib']['long_name'] = "inverse barometer"
            IS2_atl11_corr_attrs[ptx][XT]['ib']['description'] = ("Instantaneous inverse "
                "barometer effect due to atmospheric loading")
            IS2_atl11_corr_attrs[ptx][XT]['ib']['source'] = MODEL
            IS2_atl11_corr_attrs[ptx][XT]['ib']['reference'] = \
                'https://doi.org/10.1029/96RG03037'
            IS2_atl11_corr_attrs[ptx][XT]['ib']['coordinates'] = \
                "ref_pt delta_time latitude longitude"

            #-- conventional (TOPEX/POSEIDON) inverse barometer response
            IS2_atl11_corr[ptx][XT]['tpx'] = TPX['XT'].copy()
            IS2_atl11_fill[ptx][XT]['tpx'] = TPX['XT'].fill_value
            IS2_atl11_dims[ptx][XT]['tpx'] = ['ref_pt']
            IS2_atl11_corr_attrs[ptx][XT]['tpx'] = collections.OrderedDict()
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['units'] = "meters"
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['contentType'] = "referenceInformation"
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['long_name'] = "inverse barometer"
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['description'] = ("Conventional "
                "(TOPEX/POSEIDON) instantaneous inverse barometer effect due to "
                "atmospheric loading")
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['source'] = MODEL
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['reference'] = \
                'TOPEX/POSEIDON Project GDR Users Handbook'
            IS2_atl11_corr_attrs[ptx][XT]['tpx']['coordinates'] = \
                "ref_pt delta_time latitude longitude"

    #-- output HDF5 files with interpolated inverse barometer data
    fargs = (PRD,MODEL,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX)
    file_format = '{0}_{1}_IB_{2}{3}_{4}{5}_{6}_{7}{8}.h5'
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
        description="""Calculates and interpolates inverse-barometer
            responses to times and locations of ICESat-2 ATL11 annual
            land ice height data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        icesat2_toolkit.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat-2 ATL11 file to run')
    #-- directory with reanalysis data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    choices = ['ERA-Interim','ERA5','MERRA-2']
    parser.add_argument('--reanalysis','-R',
        metavar='REANALYSIS', type=str,
        default='ERA5', choices=choices,
        help='Reanalysis Model')
    #-- start and end years to run for mean
    parser.add_argument('--mean','-m',
        metavar=('START','END'), type=int, nargs=2,
        default=[2000,2020],
        help='Start and end year range for mean')
    #-- ocean fluidic density [kg/m^3]
    parser.add_argument('--density','-d',
        metavar='RHO', type=float, default=1030.0,
        help='Density of seawater in kg/m^3')
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
    args,_ = parser.parse_known_args()()

    #-- run for each input ATL11 file
    for FILE in args.infile:
        interp_IB_response_ICESat2(args.directory, FILE, args.reanalysis,
            RANGE=args.mean, DENSITY=args.density, CROSSOVERS=args.crossovers,
            VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
