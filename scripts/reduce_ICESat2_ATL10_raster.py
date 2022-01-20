#!/usr/bin/env python
u"""
reduce_ICESat2_ATL10_raster.py
Written by Tyler Sutterley (12/2021)

Create masks for reducing ICESat-2 ATL10 data using raster imagery

COMMAND LINE OPTIONS:
    -R X, --raster X: Input raster file
    -F X, --format X: Input raster file format
        netCDF4
        HDF5
        geotiff
    -v X, --variables X: variable names of data in HDF5 or netCDF4 file
        x, y and data variable names
    -P X, --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -O X, --output X: Output mask file name
    -V, --verbose: Output information about each created file
    -M X, --mode X: Permission mode of directories and files created

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    read_ICESat2_ATL10.py: reads ICESat-2 sea ice freeboard data files
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    spatial.py: utilities for reading and writing spatial data
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Written 12/2021
"""
from __future__ import print_function

import sys
import os
import re
import h5py
import pyproj
import logging
import argparse
import datetime
import warnings
import numpy as np
import scipy.spatial
import scipy.interpolate
import icesat2_toolkit.spatial
import icesat2_toolkit.time
import icesat2_toolkit.utilities
from icesat2_toolkit.convert_delta_time import convert_delta_time
from icesat2_toolkit.read_ICESat2_ATL10 import read_HDF5_ATL10
warnings.filterwarnings("ignore")

#-- PURPOSE: try to get the projection information for the input file
def get_projection(attributes, PROJECTION):
    #-- coordinate reference system string from file
    try:
        crs = pyproj.CRS.from_string(attributes['projection'])
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- EPSG projection code
    try:
        crs = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    #-- no projection can be made
    raise pyproj.exceptions.CRSError

#-- PURPOSE: find a valid Delaunay triangulation for coordinates x0 and y0
#-- http://www.qhull.org/html/qhull.htm#options
#-- Attempt 1: standard qhull options Qt Qbb Qc Qz
#-- Attempt 2: rescale and center the inputs with option QbB
#-- Attempt 3: joggle the inputs to find a triangulation with option QJ
#-- if no passing triangulations: exit with empty list
def find_valid_triangulation(x0,y0,max_points=1e6):
    #-- don't attempt triangulation if there are a large number of points
    if (len(x0) > max_points):
        #-- if too many points: set triangle as an empty list
        logging.info('Too many points for triangulation')
        return (None,[])

    #-- Attempt 1: try with standard options Qt Qbb Qc Qz
    #-- Qt: triangulated output, all facets will be simplicial
    #-- Qbb: scale last coordinate to [0,m] for Delaunay triangulations
    #-- Qc: keep coplanar points with nearest facet
    #-- Qz: add point-at-infinity to Delaunay triangulation

    #-- Attempt 2 in case of qhull error from Attempt 1 try Qt Qc QbB
    #-- Qt: triangulated output, all facets will be simplicial
    #-- Qc: keep coplanar points with nearest facet
    #-- QbB: scale input to unit cube centered at the origin

    #-- Attempt 3 in case of qhull error from Attempt 2 try QJ QbB
    #-- QJ: joggle input instead of merging facets
    #-- QbB: scale input to unit cube centered at the origin

    #-- try each set of qhull_options
    points = np.concatenate((x0[:,None],y0[:,None]),axis=1)
    for i,opt in enumerate(['Qt Qbb Qc Qz','Qt Qc QbB','QJ QbB']):
        logging.info('qhull option: {0}'.format(opt))
        try:
            triangle = scipy.spatial.Delaunay(points.data, qhull_options=opt)
        except scipy.spatial.qhull.QhullError:
            pass
        else:
            return (i+1,triangle)

    #-- if still errors: set triangle as an empty list
    return (None,[])

#-- PURPOSE: read ICESat-2 land ice data (ATL10) from NSIDC
#-- reduce to a masked region using raster imagery
def reduce_ICESat2_ATL10_raster(FILE,
    MASK=None,
    FORMAT=None,
    VARIABLES=[],
    OUTPUT=None,
    PROJECTION=None,
    VERBOSE=False,
    MODE=0o775):

    #-- create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- read data from input file
    logging.info('{0} -->'.format(os.path.basename(FILE)))
    IS2_atl10_mds,IS2_atl10_attrs,IS2_atl10_beams = read_HDF5_ATL10(FILE,
        ATTRIBUTES=True)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 sea ice freeboard file name
    rx = re.compile(r'(processed_)?(ATL\d{2})-(\d{2})_(\d{4})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,HMN,YY,MM,DD,HH,MN,SS,TRK,CYCL,SN,RL,VERS,AUX=rx.findall(FILE).pop()

    #-- read raster image for spatial coordinates and data
    dinput = icesat2_toolkit.spatial.from_file(MASK, FORMAT,
        xname=VARIABLES[0], yname=VARIABLES[1], varname=VARIABLES[2])
    #-- raster extents
    xmin,xmax,ymin,ymax = np.copy(dinput['attributes']['extent'])
    #-- check that x and y are strictly increasing
    if (np.sign(dinput['attributes']['spacing'][0]) == -1):
        dinput['x'] = dinput['x'][::-1]
        dinput['data'] = dinput['data'][:,::-1]
    if (np.sign(dinput['attributes']['spacing'][1]) == -1):
        dinput['y'] = dinput['y'][::-1]
        dinput['data'] = dinput['data'][::-1,:]
    #-- find valid points within mask
    indy,indx = np.nonzero(dinput['data'])
    #-- check that input points are within convex hull of valid model points
    gridx,gridy = np.meshgrid(dinput['x'],dinput['y'])
    v,triangle = find_valid_triangulation(gridx[indy,indx],gridy[indy,indx])
    #-- create an interpolator for input raster data
    logging.info('Building Spline Interpolator')
    SPL = scipy.interpolate.RectBivariateSpline(dinput['x'], dinput['y'],
        dinput['data'].T, kx=1, ky=1)

    #-- convert projection from input coordinates (EPSG) to data coordinates
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = get_projection(dinput['attributes'], PROJECTION)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    logging.info(crs2.to_proj4())

    #-- copy variables for outputting to HDF5 file
    IS2_atl10_mask = {}
    IS2_atl10_fill = {}
    IS2_atl10_dims = {}
    IS2_atl10_mask_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl10_mask['ancillary_data'] = {}
    IS2_atl10_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl10_mask['ancillary_data'][key] = IS2_atl10_mds['ancillary_data'][key]
        #-- Getting attributes of group and included variables
        IS2_atl10_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl10_attrs['ancillary_data'][key].items():
            IS2_atl10_mask_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl10_beams):
        #-- output data dictionaries for beam
        IS2_atl10_mask[gtx] = dict(freeboard_beam_segment={},leads={})
        IS2_atl10_fill[gtx] = dict(freeboard_beam_segment={},leads={})
        IS2_atl10_dims[gtx] = dict(freeboard_beam_segment={},leads={})
        IS2_atl10_mask_attrs[gtx] = dict(freeboard_beam_segment={},leads={})

        #-- group attributes for beam
        IS2_atl10_mask_attrs[gtx]['Description'] = IS2_atl10_attrs[gtx]['Description']
        IS2_atl10_mask_attrs[gtx]['atlas_pce'] = IS2_atl10_attrs[gtx]['atlas_pce']
        IS2_atl10_mask_attrs[gtx]['atlas_beam_type'] = IS2_atl10_attrs[gtx]['atlas_beam_type']
        IS2_atl10_mask_attrs[gtx]['groundtrack_id'] = IS2_atl10_attrs[gtx]['groundtrack_id']
        IS2_atl10_mask_attrs[gtx]['atmosphere_profile'] = IS2_atl10_attrs[gtx]['atmosphere_profile']
        IS2_atl10_mask_attrs[gtx]['atlas_spot_number'] = IS2_atl10_attrs[gtx]['atlas_spot_number']
        IS2_atl10_mask_attrs[gtx]['sc_orientation'] = IS2_atl10_attrs[gtx]['sc_orientation']

        #-- group attributes for freeboard_beam_segment
        IS2_atl10_mask_attrs[gtx]['freeboard_beam_segment']['Description'] = ("Contains freeboard "
            "estimate and associated height segment parameters for only the sea ice segments by beam.")
        IS2_atl10_mask_attrs[gtx]['freeboard_beam_segment']['data_rate'] = ("Data within this "
            "group are stored at the freeboard swath segment rate.")
        #-- group attributes for leads
        IS2_atl10_mask_attrs[gtx]['leads']['Description'] = ("Contains parameters relating "
            "to the freeboard values.")
        IS2_atl10_mask_attrs[gtx]['leads']['data_rate'] = ("Data within this "
            "group are stored at the lead index rate.")

        #-- for each ATL10 group
        for group in ['freeboard_beam_segment','leads']:
            #-- number of segments
            val = IS2_atl10_mds[gtx][group]
            n_seg = len(val['delta_time'])

            #-- convert latitude/longitude to raster image projection
            X,Y = transformer.transform(val['longitude'], val['latitude'])

            #-- check where points are within complex hull of triangulation
            #-- or within the bounds of the input raster image
            if v:
                interp_points = np.concatenate((X[:,None],Y[:,None]),axis=1)
                valid = (triangle.find_simplex(interp_points) >= 0)
            else:
                valid = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)

            #-- interpolate raster mask to points
            interp_mask = np.zeros((n_seg),dtype=bool)
            #-- skip beam interpolation if no data within bounds of raster image
            if np.any(valid):
                interp_mask[valid] = SPL.ev(X[valid], Y[valid])

            #-- delta time
            IS2_atl10_mask[gtx][group]['delta_time'] = val['delta_time'].copy()
            IS2_atl10_fill[gtx][group]['delta_time'] = None
            IS2_atl10_dims[gtx][group]['delta_time'] = None
            IS2_atl10_mask_attrs[gtx][group]['delta_time'] = {}
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['units'] = "seconds since 2018-01-01"
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['long_name'] = "Elapsed GPS seconds"
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['standard_name'] = "time"
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['source'] = "telemetry"
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['calendar'] = "standard"
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['description'] = ("Number of "
                "GPS seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch "
                "offset is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS "
                "seconds between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP "
                "epoch. By adding the offset contained within atlas_sdp_gps_epoch to delta time "
                "parameters, the time in gps_seconds relative to the GPS epoch can be computed.")
            IS2_atl10_mask_attrs[gtx][group]['delta_time']['coordinates'] = \
                "latitude longitude"
            #-- latitude
            IS2_atl10_mask[gtx][group]['latitude'] = val['latitude'].copy()
            IS2_atl10_fill[gtx][group]['latitude'] = None
            IS2_atl10_dims[gtx][group]['latitude'] = ['delta_time']
            IS2_atl10_mask_attrs[gtx][group]['latitude'] = {}
            IS2_atl10_mask_attrs[gtx][group]['latitude']['units'] = "degrees_north"
            IS2_atl10_mask_attrs[gtx][group]['latitude']['contentType'] = "physicalMeasurement"
            IS2_atl10_mask_attrs[gtx][group]['latitude']['long_name'] = "Latitude"
            IS2_atl10_mask_attrs[gtx][group]['latitude']['standard_name'] = "latitude"
            IS2_atl10_mask_attrs[gtx][group]['latitude']['description'] = ("Latitude of "
                "segment center")
            IS2_atl10_mask_attrs[gtx][group]['latitude']['valid_min'] = -90.0
            IS2_atl10_mask_attrs[gtx][group]['latitude']['valid_max'] = 90.0
            IS2_atl10_mask_attrs[gtx][group]['latitude']['coordinates'] = \
                "delta_time longitude"
            #-- longitude
            IS2_atl10_mask[gtx][group]['longitude'] = val['longitude'].copy()
            IS2_atl10_fill[gtx][group]['longitude'] = None
            IS2_atl10_dims[gtx][group]['longitude'] = ['delta_time']
            IS2_atl10_mask_attrs[gtx][group]['longitude'] = {}
            IS2_atl10_mask_attrs[gtx][group]['longitude']['units'] = "degrees_east"
            IS2_atl10_mask_attrs[gtx][group]['longitude']['contentType'] = "physicalMeasurement"
            IS2_atl10_mask_attrs[gtx][group]['longitude']['long_name'] = "Longitude"
            IS2_atl10_mask_attrs[gtx][group]['longitude']['standard_name'] = "longitude"
            IS2_atl10_mask_attrs[gtx][group]['longitude']['description'] = ("Longitude of "
                "segment center")
            IS2_atl10_mask_attrs[gtx][group]['longitude']['valid_min'] = -180.0
            IS2_atl10_mask_attrs[gtx][group]['longitude']['valid_max'] = 180.0
            IS2_atl10_mask_attrs[gtx][group]['longitude']['coordinates'] = \
                "delta_time latitude"

            #-- subsetting variables
            IS2_atl10_mask[gtx][group]['subsetting'] = {}
            IS2_atl10_fill[gtx][group]['subsetting'] = {}
            IS2_atl10_dims[gtx][group]['subsetting'] = {}
            IS2_atl10_mask_attrs[gtx][group]['subsetting'] = {}
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['Description'] = ("The "
                "subsetting group contains parameters used to reduce sea ice "
                "segments to specific regions of interest.")
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['data_rate'] = ("Data "
                "within this group are stored at the variable segment rate.")

            #-- output mask to HDF5
            IS2_atl10_mask[gtx][group]['subsetting']['mask'] = interp_mask.copy()
            IS2_atl10_fill[gtx][group]['subsetting']['mask'] = None
            IS2_atl10_dims[gtx][group]['subsetting']['mask'] = ['delta_time']
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask'] = {}
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask']['contentType'] = \
                "referenceInformation"
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask']['long_name'] = 'Mask'
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask']['description'] = \
                'Mask calculated using raster image'
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask']['source'] = \
                os.path.basename(MASK)
            IS2_atl10_mask_attrs[gtx][group]['subsetting']['mask']['coordinates'] = \
                "../delta_time ../latitude ../longitude"

    #-- use default output file name and path
    if OUTPUT:
        output_file = os.path.expanduser(OUTPUT)
    else:
        fargs = (PRD,HMN,'MASK',YY,MM,DD,HH,MN,SS,TRK,CYCL,SN,RL,VERS,AUX)
        file_format = '{0}-{1}_{2}_{3}{4}{5}{6}{7}{8}_{9}{10}{11}_{12}_{13}{14}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
    #-- print file information
    logging.info('\t{0}'.format(output_file))
    #-- write to output HDF5 file
    HDF5_ATL10_mask_write(IS2_atl10_mask, IS2_atl10_mask_attrs,
        CLOBBER=True, INPUT=os.path.basename(FILE),
        FILL_VALUE=IS2_atl10_fill, DIMENSIONS=IS2_atl10_dims,
        FILENAME=output_file)
    #-- change the permissions mode
    os.chmod(output_file, MODE)

#-- PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL10_mask_write(IS2_atl10_mask, IS2_atl10_attrs, INPUT=None,
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
    for k,v in IS2_atl10_mask['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl10_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl10_mask.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        #-- add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl10_attrs[gtx][att_name]
        #-- create freeboard_beam_segment and leads groups
        h5[gtx] = dict(freeboard_beam_segment={},leads={})
        for group in ['freeboard_beam_segment','leads']:
            fileID[gtx].create_group(group)
            for att_name in ['Description','data_rate']:
                att_val = IS2_atl10_attrs[gtx][group][att_name]
                fileID[gtx][group].attrs[att_name] = att_val

            #-- delta_time and geolocation variables
            for k in ['delta_time','latitude','longitude']:
                #-- values and attributes
                v = IS2_atl10_mask[gtx][group][k]
                attrs = IS2_atl10_attrs[gtx][group][k]
                fillvalue = FILL_VALUE[gtx][group][k]
                #-- Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}'.format(gtx,group,k)
                if fillvalue:
                    h5[gtx][group][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                        compression='gzip')
                else:
                    h5[gtx][group][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype, compression='gzip')
                #-- create or attach dimensions for HDF5 variable
                if DIMENSIONS[gtx][group][k]:
                    #-- attach dimensions
                    for i,dim in enumerate(DIMENSIONS[gtx][group][k]):
                        h5[gtx][group][k].dims[i].attach_scale(
                            h5[gtx][group][dim])
                else:
                    #-- make dimension
                    h5[gtx][group][k].make_scale(k)
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx][group][k].attrs[att_name] = att_val

            #-- add to subsetting variables
            key = 'subsetting'
            fileID[gtx][group].create_group(key)
            h5[gtx][group][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl10_attrs[gtx][group][key][att_name]
                fileID[gtx][group][key].attrs[att_name] = att_val
            for k,v in IS2_atl10_mask[gtx][group][key].items():
                #-- attributes
                attrs = IS2_atl10_attrs[gtx][group][key][k]
                fillvalue = FILL_VALUE[gtx][group][key][k]
                #-- Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}/{3}'.format(gtx,group,key,k)
                if fillvalue:
                    h5[gtx][group][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
                else:
                    h5[gtx][group][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, compression='gzip')
                #-- attach dimensions
                for i,dim in enumerate(DIMENSIONS[gtx][group][key][k]):
                    h5[gtx][group][key][k].dims[i].attach_scale(
                        h5[gtx][group][dim])
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx][group][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 L3A Sea Ice Freeboard'
    fileID.attrs['summary'] = ('Subsetting masks for sea ice freeboard segments '
        'needed to interpret and assess the quality of the freeboard estimates.')
    fileID.attrs['description'] = ('The data set (ATL10) contains estimates '
        'of sea ice freeboard, calculated using three different approaches. '
        'Sea ice leads used to establish the reference sea surface and '
        'descriptive statistics used in the height estimates are also provided')
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
    #-- add attributes for input ATL10 file
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        #-- for each ATL10 group
        for group in ['freeboard_beam_segment','leads']:
            lon = IS2_atl10_mask[gtx][group]['longitude']
            lat = IS2_atl10_mask[gtx][group]['latitude']
            delta_time = IS2_atl10_mask[gtx][group]['delta_time']
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
    atlas_sdp_gps_epoch=IS2_atl10_mask['ancillary_data']['atlas_sdp_gps_epoch']
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
    tcs = datetime.datetime(int(YY[0]), int(MM[0]), int(DD[0]),
        int(HH[0]), int(MN[0]), int(SS[0]), int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(int(YY[1]), int(MM[1]), int(DD[1]),
        int(HH[1]), int(MN[1]), int(SS[1]), int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    #-- Closing the HDF5 file
    fileID.close()

#-- Main program that calls reduce_ICESat2_ATL10_raster()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Create masks for reducing ICESat-2 data
            using raster imagery
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        icesat2_toolkit.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL10 file to run')
    #-- use default output file name
    parser.add_argument('--output','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Name and path of output file')
    #-- input raster file and file format
    parser.add_argument('--raster','-R',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Input raster file')
    parser.add_argument('--format','-F',
        type=str, default='geotiff', choices=('netCDF4','HDF5','geotiff'),
        help='Input raster file format')
    #-- variable names of data in HDF5 or netCDF4 file
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['x','y','data'],
        help='Variable names of data in HDF5 or netCDF4 files')
    #-- spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args,_ = parser.parse_known_args()

    #-- run raster mask program with parameters
    reduce_ICESat2_ATL10_raster(args.file,
        MASK=args.raster,
        FORMAT=args.format,
        VARIABLES=args.variables,
        PROJECTION=args.projection,
        OUTPUT=args.output,
        VERBOSE=args.verbose,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
