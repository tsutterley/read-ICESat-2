#!/usr/bin/env python
u"""
reduce_ICESat2_ATL11_raster.py
Written by Tyler Sutterley (06/2022)

Create masks for reducing ICESat-2 ATL11 data using raster imagery

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
    -S X, --sigma X: Standard deviation for Gaussian kernel
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
    read_ICESat2_ATL11.py: reads ICESat-2 annual land ice height data files
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    spatial.py: utilities for reading and writing spatial data
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 06/2022: added option sigma to Gaussian filter raster images
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Written 11/2021
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
import collections
import scipy.ndimage
import scipy.spatial
import scipy.interpolate
import icesat2_toolkit.spatial
import icesat2_toolkit.time
import icesat2_toolkit.utilities
from icesat2_toolkit.convert_delta_time import convert_delta_time
from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11
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

#-- PURPOSE: read ICESat-2 annual land ice height data (ATL11) from NSIDC
#-- reduce to a masked region using raster imagery
def reduce_ICESat2_ATL11_raster(FILE,
    MASK=None,
    FORMAT=None,
    VARIABLES=[],
    OUTPUT=None,
    PROJECTION=None,
    SIGMA=0.0,
    VERBOSE=False,
    MODE=0o775):

    #-- create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- read data from FILE
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = read_HDF5_ATL11(FILE,
        ATTRIBUTES=True)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_'
        r'(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()

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
    #-- gaussian filter mask to increase coverage
    if (SIGMA > 0):
        # convert nan values to 0
        ii,jj = np.nonzero(np.isfinite(dinput['data']))
        dinput['data'] = np.nan_to_num(dinput['data'], nan=0.0)
        dinput['data'] = scipy.ndimage.gaussian_filter(dinput['data'],
            SIGMA, mode='constant', cval=0)
        # return original mask values to true
        dinput['data'][ii,jj] = 1.0
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
    IS2_atl11_mask = {}
    IS2_atl11_fill = {}
    IS2_atl11_dims = {}
    IS2_atl11_mask_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl11_mask['ancillary_data'] = {}
    IS2_atl11_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl11_mask['ancillary_data'][key] = IS2_atl11_mds['ancillary_data'][key]
        #-- Getting attributes of group and included variables
        IS2_atl11_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][key].items():
            IS2_atl11_mask_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam pair within the file
    for ptx in sorted(IS2_atl11_pairs):
        #-- output data dictionaries for beam pair
        IS2_atl11_mask[ptx] = dict(subsetting=collections.OrderedDict())
        IS2_atl11_fill[ptx] = dict(subsetting={})
        IS2_atl11_dims[ptx] = dict(subsetting={})
        IS2_atl11_mask_attrs[ptx] = dict(subsetting={})

        #-- number of average segments and number of included cycles
        #-- shape of along-track data
        n_points,n_cycles = IS2_atl11_mds[ptx]['delta_time'].shape
        #-- along-track latitude, longitude and time
        latitude = np.ma.array(IS2_atl11_mds[ptx]['latitude'],
            fill_value=IS2_atl11_attrs[ptx]['latitude']['_FillValue'])
        latitude.mask = (latitude == latitude.fill_value)
        longitude = np.ma.array(IS2_atl11_mds[ptx]['longitude'],
            fill_value=IS2_atl11_attrs[ptx]['longitude']['_FillValue'])
        longitude.mask = (longitude == longitude.fill_value)
        delta_time = np.ma.array(IS2_atl11_mds[ptx]['delta_time'],
            fill_value=IS2_atl11_attrs[ptx]['delta_time']['_FillValue'])
        delta_time.mask = (delta_time == delta_time.fill_value)

        #-- convert latitude/longitude to raster image projection
        X,Y = transformer.transform(longitude, latitude)

        #-- check where points are within complex hull of triangulation
        #-- or within the bounds of the input raster image
        if v:
            interp_points = np.concatenate((X[:,None],Y[:,None]),axis=1)
            valid = (triangle.find_simplex(interp_points) >= 0)
        else:
            valid = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)

        #-- interpolate raster mask to points
        interp_mask = np.zeros((n_points),dtype=bool)
        #-- skip beam pair interpolation if no data within bounds of raster image
        if np.any(valid):
            interp_mask[valid] = SPL.ev(X[valid], Y[valid])

        #-- group attributes for beam pair
        IS2_atl11_mask_attrs[ptx]['description'] = ('Contains the primary science parameters for this '
            'data set')
        IS2_atl11_mask_attrs[ptx]['beam_pair'] = IS2_atl11_attrs[ptx]['beam_pair']
        IS2_atl11_mask_attrs[ptx]['ReferenceGroundTrack'] = IS2_atl11_attrs[ptx]['ReferenceGroundTrack']
        IS2_atl11_mask_attrs[ptx]['first_cycle'] = IS2_atl11_attrs[ptx]['first_cycle']
        IS2_atl11_mask_attrs[ptx]['last_cycle'] = IS2_atl11_attrs[ptx]['last_cycle']
        IS2_atl11_mask_attrs[ptx]['equatorial_radius'] = IS2_atl11_attrs[ptx]['equatorial_radius']
        IS2_atl11_mask_attrs[ptx]['polar_radius'] = IS2_atl11_attrs[ptx]['polar_radius']

        #-- geolocation, time and reference point
        #-- reference point
        IS2_atl11_mask[ptx]['ref_pt'] = IS2_atl11_mds[ptx]['ref_pt'].copy()
        IS2_atl11_fill[ptx]['ref_pt'] = None
        IS2_atl11_dims[ptx]['ref_pt'] = None
        IS2_atl11_mask_attrs[ptx]['ref_pt'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['ref_pt']['units'] = "1"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['contentType'] = "referenceInformation"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['long_name'] = "Reference point number"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['description'] = ("The reference point is the "
            "7 digit segment_id number corresponding to the center of the ATL06 data used "
            "for each ATL11 point.  These are sequential, starting with 1 for the first "
            "segment after an ascending equatorial crossing node.")
        IS2_atl11_mask_attrs[ptx]['ref_pt']['coordinates'] = \
            "delta_time latitude longitude"
        #-- cycle_number
        IS2_atl11_mask[ptx]['cycle_number'] = IS2_atl11_mds[ptx]['cycle_number'].copy()
        IS2_atl11_fill[ptx]['cycle_number'] = None
        IS2_atl11_dims[ptx]['cycle_number'] = None
        IS2_atl11_mask_attrs[ptx]['cycle_number'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['cycle_number']['units'] = "1"
        IS2_atl11_mask_attrs[ptx]['cycle_number']['long_name'] = "Orbital cycle number"
        IS2_atl11_mask_attrs[ptx]['cycle_number']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['cycle_number']['description'] = ("Number of 91-day periods "
            "that have elapsed since ICESat-2 entered the science orbit. Each of the 1,387 "
            "reference ground track (RGTs) is targeted in the polar regions once "
            "every 91 days.")
        #-- delta time
        IS2_atl11_mask[ptx]['delta_time'] = delta_time.copy()
        IS2_atl11_fill[ptx]['delta_time'] = delta_time.fill_value
        IS2_atl11_dims[ptx]['delta_time'] = ['ref_pt','cycle_number']
        IS2_atl11_mask_attrs[ptx]['delta_time'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl11_mask_attrs[ptx]['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl11_mask_attrs[ptx]['delta_time']['standard_name'] = "time"
        IS2_atl11_mask_attrs[ptx]['delta_time']['calendar'] = "standard"
        IS2_atl11_mask_attrs[ptx]['delta_time']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl11_mask_attrs[ptx]['delta_time']['coordinates'] = \
            "ref_pt cycle_number latitude longitude"
        #-- latitude
        IS2_atl11_mask[ptx]['latitude'] = latitude.copy()
        IS2_atl11_fill[ptx]['latitude'] = latitude.fill_value
        IS2_atl11_dims[ptx]['latitude'] = ['ref_pt']
        IS2_atl11_mask_attrs[ptx]['latitude'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['latitude']['units'] = "degrees_north"
        IS2_atl11_mask_attrs[ptx]['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_mask_attrs[ptx]['latitude']['long_name'] = "Latitude"
        IS2_atl11_mask_attrs[ptx]['latitude']['standard_name'] = "latitude"
        IS2_atl11_mask_attrs[ptx]['latitude']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['latitude']['description'] = ("Center latitude of "
            "selected segments")
        IS2_atl11_mask_attrs[ptx]['latitude']['valid_min'] = -90.0
        IS2_atl11_mask_attrs[ptx]['latitude']['valid_max'] = 90.0
        IS2_atl11_mask_attrs[ptx]['latitude']['coordinates'] = \
            "ref_pt delta_time longitude"
        #-- longitude
        IS2_atl11_mask[ptx]['longitude'] = longitude.copy()
        IS2_atl11_fill[ptx]['longitude'] = longitude.fill_value
        IS2_atl11_dims[ptx]['longitude'] = ['ref_pt']
        IS2_atl11_mask_attrs[ptx]['longitude'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['longitude']['units'] = "degrees_east"
        IS2_atl11_mask_attrs[ptx]['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl11_mask_attrs[ptx]['longitude']['long_name'] = "Longitude"
        IS2_atl11_mask_attrs[ptx]['longitude']['standard_name'] = "longitude"
        IS2_atl11_mask_attrs[ptx]['longitude']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['longitude']['description'] = ("Center longitude of "
            "selected segments")
        IS2_atl11_mask_attrs[ptx]['longitude']['valid_min'] = -180.0
        IS2_atl11_mask_attrs[ptx]['longitude']['valid_max'] = 180.0
        IS2_atl11_mask_attrs[ptx]['longitude']['coordinates'] = \
            "ref_pt delta_time latitude"

        #-- subsetting variables
        IS2_atl11_mask_attrs[ptx]['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce annual land ice height segments to specific "
            "regions of interest.")
        IS2_atl11_mask_attrs[ptx]['subsetting']['data_rate'] = ("Data within this group "
            "are stored at the average segment rate.")

        #-- interpolated raster mask
        IS2_atl11_mask[ptx]['subsetting']['mask'] = interp_mask.copy()
        IS2_atl11_fill[ptx]['subsetting']['mask'] = None
        IS2_atl11_dims[ptx]['subsetting']['mask'] = ['ref_pt']
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask']['contentType'] = "referenceInformation"
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask']['long_name'] = 'Mask'
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask']['description'] = ('Mask calculated '
            'using raster image')
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask']['source'] = os.path.basename(MASK)
        IS2_atl11_mask_attrs[ptx]['subsetting']['mask']['coordinates'] = \
            "../ref_pt ../delta_time ../latitude ../longitude"

    #-- use default output file name and path
    if OUTPUT:
        output_file = os.path.expanduser(OUTPUT)
    else:
        fargs = (PRD,'MASK',TRK,GRAN,SCYC,ECYC,RL,VERS,AUX)
        file_format = '{0}_{1}_{2}{3}_{4}{5}_{6}_{7}{8}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
    #-- print file information
    logging.info('\t{0}'.format(output_file))
    #-- write to output HDF5 file
    HDF5_ATL11_mask_write(IS2_atl11_mask, IS2_atl11_mask_attrs,
        CLOBBER=True, INPUT=os.path.basename(FILE),
        FILL_VALUE=IS2_atl11_fill, DIMENSIONS=IS2_atl11_dims,
        FILENAME=output_file)
    #-- change the permissions mode
    os.chmod(output_file, MODE)

#-- PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL11_mask_write(IS2_atl11_mask, IS2_atl11_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CLOBBER=True):
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
    for k,v in IS2_atl11_mask['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl11_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam pair
    pairs = [k for k in IS2_atl11_mask.keys() if bool(re.match(r'pt\d',k))]
    for ptx in pairs:
        fileID.create_group(ptx)
        h5[ptx] = {}
        #-- add HDF5 group attributes for beam pair
        for att_name in ['description','beam_pair','ReferenceGroundTrack',
            'first_cycle','last_cycle','equatorial_radius','polar_radius']:
            fileID[ptx].attrs[att_name] = IS2_atl11_attrs[ptx][att_name]

        #-- ref_pt, cycle number, geolocation and delta_time variables
        for k in ['ref_pt','cycle_number','delta_time','latitude','longitude']:
            #-- values and attributes
            v = IS2_atl11_mask[ptx][k]
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

        #-- add to subsetting variables
        fileID[ptx].create_group('subsetting')
        h5[ptx]['subsetting'] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl11_attrs[ptx]['subsetting'][att_name]
            fileID[ptx]['subsetting'].attrs[att_name] = att_val
        for k,v in IS2_atl11_mask[ptx]['subsetting'].items():
            #-- attributes
            attrs = IS2_atl11_attrs[ptx]['subsetting'][k]
            fillvalue = FILL_VALUE[ptx]['subsetting'][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(ptx,'subsetting',k)
            if fillvalue:
                h5[ptx]['subsetting'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                    compression='gzip')
            else:
                h5[ptx]['subsetting'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, compression='gzip')
            #-- attach dimensions
            for i,dim in enumerate(DIMENSIONS[ptx]['subsetting'][k]):
                h5[ptx]['subsetting'][k].dims[i].attach_scale(h5[ptx][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[ptx]['subsetting'][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Land Ice Height'
    fileID.attrs['summary'] = ('Subsetting masks and geophysical parameters '
        'for land ice segments needed to interpret and assess the quality '
        'of the height estimates.')
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
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for ptx in pairs:
        lon = IS2_atl11_mask[ptx]['longitude']
        lat = IS2_atl11_mask[ptx]['latitude']
        delta_time = IS2_atl11_mask[ptx]['delta_time']
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
    #-- convert start and end time from ATLAS SDP seconds into UTC time
    time_utc = convert_delta_time(np.array([tmn,tmx]))
    #-- convert to calendar date
    YY,MM,DD,HH,MN,SS = icesat2_toolkit.time.convert_julian(time_utc['julian'],
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
        description="""Create masks for reducing ICESat-2 ATL11 annual land
            ice height data using raster imagery
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        icesat2_toolkit.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL11 file to run')
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
    #-- Gaussian filter raster image to increase coverage
    parser.add_argument('--sigma','-S',
        type=float, default=0.0,
        help='Standard deviation for Gaussian kernel')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- run raster mask program with parameters
    reduce_ICESat2_ATL11_raster(args.file,
        MASK=args.raster,
        FORMAT=args.format,
        VARIABLES=args.variables,
        PROJECTION=args.projection,
        SIGMA=args.sigma,
        OUTPUT=args.output,
        VERBOSE=args.verbose,
        MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
