#!/usr/bin/env python
u"""
reduce_ICESat2_ATL07_raster.py
Written by Tyler Sutterley (05/2024)

Create masks for reducing ICESat-2 ATL07 data using raster imagery

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
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/

PROGRAM DEPENDENCIES:
    io/ATL07.py: reads ICESat-2 sea ice height data files
    spatial.py: utilities for reading and writing spatial data
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2024: use wrapper to importlib for optional dependencies
    Updated 04/2024: use timescale for temporal operations
    Updated 03/2024: use pathlib to define and operate on paths
    Updated 12/2022: single implicit import of altimetry tools
        refactored ICESat-2 data product read programs under io
    Updated 06/2022: added option sigma to Gaussian filter raster images
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Written 11/2021
"""
from __future__ import print_function

import re
import logging
import pathlib
import argparse
import datetime
import warnings
import numpy as np
import scipy.ndimage
import scipy.spatial
import scipy.interpolate
import icesat2_toolkit as is2tk
import timescale.time

# attempt imports
h5py = is2tk.utilities.import_dependency('h5py')
pyproj = is2tk.utilities.import_dependency('pyproj')

# PURPOSE: try to get the projection information for the input file
def get_projection(attributes, PROJECTION):
    # coordinate reference system string from file
    try:
        crs = pyproj.CRS.from_string(attributes['projection'])
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # EPSG projection code
    try:
        crs = pyproj.CRS.from_epsg(int(PROJECTION))
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # no projection can be made
    raise pyproj.exceptions.CRSError

# PURPOSE: find a valid Delaunay triangulation for coordinates x0 and y0
# http://www.qhull.org/html/qhull.htm#options
# Attempt 1: standard qhull options Qt Qbb Qc Qz
# Attempt 2: rescale and center the inputs with option QbB
# Attempt 3: joggle the inputs to find a triangulation with option QJ
# if no passing triangulations: exit with empty list
def find_valid_triangulation(x0, y0, max_points=1e6):
    """
    Attempt to find a valid Delaunay triangulation for coordinates

    - Attempt 1: ``Qt Qbb Qc Qz``
    - Attempt 2: ``Qt Qc QbB``
    - Attempt 3: ``QJ QbB``

    Parameters
    ----------
    x0: float
        x-coordinates
    y0: float
        y-coordinates
    max_points: int or float, default 1e6
        Maximum number of coordinates to attempt to triangulate
    """
    # don't attempt triangulation if there are a large number of points
    if (len(x0) > max_points):
        # if too many points: set triangle as an empty list
        logging.info('Too many points for triangulation')
        return (None,[])

    # Attempt 1: try with standard options Qt Qbb Qc Qz
    # Qt: triangulated output, all facets will be simplicial
    # Qbb: scale last coordinate to [0,m] for Delaunay triangulations
    # Qc: keep coplanar points with nearest facet
    # Qz: add point-at-infinity to Delaunay triangulation

    # Attempt 2 in case of qhull error from Attempt 1 try Qt Qc QbB
    # Qt: triangulated output, all facets will be simplicial
    # Qc: keep coplanar points with nearest facet
    # QbB: scale input to unit cube centered at the origin

    # Attempt 3 in case of qhull error from Attempt 2 try QJ QbB
    # QJ: joggle input instead of merging facets
    # QbB: scale input to unit cube centered at the origin

    # try each set of qhull_options
    points = np.concatenate((x0[:,None],y0[:,None]),axis=1)
    for i,opt in enumerate(['Qt Qbb Qc Qz','Qt Qc QbB','QJ QbB']):
        logging.info(f'qhull option: {opt}')
        try:
            triangle = scipy.spatial.Delaunay(points.data, qhull_options=opt)
        except scipy.spatial.qhull.QhullError:
            pass
        else:
            return (i+1,triangle)

    # if still errors: set triangle as an empty list
    return (None,[])

# PURPOSE: read ICESat-2 sea ice data (ATL07) from NSIDC
# reduce to a masked region using raster imagery
def reduce_ICESat2_ATL07_raster(FILE,
    MASK=None,
    FORMAT=None,
    VARIABLES=[],
    OUTPUT=None,
    PROJECTION=None,
    SIGMA=0.0,
    TOLERANCE=0.5,
    VERBOSE=False,
    MODE=0o775):

    # create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # read data from input file
    FILE = pathlib.Path(FILE).expanduser().absolute()
    logging.info(f'{str(FILE)} -->')
    IS2_atl07_mds,IS2_atl07_attrs,IS2_atl07_beams = \
        is2tk.io.ATL07.read_granule(FILE, ATTRIBUTES=True)
    # extract parameters from ICESat-2 ATLAS HDF5 sea ice file name
    rx = re.compile(r'(processed_)?(ATL\d{2})-(\d{2})_(\d{4})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,HMN,YY,MM,DD,HH,MN,SS,TRK,CYCL,SN,RL,VERS,AUX = \
        rx.findall(FILE.name).pop()

    # read raster image for spatial coordinates and data
    MASK = pathlib.Path(MASK).expanduser().absolute()
    dinput = is2tk.spatial.from_file(MASK, FORMAT,
        xname=VARIABLES[0], yname=VARIABLES[1], varname=VARIABLES[2])
    # raster extents
    xmin,xmax,ymin,ymax = np.copy(dinput['attributes']['extent'])
    # check that x and y are strictly increasing
    if (np.sign(dinput['attributes']['spacing'][0]) == -1):
        dinput['x'] = dinput['x'][::-1]
        dinput['data'] = dinput['data'][:,::-1]
    if (np.sign(dinput['attributes']['spacing'][1]) == -1):
        dinput['y'] = dinput['y'][::-1]
        dinput['data'] = dinput['data'][::-1,:]
    # find valid points within mask
    indy,indx = np.nonzero(dinput['data'])
    # check that input points are within convex hull of valid model points
    gridx,gridy = np.meshgrid(dinput['x'],dinput['y'])
    v,triangle = find_valid_triangulation(gridx[indy,indx],gridy[indy,indx])
    # gaussian filter mask to increase coverage
    if (SIGMA > 0):
        # convert nan values to 0
        dinput['data'] = np.nan_to_num(dinput['data'], nan=0.0)
        ii,jj = np.nonzero(np.logical_not(dinput['data'].mask) &
            (dinput['data'] != 0.0))
        # gaussian filter image
        dinput['data'] = scipy.ndimage.gaussian_filter(dinput['data'],
            SIGMA, mode='constant', cval=0)
        # return original mask values to true
        dinput['data'][ii,jj] = 1.0
    # create an interpolator for input raster data
    logging.info('Building Spline Interpolator')
    SPL = scipy.interpolate.RectBivariateSpline(dinput['x'], dinput['y'],
        dinput['data'].T, kx=1, ky=1)

    # convert projection from input coordinates (EPSG) to data coordinates
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = get_projection(dinput['attributes'], PROJECTION)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    logging.info(crs2.to_proj4())

    # copy variables for outputting to HDF5 file
    IS2_atl07_mask = {}
    IS2_atl07_fill = {}
    IS2_atl07_dims = {}
    IS2_atl07_mask_attrs = {}
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    IS2_atl07_mask['ancillary_data'] = {}
    IS2_atl07_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl07_mask['ancillary_data'][key] = IS2_atl07_mds['ancillary_data'][key]
        # Getting attributes of group and included variables
        IS2_atl07_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in IS2_atl07_attrs['ancillary_data'][key].items():
            IS2_atl07_mask_attrs['ancillary_data'][key][att_name] = att_val

    # for each input beam within the file
    for gtx in sorted(IS2_atl07_beams):
        # output data dictionaries for beam
        IS2_atl07_mask[gtx] = dict(sea_ice_segments={})
        IS2_atl07_fill[gtx] = dict(sea_ice_segments={})
        IS2_atl07_dims[gtx] = dict(sea_ice_segments={})
        IS2_atl07_mask_attrs[gtx] = dict(sea_ice_segments={})

        # number of segments
        val = IS2_atl07_mds[gtx]['sea_ice_segments']
        n_seg = len(val['height_segment_id'])

        # convert latitude/longitude to raster image projection
        X,Y = transformer.transform(val['longitude'], val['latitude'])

        # check where points are within complex hull of triangulation
        # or within the bounds of the input raster image
        if v:
            interp_points = np.concatenate((X[:,None],Y[:,None]),axis=1)
            valid = (triangle.find_simplex(interp_points) >= 0)
        else:
            valid = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)

        # interpolate raster mask to points
        interp_mask = np.zeros((n_seg),dtype=bool)
        # skip beam interpolation if no data within bounds of raster image
        if np.any(valid):
            interp_mask[valid] = (SPL.ev(X[valid], Y[valid]) >= TOLERANCE)

        # group attributes for beam
        IS2_atl07_mask_attrs[gtx]['Description'] = IS2_atl07_attrs[gtx]['Description']
        IS2_atl07_mask_attrs[gtx]['atlas_pce'] = IS2_atl07_attrs[gtx]['atlas_pce']
        IS2_atl07_mask_attrs[gtx]['atlas_beam_type'] = IS2_atl07_attrs[gtx]['atlas_beam_type']
        IS2_atl07_mask_attrs[gtx]['groundtrack_id'] = IS2_atl07_attrs[gtx]['groundtrack_id']
        IS2_atl07_mask_attrs[gtx]['atmosphere_profile'] = IS2_atl07_attrs[gtx]['atmosphere_profile']
        IS2_atl07_mask_attrs[gtx]['atlas_spot_number'] = IS2_atl07_attrs[gtx]['atlas_spot_number']
        IS2_atl07_mask_attrs[gtx]['sc_orientation'] = IS2_atl07_attrs[gtx]['sc_orientation']
        # group attributes for sea_ice_segments
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['Description'] = ("Top group for sea "
            "ice segments as computed by the ATBD algorithm.")
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['data_rate'] = ("Data within this "
            "group are stored at the variable segment rate.")

        # geolocation, time and segment ID
        # delta time
        IS2_atl07_mask[gtx]['sea_ice_segments']['delta_time'] = val['delta_time'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['delta_time'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['delta_time'] = None
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['standard_name'] = "time"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['source'] = "telemetry"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['calendar'] = "standard"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['description'] = ("Number of "
            "GPS seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch "
            "offset is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS "
            "seconds between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP "
            "epoch. By adding the offset contained within atlas_sdp_gps_epoch to delta time "
            "parameters, the time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['delta_time']['coordinates'] = \
            "height_segment_id latitude longitude"
        # latitude
        IS2_atl07_mask[gtx]['sea_ice_segments']['latitude'] = val['latitude'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['latitude'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['latitude'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['units'] = "degrees_north"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['valid_min'] = -90.0
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['valid_max'] = 90.0
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['latitude']['coordinates'] = \
            "height_segment_id delta_time longitude"
        # longitude
        IS2_atl07_mask[gtx]['sea_ice_segments']['longitude'] = val['longitude'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['longitude'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['longitude'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['units'] = "degrees_east"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['valid_min'] = -180.0
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['valid_max'] = 180.0
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['longitude']['coordinates'] = \
            "height_segment_id delta_time latitude"
        # segment ID
        IS2_atl07_mask[gtx]['sea_ice_segments']['height_segment_id'] = val['height_segment_id']
        IS2_atl07_fill[gtx]['sea_ice_segments']['height_segment_id'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['height_segment_id'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id']['units'] = "1"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id']['contentType'] = "referenceInformation"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id']['long_name'] = \
            "Identifier of each height segment"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id']['description'] = \
            "Identifier of each height segment"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['height_segment_id']['coordinates'] = \
            "delta_time latitude longitude"
        # geolocation segment beginning
        IS2_atl07_mask[gtx]['sea_ice_segments']['geoseg_beg'] = val['geoseg_beg'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['geoseg_beg'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['geoseg_beg'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg']['units'] = "1"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg']['contentType'] = "referenceInformation"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg']['long_name'] = "Beginning GEOSEG"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg']['description'] = \
            "Geolocation segment (geoseg) ID associated with the first photon used in this sea ice segment"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_beg']['coordinates'] = \
            "height_segment_id delta_time latitude longitude"
        # geolocation segment ending
        IS2_atl07_mask[gtx]['sea_ice_segments']['geoseg_end'] = val['geoseg_end'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['geoseg_end'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['geoseg_end'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end']['units'] = "1"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end']['contentType'] = "referenceInformation"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end']['long_name'] = "Ending GEOSEG"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end']['description'] = \
            "Geolocation segment (geoseg) ID associated with the last photon used in this sea ice segment"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['geoseg_end']['coordinates'] = \
            "height_segment_id delta_time latitude longitude"
        # along track distance
        IS2_atl07_mask[gtx]['sea_ice_segments']['seg_dist_x'] = val['seg_dist_x'].copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['seg_dist_x'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['seg_dist_x'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x']['units'] = "meters"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x']['contentType'] = "referenceInformation"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x']['long_name'] = "Along track distance"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x']['description'] = \
            "Along-track distance from the equator crossing to the segment center."
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['seg_dist_x']['coordinates'] = \
            "height_segment_id delta_time latitude longitude"

        # subsetting variables
        IS2_atl07_mask[gtx]['sea_ice_segments']['subsetting'] = {}
        IS2_atl07_fill[gtx]['sea_ice_segments']['subsetting'] = {}
        IS2_atl07_dims[gtx]['sea_ice_segments']['subsetting'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce sea ice segments to specific regions of interest.")
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['data_rate'] = ("Data within this "
            "group are stored at the variable segment rate.")

        # output mask to HDF5
        IS2_atl07_mask[gtx]['sea_ice_segments']['subsetting']['mask'] = interp_mask.copy()
        IS2_atl07_fill[gtx]['sea_ice_segments']['subsetting']['mask'] = None
        IS2_atl07_dims[gtx]['sea_ice_segments']['subsetting']['mask'] = ['delta_time']
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask'] = {}
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['contentType'] = \
            "referenceInformation"
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['long_name'] = 'Mask'
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['description'] = \
            'Mask calculated using raster image'
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['source'] = MASK.name
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['sigma'] = SIGMA
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['tolerance'] = TOLERANCE
        IS2_atl07_mask_attrs[gtx]['sea_ice_segments']['subsetting']['mask']['coordinates'] = \
            "../height_segment_id ../delta_time ../latitude ../longitude"

    # use default output file name and path
    if OUTPUT:
        output_file = pathlib.Path(OUTPUT).expanduser().absolute()
    else:
        fargs = (PRD,HMN,'MASK',YY,MM,DD,HH,MN,SS,TRK,CYCL,SN,RL,VERS,AUX)
        file_format = '{0}-{1}_{2}_{3}{4}{5}{6}{7}{8}_{9}{10}{11}_{12}_{13}{14}.h5'
        output_file = FILE.with_name(file_format.format(*fargs))
    # print file information
    logging.info(f'\t{str(output_file)}')
    # write to output HDF5 file
    HDF5_ATL07_mask_write(IS2_atl07_mask, IS2_atl07_mask_attrs,
        CLOBBER=True, INPUT=FILE.name,
        FILL_VALUE=IS2_atl07_fill, DIMENSIONS=IS2_atl07_dims,
        FILENAME=output_file)
    # change the permissions mode
    output_file.chmod(mode=MODE)

# PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL07_mask_write(IS2_atl07_mask, IS2_atl07_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CLOBBER=False):
    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    FILENAME = pathlib.Path(FILENAME).expanduser().absolute()
    fileID = h5py.File(FILENAME, clobber)

    # create HDF5 records
    h5 = {}

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k,v in IS2_atl07_mask['ancillary_data'].items():
        # Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl07_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # write each output beam
    beams = [k for k in IS2_atl07_mask.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        # add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl07_attrs[gtx][att_name]
        # create sea_ice_segments group
        fileID[gtx].create_group('sea_ice_segments')
        h5[gtx] = dict(sea_ice_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl07_attrs[gtx]['sea_ice_segments'][att_name]
            fileID[gtx]['sea_ice_segments'].attrs[att_name] = att_val

        # delta_time, geolocation and segment identification variables
        for k in ['delta_time','latitude','longitude','height_segment_id',
            'geoseg_beg','geoseg_end','seg_dist_x']:
            # values and attributes
            v = IS2_atl07_mask[gtx]['sea_ice_segments'][k]
            attrs = IS2_atl07_attrs[gtx]['sea_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['sea_ice_segments'][k]
            # Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'sea_ice_segments',k)
            if fillvalue:
                h5[gtx]['sea_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                    compression='gzip')
            else:
                h5[gtx]['sea_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, compression='gzip')
            # create or attach dimensions for HDF5 variable
            if DIMENSIONS[gtx]['sea_ice_segments'][k]:
                # attach dimensions
                for i,dim in enumerate(DIMENSIONS[gtx]['sea_ice_segments'][k]):
                    h5[gtx]['sea_ice_segments'][k].dims[i].attach_scale(
                        h5[gtx]['sea_ice_segments'][dim])
            else:
                # make dimension
                h5[gtx]['sea_ice_segments'][k].make_scale(k)
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['sea_ice_segments'][k].attrs[att_name] = att_val

        # add to subsetting variables
        key = 'subsetting'
        fileID[gtx]['sea_ice_segments'].create_group(key)
        h5[gtx]['sea_ice_segments'][key] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl07_attrs[gtx]['sea_ice_segments'][key][att_name]
            fileID[gtx]['sea_ice_segments'][key].attrs[att_name] = att_val
        for k,v in IS2_atl07_mask[gtx]['sea_ice_segments'][key].items():
            # attributes
            attrs = IS2_atl07_attrs[gtx]['sea_ice_segments'][key][k]
            fillvalue = FILL_VALUE[gtx]['sea_ice_segments'][key][k]
            # Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}/{3}'.format(gtx,'sea_ice_segments',key,k)
            if fillvalue:
                h5[gtx]['sea_ice_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
            else:
                h5[gtx]['sea_ice_segments'][key][k] = \
                    fileID.create_dataset(val, np.shape(v), data=v,
                    dtype=v.dtype, compression='gzip')
            # attach dimensions
            for i,dim in enumerate(DIMENSIONS[gtx]['sea_ice_segments'][key][k]):
                h5[gtx]['sea_ice_segments'][key][k].dims[i].attach_scale(
                    h5[gtx]['sea_ice_segments'][dim])
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['sea_ice_segments'][key][k].attrs[att_name] = att_val

    # HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 L3A Sea Ice Height'
    fileID.attrs['summary'] = ('Subsetting masks for sea ice segments needed '
        'to interpret and assess the quality of the height estimates.')
    fileID.attrs['description'] = ('The data set (ATL07) contains along-track '
        'heights for sea ice and open water leads (at varying length scales) '
        'relative to the WGS84 ellipsoid (ITRF2014 reference frame) after '
        'adjustment for geoidal and tidal variations, and inverted barometer '
        'effects.')
    date_created = datetime.datetime.today()
    fileID.attrs['date_created'] = date_created.isoformat()
    project = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = project
    platform = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = platform
    # add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    # add attributes for input ATL07 file
    fileID.attrs['input_files'] = pathlib.Path(INPUT).name
    # find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl07_mask[gtx]['sea_ice_segments']['longitude']
        lat = IS2_atl07_mask[gtx]['sea_ice_segments']['latitude']
        delta_time = IS2_atl07_mask[gtx]['sea_ice_segments']['delta_time']
        # setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time.min() if (delta_time.min() < tmn) else tmn
        tmx = delta_time.max() if (delta_time.max() > tmx) else tmx
    # add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    # convert start and end time from ATLAS SDP seconds into timescale
    ts = timescale.time.Timescale().from_deltatime(np.array([tmn,tmx]),
        epoch=timescale.time._atlas_sdp_epoch, standard='GPS')
    dt = np.datetime_as_string(ts.to_datetime(), unit='s')
    # add attributes with measurement date start, end and duration
    fileID.attrs['time_coverage_start'] = str(dt[0])
    fileID.attrs['time_coverage_end'] = str(dt[1])
    fileID.attrs['time_coverage_duration'] = f'{tmx-tmn:0.0f}'
    # add software information
    fileID.attrs['software_reference'] = is2tk.version.project_name
    fileID.attrs['software_version'] = is2tk.version.full_version
    # Closing the HDF5 file
    fileID.close()


# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Create masks for reducing ICESat-2 data
            using raster imagery
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('file',
        type=pathlib.Path,
        help='ICESat-2 ATL07 file to run')
    # use default output file name
    parser.add_argument('--output','-O',
        type=pathlib.Path,
        help='Name and path of output file')
    # input raster file and file format
    parser.add_argument('--raster','-R',
        type=pathlib.Path,
        help='Input raster file')
    parser.add_argument('--format','-F',
        type=str, default='geotiff', choices=('netCDF4','HDF5','geotiff'),
        help='Input raster file format')
    # variable names of data in HDF5 or netCDF4 file
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['x','y','data'],
        help='Variable names of data in HDF5 or netCDF4 files')
    # spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    # Gaussian filter raster image to increase coverage
    parser.add_argument('--sigma','-S',
        type=float, default=0.0,
        help='Standard deviation for Gaussian kernel')
    # tolerance in interpolated mask to set as valid
    parser.add_argument('--tolerance','-T',
        type=float, default=0.5,
        help='Tolerance to set as valid mask')
    # verbosity settings
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # run raster mask program with parameters
    reduce_ICESat2_ATL07_raster(args.file,
        MASK=args.raster,
        FORMAT=args.format,
        VARIABLES=args.variables,
        PROJECTION=args.projection,
        SIGMA=args.sigma,
        TOLERANCE=args.tolerance,
        OUTPUT=args.output,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
