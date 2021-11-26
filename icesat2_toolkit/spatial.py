#!/usr/bin/env python
u"""
spatial.py
Written by Tyler Sutterley (11/2021)

Utilities for reading and operating on spatial data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    h5py: Pythonic interface to the HDF5 binary data format
        https://www.h5py.org/
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL

UPDATE HISTORY:
    Written 11/2021
"""
import os
import re
import io
import gzip
import uuid
import h5py
import logging
import netCDF4
import numpy as np
import osgeo.gdal, osgeo.osr

def case_insensitive_filename(filename):
    """
    Searches a directory for a filename without case dependence
    """
    #-- check if file presently exists with input case
    if not os.access(os.path.expanduser(filename),os.F_OK):
        #-- search for filename without case dependence
        basename = os.path.basename(filename)
        directory = os.path.dirname(os.path.expanduser(filename))
        f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
        if not f:
            raise IOError('{0} not found in file system'.format(filename))
        filename = os.path.join(directory,f.pop())
    return os.path.expanduser(filename)

def from_file(filename, format, **kwargs):
    """
    Wrapper function for reading data from an input format
    """
    #-- read input file to extract spatial coordinates and data
    if (format == 'netCDF4'):
        dinput = from_netCDF4(filename, **kwargs)
    elif (format == 'HDF5'):
        dinput = from_HDF5(filename, **kwargs)
    elif (format == 'geotiff'):
        dinput = from_geotiff(filename, **kwargs)
    else:
        raise ValueError('Invalid format {0}'.format(format))
    return dinput

def from_netCDF4(filename, **kwargs):
    """
    Read data from a netCDF4 file
    Inputs: full path of input netCDF4 file
    Options:
        netCDF4 file is compressed or streamed from memory
        netCDF4 variable names of x, y, and data
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    kwargs.setdefault('xname','x')
    kwargs.setdefault('yname','y')
    kwargs.setdefault('varname','data')
    #-- read data from netCDF4 file
    #-- Open the NetCDF4 file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read as in-memory (diskless) netCDF4 dataset
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=f.read())
    elif (kwargs['compression'] == 'bytes'):
        #-- read as in-memory (diskless) netCDF4 dataset
        fileID = netCDF4.Dataset(uuid.uuid4().hex,memory=filename.read())
    else:
        #-- read netCDF4 dataset
        fileID = netCDF4.Dataset(case_insensitive_filename(filename), 'r')
    #-- Output NetCDF file information
    logging.info(fileID.filepath())
    logging.info(list(fileID.variables.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            ncattr, = [s for s in fileID.ncattrs() if re.match(attr,s,re.I)]
            dinput['attributes'][attr] = fileID.getncattr(ncattr)
        except (ValueError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','grid_mapping','_FillValue']
    #-- mapping between netCDF4 variable names and output names
    variable_mapping = dict(x=kwargs['xname'],y=kwargs['yname'],
        data=kwargs['varname'])
    #-- for each variable
    for key,nc in variable_mapping.items():
        #-- Getting the data from each NetCDF variable
        dinput[key] = fileID.variables[nc][:]
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                ncattr, = [s for s in fileID.variables[nc].ncattrs()
                    if re.match(attr,s,re.I)]
                dinput['attributes'][key][attr] = \
                    fileID.variables[nc].getncattr(ncattr)
            except (ValueError,AttributeError):
                pass
    #-- get projection information if there is a grid_mapping attribute
    if 'grid_mapping' in dinput['attributes']['data'].keys():
        #-- try getting the attribute
        grid_mapping = dinput['attributes']['data']['grid_mapping']
        for att_name in fileID[grid_mapping].ncattrs():
            dinput['attributes']['crs'][att_name] = \
                fileID.variables[nc].getncattr(ncattr)
        #-- get the spatial projection reference information from wkt
        #-- and overwrite the file-level projection attribute (if existing)
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(dinput['attributes']['crs']['crs_wkt'])
        dinput['attributes']['projection'] = srs.ExportToProj4()
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- add extent and spacing attributes
    xmin,xmax = np.min(dinput['x']),np.max(dinput['x'])
    ymin,ymax = np.min(dinput['y']),np.max(dinput['y'])
    dinput['attributes']['extent'] = (xmin,xmax,ymin,ymax)
    dx = dinput['x'][1] - dinput['x'][0]
    dy = dinput['y'][1] - dinput['y'][0]
    dinput['attributes']['spacing'] = (dx,dy)
    #-- Closing the NetCDF file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_HDF5(filename, **kwargs):
    """
    Read data from a HDF5 file
    Inputs: full path of input HDF5 file
    Options:
        HDF5 file is compressed or streamed from memory
        HDF5 variable names of x, y, and data
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    kwargs.setdefault('xname','x')
    kwargs.setdefault('yname','y')
    kwargs.setdefault('varname','data')
    #-- read data from HDF5 file
    #-- Open the HDF5 file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read gzip compressed file and extract into in-memory file object
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            fid = io.BytesIO(f.read())
        #-- set filename of BytesIO object
        fid.filename = os.path.basename(filename)
        #-- rewind to start of file
        fid.seek(0)
        #-- read as in-memory (diskless) HDF5 dataset from BytesIO object
        fileID = h5py.File(fid, 'r')
    elif (kwargs['compression'] == 'bytes'):
        #-- read as in-memory (diskless) HDF5 dataset
        fileID = h5py.File(filename, 'r')
    else:
        #-- read HDF5 dataset
        fileID = h5py.File(case_insensitive_filename(filename), 'r')
    #-- Output HDF5 file information
    logging.info(fileID.filename)
    logging.info(list(fileID.keys()))
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {}
    #-- get attributes for the file
    for attr in ['title','description','projection']:
        #-- try getting the attribute
        try:
            dinput['attributes'][attr] = fileID.attrs[attr]
        except (KeyError,AttributeError):
            pass
    #-- list of attributes to attempt to retrieve from included variables
    attributes_list = ['description','units','long_name','calendar',
        'standard_name','grid_mapping','_FillValue']
    #-- mapping between HDF5 variable names and output names
    variable_mapping = dict(x=kwargs['xname'],y=kwargs['yname'],
        data=kwargs['varname'])
    #-- for each variable
    for key,h5 in variable_mapping.items():
        #-- Getting the data from each HDF5 variable
        dinput[key] = np.copy(fileID[h5][:])
        #-- get attributes for the included variables
        dinput['attributes'][key] = {}
        for attr in attributes_list:
            #-- try getting the attribute
            try:
                dinput['attributes'][key][attr] = fileID[h5].attrs[attr]
            except (KeyError,AttributeError):
                pass
    #-- get projection information if there is a grid_mapping attribute
    if 'grid_mapping' in dinput['attributes']['data'].keys():
        #-- try getting the attribute
        grid_mapping = dinput['attributes']['data']['grid_mapping']
        for att_name,att_val in fileID[grid_mapping].attrs.items():
            dinput['attributes']['crs'][att_name] = att_val
        #-- get the spatial projection reference information from wkt
        #-- and overwrite the file-level projection attribute (if existing)
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(dinput['attributes']['crs']['crs_wkt'])
        dinput['attributes']['projection'] = srs.ExportToProj4()
    #-- convert to masked array if fill values
    if '_FillValue' in dinput['attributes']['data'].keys():
        dinput['data'] = np.ma.asarray(dinput['data'])
        dinput['data'].fill_value = dinput['attributes']['data']['_FillValue']
        dinput['data'].mask = (dinput['data'].data == dinput['data'].fill_value)
    #-- add extent and spacing attributes
    xmin,xmax = np.min(dinput['x']),np.max(dinput['x'])
    ymin,ymax = np.min(dinput['y']),np.max(dinput['y'])
    dinput['attributes']['extent'] = (xmin,xmax,ymin,ymax)
    dx = dinput['x'][1] - dinput['x'][0]
    dy = dinput['y'][1] - dinput['y'][0]
    dinput['attributes']['spacing'] = (dx,dy)
    #-- Closing the HDF5 file
    fileID.close()
    #-- return the spatial variables
    return dinput

def from_geotiff(filename, **kwargs):
    """
    Read data from a geotiff file
    Inputs: full path of input geotiff file
    Options:
        geotiff file is compressed or streamed from memory
    """
    #-- set default keyword arguments
    kwargs.setdefault('compression',None)
    #-- Open the geotiff file for reading
    if (kwargs['compression'] == 'gzip'):
        #-- read gzip compressed file and extract into memory-mapped object
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        with gzip.open(case_insensitive_filename(filename),'r') as f:
            osgeo.gdal.FileFromMemBuffer(mmap_name, f.read())
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        ds = osgeo.gdal.Open(mmap_name)
    elif (kwargs['compression'] == 'bytes'):
        #-- read as GDAL memory-mapped (diskless) geotiff dataset
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        osgeo.gdal.FileFromMemBuffer(mmap_name, filename.read())
        ds = osgeo.gdal.Open(mmap_name)
    else:
        #-- read geotiff dataset
        ds = osgeo.gdal.Open(case_insensitive_filename(filename))
    #-- print geotiff file if verbose
    logging.info(filename)
    #-- create python dictionary for output variables and attributes
    dinput = {}
    dinput['attributes'] = {c:dict() for c in ['x','y','data']}
    #-- get the spatial projection reference information
    srs = ds.GetSpatialRef()
    dinput['attributes']['projection'] = srs.ExportToProj4()
    dinput['attributes']['wkt'] = srs.ExportToWkt()
    #-- get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    #-- get geotiff info
    info_geotiff = ds.GetGeoTransform()
    dinput['attributes']['spacing'] = (info_geotiff[1],info_geotiff[5])
    #-- calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    dinput['attributes']['extent'] = (xmin,xmax,ymin,ymax)
    #-- x and y pixel center coordinates (converted from upper left)
    dinput['x'] = xmin + info_geotiff[1]/2.0 + np.arange(xsize)*info_geotiff[1]
    dinput['y'] = ymax + info_geotiff[5]/2.0 + np.arange(ysize)*info_geotiff[5]
    #-- read full image with GDAL
    dinput['data'] = ds.ReadAsArray()
    #-- check if image has fill values
    dinput['data'] = np.ma.asarray(dinput['data'])
    dinput['data'].mask = np.zeros_like(dinput['data'],dtype=bool)
    if ds.GetRasterBand(1).GetNoDataValue():
        #-- mask invalid values
        dinput['data'].fill_value = ds.GetRasterBand(1).GetNoDataValue()
        #-- create mask array for bad values
        dinput['data'].mask[:] = (dinput['data'].data == dinput['data'].fill_value)
        #-- set attribute for fill value
        dinput['attributes']['data']['_FillValue'] = dinput['data'].fill_value
    #-- close the dataset
    ds = None
    #-- return the spatial variables
    return dinput

def convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10):
    """
    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Inputs:
        phi1: latitude of input ellipsoid in degrees
        h1: height above input ellipsoid in meters
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid

    Options:
        eps: tolerance to prevent division by small numbers
            and to determine convergence
        itmax: maximum number of iterations to use in Newton-Raphson

    Returns:
        phi2: latitude of output ellipsoid in degrees
        h2: height above output ellipsoid in meters

    References:
        Astronomical Algorithms, Jean Meeus, 1991, Willmann-Bell, Inc.
            pp. 77-82
    """
    if (len(phi1) != len(h1)):
        raise ValueError('phi and h have incompatable dimensions')
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- initialize output arrays
    npts = len(phi1)
    phi2 = np.zeros((npts))
    h2 = np.zeros((npts))
    #-- for each point
    for N in range(npts):
        #-- force phi1 into range -90 <= phi1 <= 90
        if (np.abs(phi1[N]) > 90.0):
            phi1[N] = np.sign(phi1[N])*90.0
        #-- handle special case near the equator
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + a1 - a2
        if (np.abs(phi1[N]) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + a1 - a2
        #-- handle special case near the poles
        #-- phi2 = phi1 (latitudes congruent)
        #-- h2 = h1 + b1 - b2
        elif ((90.0 - np.abs(phi1[N])) < eps):
            phi2[N] = np.copy(phi1[N])
            h2[N] = h1[N] + b1 - b2
        #-- handle case if latitude is within 45 degrees of equator
        elif (np.abs(phi1[N]) <= 45):
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = b2 * b2 - a2 * a2
            k1 = a2 * hpr1cos
            k2 = b2 * hpr1sin
            #-- perform newton-raphson iteration to solve for u2
            #-- cos(u2) will not be close to zero since abs(phi1) <= 45
            for i in range(0, itmax+1):
                cosu2 = np.cos(u2)
                fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(phi2r)
        #-- handle final case where latitudes are between 45 degrees and pole
        else:
            #-- convert phi1 to radians
            phi1r = phi1[N] * np.pi/180.0
            sinphi1 = np.sin(phi1r)
            cosphi1 = np.cos(phi1r)
            #-- prevent division by very small numbers
            cosphi1 = np.copy(eps) if (cosphi1 < eps) else cosphi1
            #-- calculate tangent
            tanphi1 = sinphi1 / cosphi1
            u1 = np.arctan(b1 / a1 * tanphi1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinphi1
            hpr1cos = a1 * np.cos(u1) + h1[N] * cosphi1
            #-- set initial value for u2
            u2 = np.copy(u1)
            #-- setup constants
            k0 = a2 * a2 - b2 * b2
            k1 = b2 * hpr1sin
            k2 = a2 * hpr1cos
            #-- perform newton-raphson iteration to solve for u2
            #-- sin(u2) will not be close to zero since abs(phi1) > 45
            for i in range(0, itmax+1):
                sinu2 = np.sin(u2)
                fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                fu2p =  -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                if (np.abs(fu2p) < eps):
                    i = np.copy(itmax)
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if (np.abs(delta) < eps):
                        i = np.copy(itmax)
            #-- convert latitude to degrees and verify values between +/- 90
            phi2r = np.arctan(a2 / b2 * np.tan(u2))
            phi2[N] = phi2r*180.0/np.pi
            if (np.abs(phi2[N]) > 90.0):
                phi2[N] = np.sign(phi2[N])*90.0
            #-- calculate height
            h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(phi2r)

    #-- return the latitude and height
    return (phi2, h2)

def compute_delta_h(a1, f1, a2, f2, lat):
    """
    Compute difference in elevation for two ellipsoids at a given
        latitude using a simplified empirical equation

    Inputs:
        a1: semi-major axis of input ellipsoid
        f1: flattening of input ellipsoid
        a2: semi-major axis of output ellipsoid
        f2: flattening of output ellipsoid
        lat: array of latitudes in degrees

    Returns:
        delta_h: difference in elevation for two ellipsoids

    Reference:
        J Meeus, Astronomical Algorithms, pp. 77-82 (1991)
    """
    #-- force phi into range -90 <= phi <= 90
    gt90, = np.nonzero((lat < -90.0) | (lat > 90.0))
    lat[gt90] = np.sign(lat[gt90])*90.0
    #-- semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1)*a1
    b2 = (1.0 - f2)*a2
    #-- compute delta_a and delta_b coefficients
    delta_a = a2 - a1
    delta_b = b2 - b1
    #-- compute differences between ellipsoids
    #-- delta_h = -(delta_a * cos(phi)^2 + delta_b * sin(phi)^2)
    phi = lat * np.pi/180.0
    delta_h = -(delta_a*np.cos(phi)**2 + delta_b*np.sin(phi)**2)
    return delta_h

def wrap_longitudes(lon):
    """
    Wraps longitudes to range from -180 to +180

    Inputs:
        lon: longitude (degrees east)
    """
    phi = np.arctan2(np.sin(lon*np.pi/180.0),np.cos(lon*np.pi/180.0))
    #-- convert phi from radians to degrees
    return phi*180.0/np.pi

def to_cartesian(lon,lat,h=0.0,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Converts geodetic coordinates to Cartesian coordinates

    Inputs:
        lon: longitude (degrees east)
        lat: latitude (degrees north)

    Options:
        h: height above ellipsoid (or sphere)
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
            * for spherical coordinates set to radius of the Earth
        flat: ellipsoidal flattening (default: WGS84)
            * for spherical coordinates set to 0
    """
    #-- verify axes
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    #-- fix coordinates to be 0:360
    count = np.count_nonzero(lon < 0)
    if (count != 0):
        lt0, = np.nonzero(lon < 0)
        lon[lt0] += 360.0
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- convert from geodetic latitude to geocentric latitude
    dtr = np.pi/180.0
    #-- geodetic latitude in radians
    latitude_geodetic_rad = lat*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(lon*dtr)
    Y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(lon*dtr)
    Z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    #-- return the cartesian coordinates
    return (X,Y,Z)

def to_sphere(x,y,z):
    """
    Convert from cartesian coordinates to spherical coordinates

    Inputs:
        x,y,z in cartesian coordinates
    """
    #-- calculate radius
    rad = np.sqrt(x**2.0 + y**2.0 + z**2.0)
    #-- calculate angular coordinates
    #-- phi: azimuthal angle
    phi = np.arctan2(y,x)
    #-- th: polar angle
    th = np.arccos(z/rad)
    #-- convert to degrees and fix to 0:360
    lon = 180.0*phi/np.pi
    count = np.count_nonzero(lon < 0)
    if (count != 0):
        lt0 = np.nonzero(lon < 0)
        lon[lt0] = lon[lt0]+360.0
    #-- convert to degrees and fix to -90:90
    lat = 90.0 - (180.0*th/np.pi)
    #-- return latitude, longitude and radius
    return (lon,lat,rad)

def to_geodetic(x,y,z,a_axis=6378137.0,flat=1.0/298.257223563):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using a closed form solution

    Inputs:
        x,y,z in cartesian coordinates

    Options:
        a_axis: semimajor axis of the ellipsoid (default: WGS84)
        flat: ellipsoidal flattening (default: WGS84)

    References:
        J Zhu "Exact conversion of Earth-centered, Earth-fixed
            coordinates to geodetic coordinates"
        Journal of Guidance, Control, and Dynamics,
        16(2), 389--391, 1993
        https://arc.aiaa.org/doi/abs/10.2514/3.21016
    """
    #-- semiminor axis of the WGS84 ellipsoid [m]
    b_axis = (1.0 - flat)*a_axis
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- square of first numerical eccentricity
    e12 = ecc1**2
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- calculate distance
    w = np.sqrt(x**2 + y**2)
    #-- calculate longitude
    lon = np.arctan2(y,x)/dtr
    lat = np.zeros_like(lon)
    h = np.zeros_like(lon)
    if (w == 0):
        #-- special case where w == 0 (exact polar solution)
        h = np.sign(z)*z - b_axis
        lat = 90.0*np.sign(z)
    else:
        #-- all other cases
        l = e12/2.0
        m = (w/a_axis)**2.0
        n = ((1.0-e12)*z/b_axis)**2.0
        i = -(2.0*l**2 + m + n)/2.0
        k = (l**2.0 - m - n)*l**2.0
        q = (1.0/216.0)*(m + n - 4.0*l**2)**3.0 + m*n*l**2.0
        D = np.sqrt((2.0*q - m*n*l**2)*m*n*l**2)
        B = i/3.0 - (q+D)**(1.0/3.0) - (q-D)**(1.0/3.0)
        t = np.sqrt(np.sqrt(B**2-k) - (B+i)/2.0)-np.sign(m-n)*np.sqrt((B-i)/2.0)
        wi = w/(t+l)
        zi = (1.0-e12)*z/(t-l)
        #-- calculate latitude and height
        lat = np.arctan2(zi,((1.0-e12)*wi))/dtr
        h = np.sign(t-1.0+l)*np.sqrt((w-wi)**2.0 + (z-zi)**2.0)
    #-- return latitude, longitude and height
    return (lon,lat,h)

def scale_areas(lat, flat=1.0/298.257223563, ref=70.0):
    """
    Calculates area scaling factors for a polar stereographic projection
    including special case of at the exact pole

    Inputs:
        lat: latitude (degrees north)

    Options:
        flat: ellipsoidal flattening (default: WGS84)
        ref: reference latitude (true scale latitude)

    Returns:
        scale: area scaling factors at input latitudes

    References:
        Snyder, J P (1982) Map Projections used by the U.S. Geological Survey
            Forward formulas for the ellipsoid.  Geological Survey Bulletin
            1532, U.S. Government Printing Office.
        JPL Technical Memorandum 3349-85-101
    """
    #-- convert latitude from degrees to positive radians
    theta = np.abs(lat)*np.pi/180.0
    #-- convert reference latitude from degrees to positive radians
    theta_ref = np.abs(ref)*np.pi/180.0
    #-- square of the eccentricity of the ellipsoid
    #-- ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
    ecc2 = 2.0*flat - flat**2
    #-- eccentricity of the ellipsoid
    ecc = np.sqrt(ecc2)
    #-- calculate ratio at input latitudes
    m = np.cos(theta)/np.sqrt(1.0 - ecc2*np.sin(theta)**2)
    t = np.tan(np.pi/4.0 - theta/2.0)/((1.0 - ecc*np.sin(theta)) / \
        (1.0 + ecc*np.sin(theta)))**(ecc/2.0)
    #-- calculate ratio at reference latitude
    mref = np.cos(theta_ref)/np.sqrt(1.0 - ecc2*np.sin(theta_ref)**2)
    tref = np.tan(np.pi/4.0 - theta_ref/2.0)/((1.0 - ecc*np.sin(theta_ref)) / \
        (1.0 + ecc*np.sin(theta_ref)))**(ecc/2.0)
    #-- distance scaling
    k = (mref/m)*(t/tref)
    kp = 0.5*mref*np.sqrt(((1.0+ecc)**(1.0+ecc))*((1.0-ecc)**(1.0-ecc)))/tref
    #-- area scaling
    scale = np.where(np.isclose(theta,np.pi/2.0),1.0/(kp**2),1.0/(k**2))
    return scale
