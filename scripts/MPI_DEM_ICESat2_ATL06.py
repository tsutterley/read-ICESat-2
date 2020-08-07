#!/usr/bin/env python
u"""
MPI_DEM_ICESat2_ATL06.py
Written by Tyler Sutterley (10/2019)
Determines which digital elevation model tiles to read for a given ATL06 file
Reads 3x3 array of tiles for points within bounding box of central mosaic tile
Interpolates digital elevation model to locations of ICESat-2 ATL06 segments

ArcticDEM 2m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/

REMA 8m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/
    http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/

GIMP 30m digital elevation model tiles computed with nsidc_convert_GIMP_DEM.py
    https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
    --model=X: Set the digital elevation model (REMA, ArcticDEM, GIMP) to run
    -M X, --mode=X: Permission mode of directories and files created
    -V, --verbose: Output information about each created file

REQUIRES MPI PROGRAM
    MPI: standardized and portable message-passing system
        https://www.open-mpi.org/
        http://mpitutorial.com/

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    mpi4py: MPI for Python
        http://pythonhosted.org/mpi4py/
        http://mpi4py.readthedocs.org/en/stable/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/
        http://docs.h5py.org/en/stable/mpi.html
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL/
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    convert_calendar_decimal.py: converts from calendar date to decimal year
    time.py: Utilities for calculating time operations
    utilities: download and management utilities for syncing files

REFERENCES:
    https://www.pgc.umn.edu/guides/arcticdem/data-description/
    https://www.pgc.umn.edu/guides/rema/data-description/
    https://nsidc.org/data/nsidc-0645/versions/1

UPDATE HISTORY:
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: fiona for shapefile read.  pyproj for coordinate conversion
        can set the DEM model manually to use the GIMP DEM. verify DEM is finite
        round DEM fill value when creating mask as some DEM tiles are incorrect
        using date functions paralleling public repository. verify output DEM
    Updated 06/2019: assign ArcticDEM by name attribute.  buffer for sub-tiles
    Updated 05/2019: free up memory from invalid tiles. buffer by more geosegs
        assign ArcticDEM polygons by objectid. beam exist check in a try clause
        include mean vertical residuals and number of ground control points
    Updated 04/2019: buffer DEM tiles using values from neighboring tiles
        check if subsetted beam contains land ice data
    Forked 04/2019 from interp_DEM_triangulated_data.py
    Updated 02/2019: python3 compatibility updates
    Updated 01/2018: updated input regular expression to accept lagrangian files
    Updated 06/2017: use actual geospatial lat/lon min and max in attributes
    Updated 04/2017: updated for new triangulation processing chain
    Written 03/2017
"""
from __future__ import print_function

import sys
import os
import re
import uuid
import h5py
import fiona
import getopt
import pyproj
import tarfile
import datetime
import osgeo.gdal
import numpy as np
from mpi4py import MPI
import scipy.interpolate
from shapely.geometry import MultiPoint, Polygon
from icesat2_toolkit.convert_julian import convert_julian
from icesat2_toolkit.convert_delta_time import convert_delta_time

#-- digital elevation models
elevation_dir = {}
elevation_tile_index = {}
#-- ArcticDEM
elevation_dir['ArcticDEM'] = ['ArcticDEM']
elevation_tile_index['ArcticDEM'] = 'ArcticDEM_Tile_Index_Rel7.zip'
#-- GIMP DEM
elevation_dir['GIMP'] = ['GIMP','30m']
elevation_tile_index['GIMP'] = 'gimpdem_Tile_Index_Rel1.1.zip'
#-- REMA DEM
elevation_dir['REMA'] = ['REMA']
elevation_tile_index['REMA'] = 'REMA_Tile_Index_Rel1.1.zip'

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' --model=X: Set the digital elevation model (REMA,ArcticDEM,GIMP)')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' -V, --verbose\t\tOutput information about each created file\n')

#-- PURPOSE: keep track of MPI threads
def info(rank, size):
    print('Rank {0:d} of {1:d}'.format(rank+1,size))
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: set the DEM model to interpolate based on the input granule
def set_DEM_model(GRANULE):
    if GRANULE in ('10','11','12'):
        DEM_MODEL = 'REMA'
    elif GRANULE in ('02','03','04','05','06'):
        DEM_MODEL = 'ArcticDEM'
    return DEM_MODEL

#-- PURPOSE: read zip file containing index shapefiles for finding DEM tiles
def read_DEM_index(index_file, DEM_MODEL):
    #-- read the compressed shapefile and extract entities
    shape = fiona.open('zip://{0}'.format(os.path.expanduser(index_file)))
    epsg = shape.crs['init']
    #-- extract attribute indice for DEM tile (REMA,GIMP) or name (ArcticDEM)
    if (DEM_MODEL == 'REMA'):
        #-- REMA index file attributes:
        #-- name: DEM mosaic name for tile (file name without suffix)
        #-- tile: DEM tile identifier (IMy_IMx)
        #-- nd_value: fill value for elements with no data
        #-- resolution: DEM horizontal spatial resolution (meters)
        #-- creationda: creation date
        #-- raster: (empty)
        #-- fileurl: link to file on PGC server
        #-- spec_type: specific type (DEM)
        #-- qual: density of scenes within tile (0 to 1)
        #-- reg_src: DEM registration source (ICESat or neighbor align)
        #-- num_gcps: number of ground control points
        #-- meanresz: mean vertical residual (meters)
        #-- active: (1)
        #-- qc: (2)
        #-- rel_ver: release version
        #-- num_comp: number of components
        #-- st_area_sh: tile area (meters^2)
        #-- st_length_: perimeter length of tile (meters)
        field = 'tile'
    elif (DEM_MODEL == 'GIMP'):
        #-- GIMP index file attributes (from make_GIMP_tile_shapefile.py):
        #-- name: DEM mosaic name for tile (file name without suffix)
        #-- tile: DEM tile identifier (IMy_IMx)
        #-- nd_value: fill value for elements with no data
        #-- resolution: DEM horizontal spatial resolution (meters)
        #-- fileurl: link to file on NSIDC server
        #-- spec_type: specific type (DEM)
        #-- reg_src: DEM registration source (ICESat or neighbor align)
        #-- rel_ver: release version
        #-- num_comp: number of components
        #-- st_area_sh: tile area (meters^2)
        #-- st_length_: perimeter length of tile (meters)
        field = 'tile'
    elif (DEM_MODEL == 'ArcticDEM'):
        #-- ArcticDEM index file attributes:
        #-- objectid: DEM tile object identifier for sub-tile
        #-- name: DEM mosaic name for sub-tile (file name without suffix)
        #-- tile: DEM tile identifier (IMy_IMx) (non-unique for sub-tiles)
        #-- nd_value: fill value for elements with no data
        #-- resolution: DEM horizontal spatial resolution (meters)
        #-- creationda: creation date
        #-- raster: (empty)
        #-- fileurl: link to file on PGC server
        #-- spec_type: specific type (DEM)
        #-- qual: density of scenes within tile (0 to 1)
        #-- reg_src: DEM registration source (ICESat or neighbor align)
        #-- num_gcps: number of ground control points
        #-- meanresz: mean vertical residual (meters)
        #-- active: (1)
        #-- qc: (2)
        #-- rel_ver: release version
        #-- num_comp: number of components
        #-- st_area_sh: tile area (meters^2)
        #-- st_length_: perimeter length of tile (meters)
        field = 'name'
    #-- create python dictionary for each polygon object
    poly_dict = {}
    attrs_dict = {}
    #-- extract the entities and assign by tile name
    for i,ent in enumerate(shape.values()):
        #-- tile or name attributes
        if DEM_MODEL in ('REMA','GIMP'):
            tile = str(ent['properties'][field])
        else:
            tile, = re.findall(r'^(\d+_\d+_\d+_\d+)',ent['properties'][field])
        #-- extract attributes and assign by tile
        attrs_dict[tile] = {}
        for key,val in ent['properties'].items():
            attrs_dict[tile][key] = val
        #-- upper-left, upper-right, lower-right, lower-left, upper-left
        ul,ur,lr,ll,ul2 = ent['geometry']['coordinates'].pop()
        #-- tile boundaries
        attrs_dict[tile]['xmin'] = ul[0]
        attrs_dict[tile]['xmax'] = lr[0]
        attrs_dict[tile]['ymin'] = lr[1]
        attrs_dict[tile]['ymax'] = ul[1]
        #-- extract Polar Stereographic coordinates for entity
        x = [ul[0],ur[0],lr[0],ll[0],ul2[0]]
        y = [ul[1],ur[1],lr[1],ll[1],ul2[1]]
        poly_obj = Polygon(list(zip(x,y)))
        #-- Valid Polygon may not possess overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        poly_dict[tile] = poly_obj
    #-- close the file
    shape.close()
    #-- return the dictionaries of polygon objects and attributes
    return (poly_dict,attrs_dict,epsg)

#-- PURPOSE: read DEM tile file from gzipped tar files
def read_DEM_file(elevation_file, nd_value):
    #-- open file with tarfile (read)
    tar = tarfile.open(name=elevation_file, mode='r:gz')
    #-- find dem geotiff file within tar file
    member, = [m for m in tar.getmembers() if re.search(r'dem\.tif',m.name)]
    #-- use GDAL memory-mapped file to read dem
    mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
    osgeo.gdal.FileFromMemBuffer(mmap_name, tar.extractfile(member).read())
    ds = osgeo.gdal.Open(mmap_name)
    #-- read data matrix
    im = ds.GetRasterBand(1).ReadAsArray()
    fill_value = ds.GetRasterBand(1).GetNoDataValue()
    fill_value = 0.0 if (fill_value is None) else fill_value
    #-- get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    #-- create mask for finding invalid values
    mask = np.zeros((ysize,xsize),dtype=np.bool)
    indy,indx = np.nonzero((im == fill_value) | (~np.isfinite(im)) |
        (np.ceil(im) == np.ceil(fill_value)))
    mask[indy,indx] = True
    #-- verify that values are finite by replacing with nd_value
    im[indy,indx] = nd_value
    #-- get geotiff info
    info_geotiff = ds.GetGeoTransform()
    #-- calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    #-- close files
    ds = None
    osgeo.gdal.Unlink(mmap_name)
    tar.close()
    #-- create image x and y arrays
    xi = np.arange(xmin,xmax+info_geotiff[1],info_geotiff[1])
    yi = np.arange(ymax,ymin+info_geotiff[5],info_geotiff[5])
    #-- return values (flip y values to be monotonically increasing)
    return (im[::-1,:],mask[::-1,:],xi,yi[::-1])

#-- PURPOSE: read DEM tile file from gzipped tar files to buffer main tile
def read_DEM_buffer(elevation_file, xlimits, ylimits, nd_value):
    #-- open file with tarfile (read)
    tar = tarfile.open(name=elevation_file, mode='r:gz')
    #-- find dem geotiff file within tar file
    member, = [m for m in tar.getmembers() if re.search(r'dem\.tif',m.name)]
    #-- use GDAL memory-mapped file to read dem
    mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
    osgeo.gdal.FileFromMemBuffer(mmap_name, tar.extractfile(member).read())
    ds = osgeo.gdal.Open(mmap_name)
    #-- get geotiff info
    info_geotiff = ds.GetGeoTransform()
    #-- original image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    #-- reduce input image with GDAL
    #-- Specify offset and rows and columns to read
    xoffset = np.int((xlimits[0] - xmin)/info_geotiff[1])
    yoffset = np.int((ymax - ylimits[1])/np.abs(info_geotiff[5]))
    xcount = np.int((xlimits[1] - xlimits[0])/info_geotiff[1]) + 1
    ycount = np.int((ylimits[1] - ylimits[0])/np.abs(info_geotiff[5])) + 1
    #-- read data matrix
    im = ds.GetRasterBand(1).ReadAsArray(xoffset, yoffset, xcount, ycount)
    fill_value = ds.GetRasterBand(1).GetNoDataValue()
    fill_value = 0.0 if (fill_value is None) else fill_value
    #-- create mask for finding invalid values
    mask = np.zeros((ycount,xcount),dtype=np.bool)
    indy,indx = np.nonzero((im == fill_value) | (~np.isfinite(im)) |
        (np.ceil(im) == np.ceil(fill_value)))
    mask[indy,indx] = True
    #-- verify that values are finite by replacing with nd_value
    im[indy,indx] = nd_value
    #-- reduced x and y limits of image
    xmin_reduced = xmin + xoffset*info_geotiff[1]
    xmax_reduced = xmin + xoffset*info_geotiff[1] + (xcount-1)*info_geotiff[1]
    ymax_reduced = ymax + yoffset*info_geotiff[5]
    ymin_reduced = ymax + yoffset*info_geotiff[5] + (ycount-1)*info_geotiff[5]
    #-- close files
    ds = None
    osgeo.gdal.Unlink(mmap_name)
    tar.close()
    #-- create image x and y arrays
    xi = np.arange(xmin_reduced,xmax_reduced+info_geotiff[1],info_geotiff[1])
    yi = np.arange(ymax_reduced,ymin_reduced+info_geotiff[5],info_geotiff[5])
    #-- return values (flip y values to be monotonically increasing)
    return (im[::-1,:],mask[::-1,:],xi,yi[::-1])

#-- PURPOSE: read ICESat-2 data from NSIDC or MPI_ICESat2_ATL03.py
#-- interpolate DEM data to x and y coordinates
def main():
    #-- start MPI communicator
    comm = MPI.COMM_WORLD

    #-- Read the system arguments listed after the program
    long_options=['help','directory=','model=','mode=','verbose']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:M:V', long_options)

    #-- working data directory for location of DEM files
    base_dir = os.getcwd()
    #-- set the DEM model to run for a given granule (else set automatically)
    DEM_MODEL = None
    #-- verbosity settings
    VERBOSE = False
    #-- permissions mode of the local files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage() if (comm.rank==0) else None
            sys.exit()
        elif opt in ("-D","--directory"):
            base_dir = os.path.expanduser(arg)
        elif opt in ("--model"):
            DEM_MODEL = arg
        elif opt in ("-V","--verbose"):
            #-- output module information for process
            info(comm.rank,comm.size)
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter HDF5 file as system argument
    if not arglist:
        raise IOError('No input file entered as system arguments')
    #-- tilde-expansion of listed input file
    FILE = os.path.expanduser(arglist[0])

    #-- read data from input file
    print('{0} -->'.format(FILE)) if (VERBOSE and (comm.rank==0)) else None
    #-- Open the HDF5 file for reading
    fileID = h5py.File(FILE, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(FILE)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = rx.findall(FILE).pop()

    #-- set the digital elevation model based on ICESat-2 granule
    DEM_MODEL = set_DEM_model(GRAN) if (DEM_MODEL is None) else DEM_MODEL
    #-- regular expression pattern for extracting parameters from ArcticDEM name
    rx1 = re.compile(r'(\d+)_(\d+)_(\d+)_(\d+)_(\d+m)_(.*?)$', re.VERBOSE)
    #-- full path to DEM directory
    elevation_directory=os.path.join(base_dir,*elevation_dir[DEM_MODEL])
    #-- zip file containing index shapefiles for finding DEM tiles
    index_file=os.path.join(elevation_directory,elevation_tile_index[DEM_MODEL])

    #-- read data on rank 0
    if (comm.rank == 0):
        #-- read index file for determining which tiles to read
        tile_dict,tile_attrs,tile_epsg = read_DEM_index(index_file,DEM_MODEL)
    else:
        #-- create empty object for list of shapely objects
        tile_dict = None
        tile_attrs = None
        tile_epsg = None

    #-- Broadcast Shapely polygon objects
    tile_dict = comm.bcast(tile_dict, root=0)
    tile_attrs = comm.bcast(tile_attrs, root=0)
    tile_epsg = comm.bcast(tile_epsg, root=0)

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

    #-- copy variables for outputting to HDF5 file
    IS2_atl06_dem = {}
    IS2_atl06_fill = {}
    IS2_atl06_dem_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl06_dem['ancillary_data'] = {}
    IS2_atl06_dem_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl06_dem['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        IS2_atl06_dem_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl06_dem_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl06_beams):
        #-- output data dictionaries for beam
        IS2_atl06_dem[gtx] = dict(land_ice_segments={})
        IS2_atl06_fill[gtx] = dict(land_ice_segments={})
        IS2_atl06_dem_attrs[gtx] = dict(land_ice_segments={})

        #-- number of segments
        segment_id = fileID[gtx]['land_ice_segments']['segment_id'][:]
        n_seg, = fileID[gtx]['land_ice_segments']['segment_id'].shape
        #-- invalid value
        fv = fileID[gtx]['land_ice_segments']['h_li'].fillvalue

        #-- define indices to run for specific process
        ind = np.arange(comm.Get_rank(), n_seg, comm.Get_size(), dtype=np.int)

        #-- extract delta time
        delta_time = np.ma.array(fileID[gtx]['land_ice_segments']['delta_time'][:],
            mask=(fileID[gtx]['land_ice_segments']['delta_time'][:]==fv),
            fill_value=fv)
        #-- extract lat/lon
        longitude = np.ma.array(fileID[gtx]['land_ice_segments']['longitude'][:],
            mask=(fileID[gtx]['land_ice_segments']['longitude'][:]==fv),
            fill_value=fv)
        latitude = np.ma.array(fileID[gtx]['land_ice_segments']['latitude'][:],
            mask=(fileID[gtx]['land_ice_segments']['latitude'][:]==fv),
            fill_value=fv)
        #-- output interpolated digital elevation model
        distributed_dem = np.ma.zeros((n_seg),fill_value=fv,dtype=np.float32)
        distributed_dem.mask = np.ones((n_seg),dtype=np.bool)
        dem_h = np.ma.zeros((n_seg),fill_value=fv,dtype=np.float32)
        dem_h.mask = np.ones((n_seg),dtype=np.bool)
        #-- convert tile projection from latitude longitude to tile EPSG
        proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
        proj2 = pyproj.Proj("+init={0}".format(tile_epsg))
        X,Y = pyproj.transform(proj1, proj2, longitude, latitude)

        #-- convert reduced x and y to shapely multipoint object
        xy_point = MultiPoint(list(zip(X[ind], Y[ind])))

        #-- create complete masks for each DEM tile
        associated_map = {}
        for key,poly_obj in tile_dict.items():
            #-- create empty intersection map array for distributing
            distributed_map = np.zeros((n_seg),dtype=np.int)
            #-- create empty intersection map array for receiving
            associated_map[key] = np.zeros((n_seg),dtype=np.int)
            #-- finds if points are encapsulated (within tile)
            int_test = poly_obj.intersects(xy_point)
            if int_test:
                #-- extract intersected points
                int_map = list(map(poly_obj.intersects,xy_point))
                int_indices, = np.nonzero(int_map)
                #-- set distributed_map indices to True for intersected points
                distributed_map[ind[int_indices]] = True
            #-- communicate output MPI matrices between ranks
            #-- operation is a logical "or" across the elements.
            comm.Allreduce(sendbuf=[distributed_map, MPI.BOOL], \
                recvbuf=[associated_map[key], MPI.BOOL], op=MPI.LOR)
            distributed_map = None
        #-- wait for all processes to finish calculation
        comm.Barrier()
        #-- find valid tiles and free up memory from invalid tiles
        valid_tiles = [k for k,v in associated_map.items() if v.any()]
        invalid_tiles = sorted(set(associated_map.keys()) - set(valid_tiles))
        for key in invalid_tiles:
            associated_map[key] = None

        #-- group attributes for beam
        IS2_atl06_dem_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl06_dem_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl06_dem_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl06_dem_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl06_dem_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl06_dem_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl06_dem_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        #-- group attributes for land_ice_segments
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        #-- geolocation, time and segment ID
        #-- delta time
        IS2_atl06_dem[gtx]['land_ice_segments']['delta_time'] = delta_time
        IS2_atl06_fill[gtx]['land_ice_segments']['delta_time'] = delta_time.fill_value
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['standard_name'] = "time"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['calendar'] = "standard"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['delta_time']['coordinates'] = \
            "segment_id latitude longitude"
        #-- latitude
        IS2_atl06_dem[gtx]['land_ice_segments']['latitude'] = latitude
        IS2_atl06_fill[gtx]['land_ice_segments']['latitude'] = latitude.fill_value
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['units'] = "degrees_north"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['valid_min'] = -90.0
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['valid_max'] = 90.0
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['latitude']['coordinates'] = \
            "segment_id delta_time longitude"
        #-- longitude
        IS2_atl06_dem[gtx]['land_ice_segments']['longitude'] = longitude
        IS2_atl06_fill[gtx]['land_ice_segments']['longitude'] = longitude.fill_value
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['units'] = "degrees_east"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['valid_min'] = -180.0
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['valid_max'] = 180.0
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['longitude']['coordinates'] = \
            "segment_id delta_time latitude"
        #-- segment ID
        IS2_atl06_dem[gtx]['land_ice_segments']['segment_id'] = segment_id
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id']['units'] = "1"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id']['contentType'] = "referenceInformation"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id']['long_name'] = "Along-track segment ID number"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id']['description'] = ("A 7 digit number "
            "identifying the along-track geolocation segment number.  These are sequential, starting with "
            "1 for the first segment after an ascending equatorial crossing node. Equal to the segment_id for "
            "the second of the two 20m ATL03 segments included in the 40m ATL06 segment")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['segment_id']['coordinates'] = \
            "delta_time latitude longitude"

        #-- dem variables
        IS2_atl06_dem[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['Description'] = ("The dem group "
            "contains the reference digital elevation model and geoid heights.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        #-- subsetting variables
        IS2_atl06_dem[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce land ice segments to specific regions of interest.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        #-- for each valid tile
        for key in valid_tiles:
            #-- output mask to HDF5
            IS2_atl06_dem[gtx]['land_ice_segments']['subsetting'][key] = associated_map[key]
            IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'][key] = None
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key] = {}
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['contentType'] = "referenceInformation"
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['long_name'] = '{0} Mask'.format(key)
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['description'] = ('Name '
                'of DEM tile {0} encapsulating the land ice segments.').format(tile_attrs[key]['tile'])
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['source'] = DEM_MODEL
            #-- add DEM attributes
            if DEM_MODEL in ('REMA','ArcticDEM'):
                IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['sigma_h'] = tile_attrs[key]['meanresz']
                IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['n_gcp'] = tile_attrs[key]['num_gcps']
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['coordinates'] = \
                "../segment_id ../delta_time ../latitude ../longitude"

        #-- read and interpolate DEM to coordinates in parallel
        for t in range(comm.Get_rank(), len(valid_tiles), comm.Get_size()):
            key = valid_tiles[t]
            sub = tile_attrs[key]['tile']
            name = tile_attrs[key]['name']
            #-- read central DEM file (geotiff within gzipped tar file)
            tar = '{0}.tar.gz'.format(name)
            elevation_file = os.path.join(elevation_directory,sub,tar)
            DEM,MASK,xi,yi = read_DEM_file(elevation_file,fv)
            #-- buffer DEM using values from adjacent tiles
            #-- use 400m (10 geosegs and divisible by ArcticDEM and REMA pixels)
            #-- use 1500m for GIMP
            bf = 1500 if (DEM_MODEL == 'GIMP') else 400
            ny,nx = np.shape(DEM)
            dx = np.abs(xi[1]-xi[0]).astype('i')
            dy = np.abs(yi[1]-yi[0]).astype('i')
            #-- new buffered DEM and mask
            d = np.full((ny+2*bf//dy,nx+2*bf//dx),fv,dtype=np.float32)
            m = np.ones((ny+2*bf//dy,nx+2*bf//dx),dtype=np.bool)
            d[bf//dy:-bf//dy,bf//dx:-bf//dx] = DEM.copy()
            m[bf//dy:-bf//dy,bf//dx:-bf//dx] = MASK.copy()
            DEM,MASK = (None,None)
            #-- new buffered image x and y coordinates
            x = (xi[0] - bf) + np.arange((nx+2*bf//dx))*dx
            y = (yi[0] - bf) + np.arange((ny+2*bf//dy))*dy
            #-- min and max of left column, center column, right column
            XL,XC,XR = [[xi[0]-bf,xi[0]-dx],[xi[0],xi[-1]],[xi[-1]+dx,xi[-1]+bf]]
            xlimits = [XL,XL,XL,XC,XC,XR,XR,XR] #-- LLLCCRRR
            #-- min and max of bottom row, middle row, top row
            YB,YM,YT = [[yi[0]-bf,yi[0]-dy],[yi[0],yi[-1]],[yi[-1]+dy,yi[-1]+bf]]
            ylimits = [YB,YM,YT,YB,YT,YB,YM,YT] #-- BMTBTBMT

            #-- buffer using neighbor tiles (REMA/GIMP) or sub-tiles (ArcticDEM)
            if (DEM_MODEL == 'REMA'):
                #-- REMA tiles to read to buffer the image
                IMy,IMx = np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                #-- neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles = [IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] #-- LLLCCRRR
                ytiles = [IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] #-- BMTBTBMT
                for xtl,ytl,xlim,ylim in zip(xtiles,ytiles,xlimits,ylimits):
                    #-- read DEM file (geotiff within gzipped tar file)
                    bkey = '{0:02d}_{1:02d}'.format(ytl,xtl)
                    #-- if buffer file is a valid tile within the DEM
                    #-- if file doesn't exist: will be all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        btar = '{0}.tar.gz'.format(tile_attrs[bkey]['name'])
                        buffer_file = os.path.join(elevation_directory,bkey,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise IOError('{0} not found'.format(buffer_file))
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = np.int((x1[0] - x[0])//dx)
                        xmax = np.int((x1[-1] - x[0])//dx) + 1
                        ymin = np.int((y1[0] - y[0])//dy)
                        ymax = np.int((y1[-1] - y[0])//dy) + 1
                        #-- add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)
            elif (DEM_MODEL == 'GIMP'):
                #-- GIMP tiles to read to buffer the image
                IMx,IMy = np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                #-- neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles = [IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] #-- LLLCCRRR
                ytiles = [IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] #-- BMTBTBMT
                for xtl,ytl,xlim,ylim in zip(xtiles,ytiles,xlimits,ylimits):
                    #-- read DEM file (geotiff within gzipped tar file)
                    bkey = '{0:d}_{1:d}'.format(xtl,ytl)
                    #-- if buffer file is a valid tile within the DEM
                    #-- if file doesn't exist: will be all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        btar = '{0}.tar.gz'.format(tile_attrs[bkey]['name'])
                        buffer_file = os.path.join(elevation_directory,bkey,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise IOError('{0} not found'.format(buffer_file))
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = np.int((x1[0] - x[0])//dx)
                        xmax = np.int((x1[-1] - x[0])//dx) + 1
                        ymin = np.int((y1[0] - y[0])//dy)
                        ymax = np.int((y1[-1] - y[0])//dy) + 1
                        #-- add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)
            elif (DEM_MODEL == 'ArcticDEM'):
                #-- ArcticDEM sub-tiles to read to buffer the image
                #-- extract parameters from tile filename
                IMy,IMx,STx,STy,res,vers = rx1.findall(name).pop()
                IMy,IMx,STx,STy = np.array([IMy,IMx,STx,STy],dtype='i')
                #-- neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                #-- LLLCCRRR
                xtiles = [IMx+(STx-2)//2,IMx+(STx-2)//2,IMx+(STx-2)//2,IMx,IMx,
                    IMx+STx//2,IMx+STx//2,IMx+STx//2]
                xsubtiles = [(STx-2) % 2 + 1,(STx-2) % 2 + 1,(STx-2) % 2 + 1,
                    STx,STx,STx % 2 + 1,STx % 2 + 1,STx % 2 + 1]
                #-- BMTBTBMT
                ytiles = [IMy+(STy-2)//2,IMy,IMy+STy//2,IMy+(STy-2)//2,
                    IMy+STy//2,IMy+(STy-2)//2,IMy,IMy+STy//2]
                ysubtiles = [(STy-2) % 2 + 1,STy,STy % 2 + 1,(STy-2) % 2 + 1,
                    STy % 2 + 1,(STy-2) % 2 + 1,STy,STy % 2 + 1]
                #-- for each buffer tile and sub-tile
                kwargs = (xtiles,ytiles,xsubtiles,ysubtiles,xlimits,ylimits)
                for xtl,ytl,xs,ys,xlim,ylim in zip(*kwargs):
                    #-- read DEM file (geotiff within gzipped tar file)
                    args = (ytl,xtl,xs,ys,res,vers)
                    bkey = '{0:02d}_{1:02d}_{2}_{3}'.format(*args)
                    #-- if buffer file is a valid sub-tile within the DEM
                    #-- if file doesn't exist: all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        btar = '{0}.tar.gz'.format(tile_attrs[bkey]['name'])
                        buffer_file = os.path.join(elevation_directory,bsub,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise IOError('{0} not found'.format(buffer_file))
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = np.int((x1[0] - x[0])//dx)
                        xmax = np.int((x1[-1] - x[0])//dx) + 1
                        ymin = np.int((y1[0] - y[0])//dy)
                        ymax = np.int((y1[-1] - y[0])//dy) + 1
                        #-- add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)

            #-- indices of x and y coordinates within tile
            tile_indices, = np.nonzero(associated_map[key])
            #-- use spline interpolation to calculate DEM values at coordinates
            f1 = scipy.interpolate.RectBivariateSpline(x,y,d.T,kx=1,ky=1)
            f2 = scipy.interpolate.RectBivariateSpline(x,y,m.T,kx=1,ky=1)
            dataout = f1.ev(X[tile_indices],Y[tile_indices])
            maskout = f2.ev(X[tile_indices],Y[tile_indices])
            #-- save DEM to output variables
            distributed_dem.data[tile_indices] = dataout
            distributed_dem.mask[tile_indices] = maskout.astype(np.bool)
            #-- clear DEM and mask variables
            f1,f2,dataout,maskout,d,m = (None,None,None,None,None,None)

        #-- communicate output MPI matrices between ranks
        #-- operations are element summations and logical "and" across elements
        comm.Allreduce(sendbuf=[distributed_dem.data, MPI.FLOAT], \
            recvbuf=[dem_h.data, MPI.FLOAT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[distributed_dem.mask, MPI.BOOL], \
            recvbuf=[dem_h.mask, MPI.BOOL], op=MPI.LAND)
        distributed_dem = None
        #-- wait for all distributed processes to finish for beam
        comm.Barrier()

        #-- output interpolated DEM to HDF5
        dem_h.mask[np.abs(dem_h.data) >= 1e4] = True
        dem_h.data[dem_h.mask] = dem_h.fill_value
        IS2_atl06_dem[gtx]['land_ice_segments']['dem']['dem_h'] = dem_h
        IS2_atl06_fill[gtx]['land_ice_segments']['dem']['dem_h'] = dem_h.fill_value
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['units'] = "meters"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['contentType'] = "referenceInformation"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['long_name'] = "DEM Height"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['description'] = ("Height of the DEM, "
            "interpolated by cubic-spline interpolation in the DEM coordinate system to the segment location.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['source'] = DEM_MODEL
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

    #-- parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and bool(valid_tiles):
        #-- output HDF5 files with output masks
        args = (PRD,DEM_MODEL,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX)
        file_format='{0}_{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
        #-- print file information
        print('\t{0}'.format(file_format.format(*args))) if VERBOSE else None
        HDF5_ATL06_dem_write(IS2_atl06_dem, IS2_atl06_dem_attrs, CLOBBER=True,
            INPUT=os.path.basename(FILE), FILL_VALUE=IS2_atl06_fill,
            FILENAME=os.path.join(DIRECTORY,file_format.format(*args)))
        #-- change the permissions mode
        os.chmod(os.path.join(DIRECTORY,file_format.format(*args)), MODE)
    #-- close the input file
    fileID.close()

#-- PURPOSE: outputting the interpolated DEM data for ICESat-2 data to HDF5
def HDF5_ATL06_dem_write(IS2_atl06_dem, IS2_atl06_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, CLOBBER=True):
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
    for k,v in IS2_atl06_dem['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl06_dem.keys() if bool(re.match(r'gt\d[lr]',k))]
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

        #-- segment_id
        v = IS2_atl06_dem[gtx]['land_ice_segments']['segment_id']
        attrs = IS2_atl06_attrs[gtx]['land_ice_segments']['segment_id']
        #-- Defining the HDF5 dataset variables
        val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments','segment_id')
        h5[gtx]['land_ice_segments']['segment_id'] = fileID.create_dataset(val,
            np.shape(v), data=v, dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in attrs.items():
            h5[gtx]['land_ice_segments']['segment_id'].attrs[att_name] = att_val

        #-- geolocation, time and height variables
        for k in ['latitude','longitude','delta_time']:
            #-- values and attributes
            v = IS2_atl06_dem[gtx]['land_ice_segments'][k]
            attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments',k)
            h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                compression='gzip')
            #-- attach dimensions
            for dim in ['segment_id']:
                h5[gtx]['land_ice_segments'][k].dims.create_scale(
                    h5[gtx]['land_ice_segments'][dim], dim)
                h5[gtx]['land_ice_segments'][k].dims[0].attach_scale(
                    h5[gtx]['land_ice_segments'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][k].attrs[att_name] = att_val

        #-- add to output variables
        for key in ['subsetting','dem']:
            fileID[gtx]['land_ice_segments'].create_group(key)
            h5[gtx]['land_ice_segments'][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name]
                fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
            for k,v in IS2_atl06_dem[gtx]['land_ice_segments'][key].items():
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
                for dim in ['segment_id']:
                    h5[gtx]['land_ice_segments'][key][k].dims.create_scale(
                        h5[gtx]['land_ice_segments'][dim], dim)
                    h5[gtx]['land_ice_segments'][key][k].dims[0].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx]['land_ice_segments'][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Land Ice Height'
    fileID.attrs['summary'] = ('Subsetting masks and geophysical parameters '
        'for ice-sheets segments needed to interpret and assess the quality '
        'of land height estimates.')
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
    fileID.attrs['references'] = 'http://nsidc.org/data/icesat2/data.html'
    fileID.attrs['processing_level'] = '4'
    #-- add attributes for input ATL03 and ATL09 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl06_dem[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl06_dem[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl06_dem[gtx]['land_ice_segments']['delta_time']
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
    #-- convert start and end time from ATLAS SDP seconds into UTC time
    time_utc = convert_delta_time(np.array([tmn,tmx]))
    #-- convert to calendar date with convert_julian.py
    YY,MM,DD,HH,MN,SS = convert_julian(time_utc['julian'],FORMAT='tuple')
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

#-- run main program
if __name__ == '__main__':
    main()
