#!/usr/bin/env python
u"""
MPI_DEM_ICESat2_ATL03.py
Written by Tyler Sutterley (05/2021)
Determines which digital elevation model tiles to read for a given ATL03 file
Reads 3x3 array of tiles for points within bounding box of central mosaic tile
Interpolates digital elevation model to ICESat-2 ATL03 photon event locations

ArcticDEM 2m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/

REMA 8m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/
    http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/

GIMP 30m digital elevation model tiles computed with nsidc_convert_GIMP_DEM.py
    https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/

COMMAND LINE OPTIONS:
    -D X, --directory X: Working data directory
    -m X, --model X: Digital elevation model (REMA, ArcticDEM, GIMP) to run
    -V, --verbose: Output information about each created file
    -M X, --mode X: Permission mode of directories and files created

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
        https://www.h5py.org
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
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

REFERENCES:
    https://www.pgc.umn.edu/guides/arcticdem/data-description/
    https://www.pgc.umn.edu/guides/rema/data-description/
    https://nsidc.org/data/nsidc-0645/versions/1

UPDATE HISTORY:
    Updated 05/2021: print full path of output filename
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 12/2020: H5py deprecation warning change to use make_scale
        using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 10/2020: using argparse to set parameters
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 06/2020: add additional beam check within heights groups
    Updated 10/2019: using delta_time as output HDF5 variable dimensions
        changing Y/N flags to True/False
    Updated 09/2019: fiona for shapefile read.  pyproj for coordinate conversion
        can set the DEM model manually to use the GIMP DEM. verify DEM is finite
        round DEM fill value when creating mask as some DEM tiles are incorrect
        using date functions paralleling public repository. verify output DEM
    Updated 06/2019: assign ArcticDEM by name attribute.  buffer for sub-tiles
    Updated 05/2019: free up memory from invalid tiles. buffer by more geosegs
        include mean vertical residuals and number of ground control points
        assign ArcticDEM polygons by objectid
    Updated 04/2019: buffer DEM tiles using values from neighboring tiles
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
import pyproj
import tarfile
import datetime
import argparse
import osgeo.gdal
import numpy as np
from mpi4py import MPI
import scipy.interpolate
from shapely.geometry import MultiPoint, Polygon
from icesat2_toolkit.convert_delta_time import convert_delta_time
import icesat2_toolkit.time

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
    mask = np.zeros((ysize,xsize),dtype=bool)
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
    xoffset = int((xlimits[0] - xmin)/info_geotiff[1])
    yoffset = int((ymax - ylimits[1])/np.abs(info_geotiff[5]))
    xcount = int((xlimits[1] - xlimits[0])/info_geotiff[1]) + 1
    ycount = int((ylimits[1] - ylimits[0])/np.abs(info_geotiff[5])) + 1
    #-- read data matrix
    im = ds.GetRasterBand(1).ReadAsArray(xoffset, yoffset, xcount, ycount)
    fill_value = ds.GetRasterBand(1).GetNoDataValue()
    fill_value = 0.0 if (fill_value is None) else fill_value
    #-- create mask for finding invalid values
    mask = np.zeros((ycount,xcount),dtype=bool)
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

#-- PURPOSE: read ICESat-2 ATL03 data from NSIDC
#-- interpolate DEM data to x and y coordinates
def main():
    #-- start MPI communicator
    comm = MPI.COMM_WORLD

    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Interpolate DEMs to ICESat-2 ATL03 photon event locations
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL03 file to run')
    #-- working data directory for location of DEM files
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Digital elevation model (REMA, ArcticDEM, GIMP) to run
    #-- set the DEM model to run for a given granule (else set automatically)
    parser.add_argument('--model','-m',
        metavar='DEM', type=str, choices=('REMA', 'ArcticDEM', 'GIMP'),
        help='Digital Elevation Model to run')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- output module information for process
    if args.verbose:
        info(comm.rank,comm.size)
    if args.verbose and (comm.rank==0):
        print('{0} -->'.format(args.file))

    #-- read data from input file
    #-- Open the HDF5 file for reading
    fileID = h5py.File(args.file, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(args.file)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(args.file).pop()

    #-- set the digital elevation model based on ICESat-2 granule
    DEM_MODEL = set_DEM_model(GRN) if (args.model is None) else args.model
    #-- regular expression pattern for extracting parameters from ArcticDEM name
    rx1 = re.compile(r'(\d+)_(\d+)_(\d+)_(\d+)_(\d+m)_(.*?)$', re.VERBOSE)
    #-- full path to DEM directory
    elevation_directory = os.path.join(args.directory,*elevation_dir[DEM_MODEL])
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
    valid_tiles = False

    #-- pyproj transformer for converting from latitude/longitude
    #-- into DEM tile coordinates
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string(tile_epsg)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

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

    #-- copy variables for outputting to HDF5 file
    IS2_atl03_dem = {}
    IS2_atl03_fill = {}
    IS2_atl03_dims = {}
    IS2_atl03_dem_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl03_dem['ancillary_data'] = {}
    IS2_atl03_dem_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl03_dem['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        IS2_atl03_dem_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl03_dem_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl03_beams):
        #-- output data dictionaries for beam
        IS2_atl03_dem[gtx] = dict(heights={},subsetting={})
        IS2_atl03_fill[gtx] = dict(heights={},subsetting={})
        IS2_atl03_dims[gtx] = dict(heights={},subsetting={})
        IS2_atl03_dem_attrs[gtx] = dict(heights={},subsetting={})

        #-- number of photon events
        n_pe, = fileID[gtx]['heights']['h_ph'].shape
        #-- invalid value
        fv = fileID[gtx]['geolocation']['sigma_h'].fillvalue

        #-- define indices to run for specific process
        ind = np.arange(comm.Get_rank(), n_pe, comm.Get_size(), dtype=int)
        #-- extract delta time
        delta_time = fileID[gtx]['heights']['delta_time'][:]
        #-- extract lat/lon
        longitude = fileID[gtx]['heights']['lon_ph'][:]
        latitude = fileID[gtx]['heights']['lat_ph'][:]
        #-- output interpolated digital elevation model
        distributed_dem = np.ma.zeros((n_pe),fill_value=fv,dtype=np.float32)
        distributed_dem.mask = np.ones((n_pe),dtype=bool)
        dem_h = np.ma.zeros((n_pe),fill_value=fv,dtype=np.float32)
        dem_h.mask = np.ones((n_pe),dtype=bool)
        #-- convert projection from latitude/longitude to tile EPSG
        X,Y = transformer.transform(longitude, latitude)

        #-- convert reduced x and y to shapely multipoint object
        xy_point = MultiPoint(list(zip(X[ind], Y[ind])))

        #-- create complete masks for each DEM tile
        associated_map = {}
        for key,poly_obj in tile_dict.items():
            #-- create empty intersection map array for distributing
            distributed_map = np.zeros((n_pe),dtype=int)
            #-- create empty intersection map array for receiving
            associated_map[key] = np.zeros((n_pe),dtype=int)
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
        IS2_atl03_dem_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl03_dem_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl03_dem_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl03_dem_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl03_dem_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl03_dem_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl03_dem_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        #-- group attributes for heights
        IS2_atl03_dem_attrs[gtx]['heights']['Description'] = ("Contains arrays of the "
            "parameters for each received photon.")
        IS2_atl03_dem_attrs[gtx]['heights']['data_rate'] = ("Data are stored at the "
            "photon detection rate.")

        #-- geolocation and time
        #-- delta time
        IS2_atl03_dem[gtx]['heights']['delta_time'] = delta_time
        IS2_atl03_fill[gtx]['heights']['delta_time'] = None
        IS2_atl03_dims[gtx]['heights']['delta_time'] = None
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time'] = {}
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['standard_name'] = "time"
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['calendar'] = "standard"
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl03_dem_attrs[gtx]['heights']['delta_time']['coordinates'] = \
            "lat_ph lon_ph"
        #-- latitude
        IS2_atl03_dem[gtx]['heights']['latitude'] = latitude
        IS2_atl03_fill[gtx]['heights']['latitude'] = None
        IS2_atl03_dims[gtx]['heights']['latitude'] = ['delta_time']
        IS2_atl03_dem_attrs[gtx]['heights']['latitude'] = {}
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['units'] = "degrees_north"
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['long_name'] = "Latitude"
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['standard_name'] = "latitude"
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['description'] = ("Latitude of each "
            "received photon. Computed from the ECF Cartesian coordinates of the bounce point.")
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['valid_min'] = -90.0
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['valid_max'] = 90.0
        IS2_atl03_dem_attrs[gtx]['heights']['latitude']['coordinates'] = \
            "delta_time lon_ph"
        #-- longitude
        IS2_atl03_dem[gtx]['heights']['longitude'] = longitude
        IS2_atl03_fill[gtx]['heights']['longitude'] = None
        IS2_atl03_dims[gtx]['heights']['longitude'] = ['delta_time']
        IS2_atl03_dem_attrs[gtx]['heights']['longitude'] = {}
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['units'] = "degrees_east"
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['long_name'] = "Longitude"
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['standard_name'] = "longitude"
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['description'] = ("Longitude of each "
            "received photon. Computed from the ECF Cartesian coordinates of the bounce point.")
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['valid_min'] = -180.0
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['valid_max'] = 180.0
        IS2_atl03_dem_attrs[gtx]['heights']['longitude']['coordinates'] = \
            "delta_time lat_ph"

        #-- subsetting variables
        IS2_atl03_dem_attrs[gtx]['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce photon events to specific regions of interest.")
        IS2_atl03_dem_attrs[gtx]['subsetting']['data_rate'] = ("Data are stored at the photon "
            "detection rate.")

        #-- for each valid tile
        for key in valid_tiles:
            #-- output mask to HDF5
            IS2_atl03_dem[gtx]['subsetting'][key] = associated_map[key]
            IS2_atl03_fill[gtx]['subsetting'][key] = None
            IS2_atl03_dims[gtx]['subsetting'][key] = ['delta_time']
            IS2_atl03_dem_attrs[gtx]['subsetting'][key] = {}
            IS2_atl03_dem_attrs[gtx]['subsetting'][key]['contentType'] = "referenceInformation"
            IS2_atl03_dem_attrs[gtx]['subsetting'][key]['long_name'] = '{0} Mask'.format(key)
            IS2_atl03_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['description'] = ('Name '
                'of DEM tile {0} encapsulating the land ice segments.').format(tile_attrs[key]['tile'])
            IS2_atl03_dem_attrs[gtx]['subsetting'][key]['source'] = DEM_MODEL
            #-- add DEM attributes
            if DEM_MODEL in ('REMA','ArcticDEM'):
                IS2_atl03_dem_attrs[gtx]['subsetting'][key]['sigma_h'] = tile_attrs[key]['meanresz']
                IS2_atl03_dem_attrs[gtx]['subsetting'][key]['n_gcp'] = tile_attrs[key]['num_gcps']
            IS2_atl03_dem_attrs[gtx]['subsetting'][key]['coordinates'] = \
                "../heights/delta_time ../heights/lat_ph ../heights/lon_ph"

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
            #-- use 200m (10 geosegs and divisible by ArcticDEM and REMA pixels)
            #-- use 750m for GIMP
            bf = 750 if (DEM_MODEL == 'GIMP') else 200
            ny,nx = np.shape(DEM)
            dx = np.abs(xi[1]-xi[0]).astype('i')
            dy = np.abs(yi[1]-yi[0]).astype('i')
            #-- new buffered DEM and mask
            d = np.full((ny+2*bf//dy,nx+2*bf//dx),fv,dtype=np.float32)
            m = np.ones((ny+2*bf//dy,nx+2*bf//dx),dtype=bool)
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
                IMy,IMx=np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                #-- neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles=[IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] #-- LLLCCRRR
                ytiles=[IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] #-- BMTBTBMT
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
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
                        #-- add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)
            elif (DEM_MODEL == 'GIMP'):
                #-- GIMP tiles to read to buffer the image
                IMx,IMy=np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                #-- neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles=[IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] #-- LLLCCRRR
                ytiles=[IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] #-- BMTBTBMT
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
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
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
                    arg = (ytl,xtl,xs,ys,res,vers)
                    bkey = '{0:02d}_{1:02d}_{2}_{3}'.format(*arg)
                    #-- if buffer file is a valid sub-tile within the DEM
                    #-- if file doesn't exist: all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        btar = '{0}.tar.gz'.format(tile_attrs[bkey]['name'])
                        buffer_file = os.path.join(elevation_directory,bsub,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise IOError('{0} not found'.format(buffer_file))
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
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
            distributed_dem.mask[tile_indices] = maskout.astype(bool)
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
        IS2_atl03_dem[gtx]['heights']['dem_h'] = dem_h
        IS2_atl03_fill[gtx]['heights']['dem_h'] = dem_h.fill_value
        IS2_atl03_dims[gtx]['heights']['dem_h'] = ['delta_time']
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h'] = {}
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['units'] = "meters"
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['contentType'] = "referenceInformation"
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['long_name'] = "DEM Height"
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['description'] = ("Height of the DEM, "
            "interpolated by cubic-spline interpolation in the DEM coordinate system to the "
            "PE location.")
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['source'] = DEM_MODEL
        IS2_atl03_dem_attrs[gtx]['heights']['dem_h']['coordinates'] = \
            "delta_time lat_ph lon_ph"

    #-- parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and bool(valid_tiles):
        #-- output HDF5 files with output masks
        fargs = (PRD,DEM_MODEL,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX)
        file_format = '{0}_{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
        #-- print file information
        if args.verbose:
            print('\t{0}'.format(output_file))
        #-- write to output HDF5 file
        HDF5_ATL03_dem_write(IS2_atl03_dem, IS2_atl03_dem_attrs,
            CLOBBER=True, INPUT=os.path.basename(args.file),
            FILL_VALUE=IS2_atl03_fill, DIMENSIONS=IS2_atl03_dims,
            FILENAME=output_file)
        #-- change the permissions mode
        os.chmod(output_file, args.mode)
    #-- close the input file
    fileID.close()

#-- PURPOSE: outputting the interpolated DEM data for ICESat-2 data to HDF5
def HDF5_ATL03_dem_write(IS2_atl03_dem, IS2_atl03_attrs, INPUT=None,
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
    for k,v in IS2_atl03_dem['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl03_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl03_dem.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        h5[gtx] = {}
        #-- add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl03_attrs[gtx][att_name]

        #-- for each output data group
        for key in ['heights','subsetting']:
            #-- create group
            fileID[gtx].create_group(key)
            h5[gtx][key] = {}
            for att_name in ['Description','data_rate']:
                att_val = IS2_atl03_attrs[gtx][key][att_name]
                fileID[gtx][key].attrs[att_name] = att_val

            #-- all variables for group
            groupkeys=set(IS2_atl03_dem[gtx][key].keys())-set(['delta_time'])
            for k in ['delta_time',*sorted(groupkeys)]:
                #-- values and attributes
                v = IS2_atl03_dem[gtx][key][k]
                attrs = IS2_atl03_attrs[gtx][key][k]
                fillvalue = FILL_VALUE[gtx][key][k]
                #-- Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}'.format(gtx,key,k)
                if fillvalue:
                    h5[gtx][key][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype,
                        fillvalue=fillvalue, compression='gzip')
                else:
                    h5[gtx][key][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype,
                        compression='gzip')
                #-- create or attach dimensions for HDF5 variable
                if DIMENSIONS[gtx][key][k]:
                    #-- attach dimensions
                    for i,dim in enumerate(DIMENSIONS[gtx][key][k]):
                        h5[gtx][key][k].dims[i].attach_scale(
                            h5[gtx][key][dim])
                else:
                    #-- make dimension
                    h5[gtx][key][k].make_scale(k)
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 L2A Global Geolocated Photon Data'
    fileID.attrs['summary'] = ("The purpose of ATL03 is to provide along-track "
        "photon data for all 6 ATLAS beams and associated statistics.")
    fileID.attrs['description'] = ("Photon heights determined by ATBD "
        "Algorithm using POD and PPD. All photon events per transmit pulse per "
        "beam. Includes POD and PPD vectors. Classification of each photon by "
        "several ATBD Algorithms.")
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
    #-- add attributes for input ATL03 and ATL09 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl03_dem[gtx]['heights']['longitude']
        lat = IS2_atl03_dem[gtx]['heights']['latitude']
        delta_time = IS2_atl03_dem[gtx]['heights']['delta_time']
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

#-- run main program
if __name__ == '__main__':
    main()
