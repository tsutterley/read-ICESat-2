#!/usr/bin/env python
u"""
MPI_DEM_ICESat2_ATL06.py
Written by Tyler Sutterley (12/2022)
Determines which digital elevation model tiles to read for a given ATL06 file
Reads 3x3 array of tiles for points within bounding box of central mosaic tile
Interpolates digital elevation model to locations of ICESat-2 ATL06 segments

ArcticDEM 2m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/
    http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/

REMA 8m digital elevation model tiles
    http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v2.0/
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
    Updated 12/2022: single implicit import of altimetry tools
    Updated 11/2022: new ArcticDEM and REMA mosaic index shapefiles
        verify coordinate reference system attribute from shapefile
    Updated 09/2022: use tar virtual file system to extract images
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 11/2021: use delta time as output dimension for HDF5 variables
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 12/2020: H5py deprecation warning change to use make_scale
        using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 10/2020: using argparse to set parameters
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
import pyproj
import logging
import tarfile
import datetime
import argparse
import warnings
import numpy as np
import scipy.interpolate
import icesat2_toolkit as is2tk

# attempt imports
try:
    import fiona
except ModuleNotFoundError:
    warnings.warn("fiona not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import h5py
except ModuleNotFoundError:
    warnings.warn("h5py not available", ImportWarning)
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    warnings.warn("mpi4py not available", ImportWarning)
try:
    import osgeo.gdal
except ModuleNotFoundError:
    warnings.warn("GDAL not available", ImportWarning)
try:
    from shapely.geometry import MultiPoint, Polygon
except ModuleNotFoundError:
    warnings.warn("shapely not available", ImportWarning)

# digital elevation models
elevation_dir = {}
elevation_tile_index = {}
# ArcticDEM
elevation_dir['ArcticDEM'] = ['ArcticDEM']
elevation_tile_index['ArcticDEM'] = 'ArcticDEM_Mosaic_Index_v3_shp.zip'
# GIMP DEM
elevation_dir['GIMP'] = ['GIMP','30m']
elevation_tile_index['GIMP'] = 'gimpdem_Tile_Index_Rel1.1.zip'
# REMA DEM
elevation_dir['REMA'] = ['REMA']
elevation_tile_index['REMA'] = 'REMA_Mosaic_Index_v2_shp.zip'

# PURPOSE: keep track of MPI threads
def info(rank, size):
    logging.info(f'Rank {rank+1:d} of {size:d}')
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Interpolate DEMs to ICESat-2 ATL06 land ice locations
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL06 file to run')
    # working data directory for location of DEM files
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    # Digital elevation model (REMA, ArcticDEM, GIMP) to run
    # set the DEM model to run for a given granule (else set automatically)
    parser.add_argument('--model','-m',
        metavar='DEM', type=str, choices=('REMA', 'ArcticDEM', 'GIMP'),
        help='Digital Elevation Model to run')
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

# PURPOSE: set the DEM model to interpolate based on the input granule
def set_DEM_model(GRANULE):
    if GRANULE in ('10','11','12'):
        DEM_MODEL = 'REMA'
    elif GRANULE in ('02','03','04','05','06'):
        DEM_MODEL = 'ArcticDEM'
    return DEM_MODEL

# PURPOSE: read zip file containing index shapefiles for finding DEM tiles
def read_DEM_index(index_file, DEM_MODEL):
    # read the compressed shapefile and extract entities
    shape = fiona.open(f'zip://{os.path.expanduser(index_file)}')
    # extract coordinate reference system
    if ('init' in shape.crs.keys()):
        epsg = pyproj.CRS(shape.crs['init']).to_epsg()
    else:
        epsg = pyproj.CRS(shape.crs).to_epsg()
    # extract attribute indice for DEM tile (REMA,GIMP) or name (ArcticDEM)
    if (DEM_MODEL == 'REMA'):
        # REMA index file attributes:
        # name: DEM mosaic name for tile (file name without suffix)
        # tile: DEM tile identifier (IMy_IMx)
        # nd_value: fill value for elements with no data
        # resolution: DEM horizontal spatial resolution (meters)
        # creationda: creation date
        # raster: (empty)
        # fileurl: link to file on PGC server
        # spec_type: specific type (DEM)
        # qual: density of scenes within tile (0 to 1)
        # reg_src: DEM registration source (ICESat or neighbor align)
        # num_gcps: number of ground control points
        # meanresz: mean vertical residual (meters)
        # active: (1)
        # qc: (2)
        # rel_ver: release version
        # num_comp: number of components
        # st_area_sh: tile area (meters^2)
        # st_length_: perimeter length of tile (meters)
        field = 'tile'
    elif (DEM_MODEL == 'GIMP'):
        # GIMP index file attributes (from make_GIMP_tile_shapefile.py):
        # name: DEM mosaic name for tile (file name without suffix)
        # tile: DEM tile identifier (IMy_IMx)
        # nd_value: fill value for elements with no data
        # resolution: DEM horizontal spatial resolution (meters)
        # fileurl: link to file on NSIDC server
        # spec_type: specific type (DEM)
        # reg_src: DEM registration source (ICESat or neighbor align)
        # rel_ver: release version
        # num_comp: number of components
        # st_area_sh: tile area (meters^2)
        # st_length_: perimeter length of tile (meters)
        field = 'tile'
    elif (DEM_MODEL == 'ArcticDEM'):
        # ArcticDEM index file attributes:
        # objectid: DEM tile object identifier for sub-tile
        # name: DEM mosaic name for sub-tile (file name without suffix)
        # tile: DEM tile identifier (IMy_IMx) (non-unique for sub-tiles)
        # nd_value: fill value for elements with no data
        # resolution: DEM horizontal spatial resolution (meters)
        # creationda: creation date
        # raster: (empty)
        # fileurl: link to file on PGC server
        # spec_type: specific type (DEM)
        # qual: density of scenes within tile (0 to 1)
        # reg_src: DEM registration source (ICESat or neighbor align)
        # num_gcps: number of ground control points
        # meanresz: mean vertical residual (meters)
        # active: (1)
        # qc: (2)
        # rel_ver: release version
        # num_comp: number of components
        # st_area_sh: tile area (meters^2)
        # st_length_: perimeter length of tile (meters)
        field = 'name'
    # create python dictionary for each polygon object
    poly_dict = {}
    attrs_dict = {}
    # extract the entities and assign by tile name
    for i,ent in enumerate(shape.values()):
        # tile or name attributes
        if DEM_MODEL in ('REMA','GIMP'):
            tile = str(ent['properties'][field])
        else:
            tile, = re.findall(r'^(\d+_\d+_\d+_\d+)',ent['properties'][field])
        # extract attributes and assign by tile
        attrs_dict[tile] = {}
        for key,val in ent['properties'].items():
            attrs_dict[tile][key] = val
        # upper-left, upper-right, lower-right, lower-left, upper-left
        ul,ur,lr,ll,ul2 = ent['geometry']['coordinates'].pop()
        # tile boundaries
        attrs_dict[tile]['xmin'] = ul[0]
        attrs_dict[tile]['xmax'] = lr[0]
        attrs_dict[tile]['ymin'] = lr[1]
        attrs_dict[tile]['ymax'] = ul[1]
        # extract Polar Stereographic coordinates for entity
        x = [ul[0],ur[0],lr[0],ll[0],ul2[0]]
        y = [ul[1],ur[1],lr[1],ll[1],ul2[1]]
        poly_obj = Polygon(list(zip(x,y)))
        # Valid Polygon may not possess overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        poly_dict[tile] = poly_obj
    # close the file
    shape.close()
    # return the dictionaries of polygon objects and attributes
    return (poly_dict,attrs_dict,epsg)

# PURPOSE: read DEM tile file from gzipped tar files
def read_DEM_file(elevation_file, nd_value):
    # open file with tarfile (read)
    tar = tarfile.open(name=elevation_file, mode='r:gz')
    # find dem geotiff file within tar file
    member, = [m for m in tar.getmembers() if re.search(r'dem\.tif',m.name)]
    # use GDAL virtual file systems to read dem
    mmap_name = f"/vsitar/{elevation_file}/{member.name}"
    ds = osgeo.gdal.Open(mmap_name)
    # read data matrix
    im = ds.GetRasterBand(1).ReadAsArray()
    fill_value = ds.GetRasterBand(1).GetNoDataValue()
    fill_value = 0.0 if (fill_value is None) else fill_value
    # get dimensions
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    # create mask for finding invalid values
    mask = np.zeros((ysize,xsize),dtype=bool)
    indy,indx = np.nonzero((im == fill_value) | (~np.isfinite(im)) |
        (np.ceil(im) == np.ceil(fill_value)))
    mask[indy,indx] = True
    # verify that values are finite by replacing with nd_value
    im[indy,indx] = nd_value
    # get geotiff info
    info_geotiff = ds.GetGeoTransform()
    # calculate image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    xmax = xmin + (xsize-1)*info_geotiff[1]
    ymin = ymax + (ysize-1)*info_geotiff[5]
    # close files
    ds = None
    osgeo.gdal.Unlink(mmap_name)
    tar.close()
    # create image x and y arrays
    xi = np.arange(xmin,xmax+info_geotiff[1],info_geotiff[1])
    yi = np.arange(ymax,ymin+info_geotiff[5],info_geotiff[5])
    # return values (flip y values to be monotonically increasing)
    return (im[::-1,:],mask[::-1,:],xi,yi[::-1])

# PURPOSE: read DEM tile file from gzipped tar files to buffer main tile
def read_DEM_buffer(elevation_file, xlimits, ylimits, nd_value):
    # open file with tarfile (read)
    tar = tarfile.open(name=elevation_file, mode='r:gz')
    # find dem geotiff file within tar file
    member, = [m for m in tar.getmembers() if re.search(r'dem\.tif',m.name)]
    # use GDAL virtual file systems to read dem
    mmap_name = f"/vsitar/{elevation_file}/{member.name}"
    ds = osgeo.gdal.Open(mmap_name)
    # get geotiff info
    info_geotiff = ds.GetGeoTransform()
    # original image extents
    xmin = info_geotiff[0]
    ymax = info_geotiff[3]
    # reduce input image with GDAL
    # Specify offset and rows and columns to read
    xoffset = int((xlimits[0] - xmin)/info_geotiff[1])
    yoffset = int((ymax - ylimits[1])/np.abs(info_geotiff[5]))
    xcount = int((xlimits[1] - xlimits[0])/info_geotiff[1]) + 1
    ycount = int((ylimits[1] - ylimits[0])/np.abs(info_geotiff[5])) + 1
    # read data matrix
    im = ds.GetRasterBand(1).ReadAsArray(xoffset, yoffset, xcount, ycount)
    fill_value = ds.GetRasterBand(1).GetNoDataValue()
    fill_value = 0.0 if (fill_value is None) else fill_value
    # create mask for finding invalid values
    mask = np.zeros((ycount,xcount),dtype=bool)
    indy,indx = np.nonzero((im == fill_value) | (~np.isfinite(im)) |
        (np.ceil(im) == np.ceil(fill_value)))
    mask[indy,indx] = True
    # verify that values are finite by replacing with nd_value
    im[indy,indx] = nd_value
    # reduced x and y limits of image
    xmin_reduced = xmin + xoffset*info_geotiff[1]
    xmax_reduced = xmin + xoffset*info_geotiff[1] + (xcount-1)*info_geotiff[1]
    ymax_reduced = ymax + yoffset*info_geotiff[5]
    ymin_reduced = ymax + yoffset*info_geotiff[5] + (ycount-1)*info_geotiff[5]
    # close files
    ds = None
    osgeo.gdal.Unlink(mmap_name)
    tar.close()
    # create image x and y arrays
    xi = np.arange(xmin_reduced,xmax_reduced+info_geotiff[1],info_geotiff[1])
    yi = np.arange(ymax_reduced,ymin_reduced+info_geotiff[5],info_geotiff[5])
    # return values (flip y values to be monotonically increasing)
    return (im[::-1,:],mask[::-1,:],xi,yi[::-1])

# PURPOSE: read ICESat-2 land ice height data (ATL06)
# interpolate DEM data to x and y coordinates
def main():
    # start MPI communicator
    comm = MPI.COMM_WORLD

    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevel = logging.INFO if args.verbose else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # output module information for process
    info(comm.rank,comm.size)
    if (comm.rank == 0):
        logging.info(f'{args.file} -->')

    # read data from input file
    # Open the HDF5 file for reading
    fileID = h5py.File(args.file, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(args.file)
    # extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(args.file).pop()

    # set the digital elevation model based on ICESat-2 granule
    DEM_MODEL = set_DEM_model(GRN) if (args.model is None) else args.model
    # regular expression pattern for extracting parameters from ArcticDEM name
    rx1 = re.compile(r'(\d+)_(\d+)_(\d+)_(\d+)_(\d+m)_(.*?)$', re.VERBOSE)
    # full path to DEM directory
    elevation_directory = os.path.join(args.directory,*elevation_dir[DEM_MODEL])
    # zip file containing index shapefiles for finding DEM tiles
    index_file=os.path.join(elevation_directory,elevation_tile_index[DEM_MODEL])

    # read data on rank 0
    if (comm.rank == 0):
        # read index file for determining which tiles to read
        tile_dict,tile_attrs,tile_epsg = read_DEM_index(index_file,DEM_MODEL)
    else:
        # create empty object for list of shapely objects
        tile_dict = None
        tile_attrs = None
        tile_epsg = None

    # Broadcast Shapely polygon objects
    tile_dict = comm.bcast(tile_dict, root=0)
    tile_attrs = comm.bcast(tile_attrs, root=0)
    tile_epsg = comm.bcast(tile_epsg, root=0)
    valid_tiles = False

    # pyproj transformer for converting from latitude/longitude
    # into DEM tile coordinates
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = pyproj.CRS.from_epsg(tile_epsg)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    # read each input beam within the file
    IS2_atl06_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains land ice data
        try:
            fileID[gtx]['land_ice_segments']['segment_id']
        except KeyError:
            pass
        else:
            IS2_atl06_beams.append(gtx)

    # copy variables for outputting to HDF5 file
    IS2_atl06_dem = {}
    IS2_atl06_fill = {}
    IS2_atl06_dims = {}
    IS2_atl06_dem_attrs = {}
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    IS2_atl06_dem['ancillary_data'] = {}
    IS2_atl06_dem_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl06_dem['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        IS2_atl06_dem_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl06_dem_attrs['ancillary_data'][key][att_name] = att_val

    # for each input beam within the file
    for gtx in sorted(IS2_atl06_beams):
        # output data dictionaries for beam
        IS2_atl06_dem[gtx] = dict(land_ice_segments={})
        IS2_atl06_fill[gtx] = dict(land_ice_segments={})
        IS2_atl06_dims[gtx] = dict(land_ice_segments={})
        IS2_atl06_dem_attrs[gtx] = dict(land_ice_segments={})

        # number of segments
        segment_id = fileID[gtx]['land_ice_segments']['segment_id'][:]
        n_seg, = fileID[gtx]['land_ice_segments']['segment_id'].shape
        # invalid value
        fv = fileID[gtx]['land_ice_segments']['h_li'].fillvalue

        # define indices to run for specific process
        ind = np.arange(comm.Get_rank(), n_seg, comm.Get_size(), dtype=int)

        # extract delta time
        delta_time = fileID[gtx]['land_ice_segments']['delta_time'][:].copy()
        # extract lat/lon
        longitude = fileID[gtx]['land_ice_segments']['longitude'][:].copy()
        latitude = fileID[gtx]['land_ice_segments']['latitude'][:].copy()
        # output interpolated digital elevation model
        distributed_dem = np.ma.zeros((n_seg),fill_value=fv,dtype=np.float32)
        distributed_dem.mask = np.ones((n_seg),dtype=bool)
        dem_h = np.ma.zeros((n_seg),fill_value=fv,dtype=np.float32)
        dem_h.mask = np.ones((n_seg),dtype=bool)
        # convert projection from latitude/longitude to tile EPSG
        X,Y = transformer.transform(longitude, latitude)

        # convert reduced x and y to shapely multipoint object
        xy_point = MultiPoint(list(zip(X[ind], Y[ind])))

        # create complete masks for each DEM tile
        associated_map = {}
        for key,poly_obj in tile_dict.items():
            # create empty intersection map array for distributing
            distributed_map = np.zeros((n_seg),dtype=int)
            # create empty intersection map array for receiving
            associated_map[key] = np.zeros((n_seg),dtype=int)
            # finds if points are encapsulated (within tile)
            int_test = poly_obj.intersects(xy_point)
            if int_test:
                # extract intersected points
                int_map = list(map(poly_obj.intersects, xy_point.geoms))
                int_indices, = np.nonzero(int_map)
                # set distributed_map indices to True for intersected points
                distributed_map[ind[int_indices]] = True
            # communicate output MPI matrices between ranks
            # operation is a logical "or" across the elements.
            comm.Allreduce(sendbuf=[distributed_map, MPI.BOOL], \
                recvbuf=[associated_map[key], MPI.BOOL], op=MPI.LOR)
            distributed_map = None
        # wait for all processes to finish calculation
        comm.Barrier()
        # find valid tiles and free up memory from invalid tiles
        valid_tiles = [k for k,v in associated_map.items() if v.any()]
        invalid_tiles = sorted(set(associated_map.keys()) - set(valid_tiles))
        for key in invalid_tiles:
            associated_map[key] = None

        # group attributes for beam
        IS2_atl06_dem_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl06_dem_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl06_dem_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl06_dem_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl06_dem_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl06_dem_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl06_dem_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        # group attributes for land_ice_segments
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        # geolocation, time and segment ID
        # delta time
        IS2_atl06_dem[gtx]['land_ice_segments']['delta_time'] = delta_time
        IS2_atl06_fill[gtx]['land_ice_segments']['delta_time'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['delta_time'] = None
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
        # latitude
        IS2_atl06_dem[gtx]['land_ice_segments']['latitude'] = latitude
        IS2_atl06_fill[gtx]['land_ice_segments']['latitude'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['latitude'] = ['delta_time']
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
        # longitude
        IS2_atl06_dem[gtx]['land_ice_segments']['longitude'] = longitude
        IS2_atl06_fill[gtx]['land_ice_segments']['longitude'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['longitude'] = ['delta_time']
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
        # segment ID
        IS2_atl06_dem[gtx]['land_ice_segments']['segment_id'] = segment_id
        IS2_atl06_fill[gtx]['land_ice_segments']['segment_id'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['segment_id'] = ['delta_time']
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

        # dem variables
        IS2_atl06_dem[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_dims[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['Description'] = ("The dem group "
            "contains the reference digital elevation model and geoid heights.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        # subsetting variables
        IS2_atl06_dem[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dims[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce land ice segments to specific regions of interest.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        # for each valid tile
        for key in valid_tiles:
            # output mask to HDF5
            IS2_atl06_dem[gtx]['land_ice_segments']['subsetting'][key] = associated_map[key]
            IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'][key] = None
            IS2_atl06_dims[gtx]['land_ice_segments']['subsetting'][key] = ['delta_time']
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key] = {}
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['contentType'] = "referenceInformation"
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['long_name'] = f'{key} Mask'
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['description'] = ('Name '
                'of DEM tile {0} encapsulating the land ice segments.').format(tile_attrs[key]['tile'])
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['source'] = DEM_MODEL
            # add DEM attributes
            if DEM_MODEL in ('REMA','ArcticDEM'):
                IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['sigma_h'] = tile_attrs[key]['meanresz']
                IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['n_gcp'] = tile_attrs[key]['num_gcps']
            IS2_atl06_dem_attrs[gtx]['land_ice_segments']['subsetting'][key]['coordinates'] = \
                "../segment_id ../delta_time ../latitude ../longitude"

        # read and interpolate DEM to coordinates in parallel
        for t in range(comm.Get_rank(), len(valid_tiles), comm.Get_size()):
            key = valid_tiles[t]
            sub = tile_attrs[key]['tile']
            name = tile_attrs[key]['name']
            # read central DEM file (geotiff within gzipped tar file)
            tar = f'{name}.tar.gz'
            elevation_file = os.path.join(elevation_directory,sub,tar)
            DEM,MASK,xi,yi = read_DEM_file(elevation_file,fv)
            # buffer DEM using values from adjacent tiles
            # use 400m (10 geosegs and divisible by ArcticDEM and REMA pixels)
            # use 1500m for GIMP
            bf = 1500 if (DEM_MODEL == 'GIMP') else 400
            ny,nx = np.shape(DEM)
            dx = np.abs(xi[1]-xi[0]).astype('i')
            dy = np.abs(yi[1]-yi[0]).astype('i')
            # new buffered DEM and mask
            d = np.full((ny+2*bf//dy,nx+2*bf//dx),fv,dtype=np.float32)
            m = np.ones((ny+2*bf//dy,nx+2*bf//dx),dtype=bool)
            d[bf//dy:-bf//dy,bf//dx:-bf//dx] = DEM.copy()
            m[bf//dy:-bf//dy,bf//dx:-bf//dx] = MASK.copy()
            DEM,MASK = (None,None)
            # new buffered image x and y coordinates
            x = (xi[0] - bf) + np.arange((nx+2*bf//dx))*dx
            y = (yi[0] - bf) + np.arange((ny+2*bf//dy))*dy
            # min and max of left column, center column, right column
            XL,XC,XR = [[xi[0]-bf,xi[0]-dx],[xi[0],xi[-1]],[xi[-1]+dx,xi[-1]+bf]]
            xlimits = [XL,XL,XL,XC,XC,XR,XR,XR] # LLLCCRRR
            # min and max of bottom row, middle row, top row
            YB,YM,YT = [[yi[0]-bf,yi[0]-dy],[yi[0],yi[-1]],[yi[-1]+dy,yi[-1]+bf]]
            ylimits = [YB,YM,YT,YB,YT,YB,YM,YT] # BMTBTBMT

            # buffer using neighbor tiles (REMA/GIMP) or sub-tiles (ArcticDEM)
            if (DEM_MODEL == 'REMA'):
                # REMA tiles to read to buffer the image
                IMy,IMx=np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                # neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles=[IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] # LLLCCRRR
                ytiles=[IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] # BMTBTBMT
                for xtl,ytl,xlim,ylim in zip(xtiles,ytiles,xlimits,ylimits):
                    # read DEM file (geotiff within gzipped tar file)
                    bkey = f'{ytl:02d}_{xtl:02d}'
                    # if buffer file is a valid tile within the DEM
                    # if file doesn't exist: will be all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        bname = tile_attrs[bkey]['name']
                        btar = f'{bname}.tar.gz'
                        buffer_file = os.path.join(elevation_directory,bkey,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise FileNotFoundError(f'{buffer_file} not found')
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
                        # add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)
            elif (DEM_MODEL == 'GIMP'):
                # GIMP tiles to read to buffer the image
                IMx,IMy=np.array(re.findall(r'(\d+)_(\d+)',sub).pop(),dtype='i')
                # neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                xtiles=[IMx-1,IMx-1,IMx-1,IMx,IMx,IMx+1,IMx+1,IMx+1] # LLLCCRRR
                ytiles=[IMy-1,IMy,IMy+1,IMy-1,IMy+1,IMy-1,IMy,IMy+1] # BMTBTBMT
                for xtl,ytl,xlim,ylim in zip(xtiles,ytiles,xlimits,ylimits):
                    # read DEM file (geotiff within gzipped tar file)
                    bkey = f'{xtl:d}_{ytl:d}'
                    # if buffer file is a valid tile within the DEM
                    # if file doesn't exist: will be all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        bname = tile_attrs[bkey]['name']
                        btar = f'{bname}.tar.gz'
                        buffer_file = os.path.join(elevation_directory,bkey,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise FileNotFoundError(f'{buffer_file} not found')
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
                        # add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)
            elif (DEM_MODEL == 'ArcticDEM'):
                # ArcticDEM sub-tiles to read to buffer the image
                # extract parameters from tile filename
                IMy,IMx,STx,STy,res,vers = rx1.findall(name).pop()
                IMy,IMx,STx,STy = np.array([IMy,IMx,STx,STy],dtype='i')
                # neighboring tiles for buffering DEM (LB,LM,LT,CB,CT,RB,RM,RT)
                # LLLCCRRR
                xtiles = [IMx+(STx-2)//2,IMx+(STx-2)//2,IMx+(STx-2)//2,IMx,IMx,
                    IMx+STx//2,IMx+STx//2,IMx+STx//2]
                xsubtiles = [(STx-2) % 2 + 1,(STx-2) % 2 + 1,(STx-2) % 2 + 1,
                    STx,STx,STx % 2 + 1,STx % 2 + 1,STx % 2 + 1]
                # BMTBTBMT
                ytiles = [IMy+(STy-2)//2,IMy,IMy+STy//2,IMy+(STy-2)//2,
                    IMy+STy//2,IMy+(STy-2)//2,IMy,IMy+STy//2]
                ysubtiles = [(STy-2) % 2 + 1,STy,STy % 2 + 1,(STy-2) % 2 + 1,
                    STy % 2 + 1,(STy-2) % 2 + 1,STy,STy % 2 + 1]
                # for each buffer tile and sub-tile
                kwargs = (xtiles,ytiles,xsubtiles,ysubtiles,xlimits,ylimits)
                for xtl,ytl,xs,ys,xlim,ylim in zip(*kwargs):
                    # read DEM file (geotiff within gzipped tar file)
                    bkey = f'{ytl:02d}_{xtl:02d}_{xs}_{ys}'
                    # if buffer file is a valid sub-tile within the DEM
                    # if file doesn't exist: all fill value with all mask
                    if bkey in tile_attrs.keys():
                        bsub = tile_attrs[bkey]['tile']
                        bname = tile_attrs[bkey]['name']
                        btar = f'{bname}.tar.gz'
                        buffer_file = os.path.join(elevation_directory,bsub,btar)
                        if not os.access(buffer_file, os.F_OK):
                            raise FileNotFoundError(f'{buffer_file} not found')
                        DEM,MASK,x1,y1=read_DEM_buffer(buffer_file,xlim,ylim,fv)
                        xmin = int((x1[0] - x[0])//dx)
                        xmax = int((x1[-1] - x[0])//dx) + 1
                        ymin = int((y1[0] - y[0])//dy)
                        ymax = int((y1[-1] - y[0])//dy) + 1
                        # add to buffered DEM and mask
                        d[ymin:ymax,xmin:xmax] = DEM.copy()
                        m[ymin:ymax,xmin:xmax] = MASK.copy()
                        DEM,MASK = (None,None)

            # indices of x and y coordinates within tile
            tile_indices, = np.nonzero(associated_map[key])
            # use spline interpolation to calculate DEM values at coordinates
            f1 = scipy.interpolate.RectBivariateSpline(x,y,d.T,kx=1,ky=1)
            f2 = scipy.interpolate.RectBivariateSpline(x,y,m.T,kx=1,ky=1)
            dataout = f1.ev(X[tile_indices],Y[tile_indices])
            maskout = f2.ev(X[tile_indices],Y[tile_indices])
            # save DEM to output variables
            distributed_dem.data[tile_indices] = dataout
            distributed_dem.mask[tile_indices] = maskout.astype(bool)
            # clear DEM and mask variables
            f1,f2,dataout,maskout,d,m = (None,None,None,None,None,None)

        # communicate output MPI matrices between ranks
        # operations are element summations and logical "and" across elements
        comm.Allreduce(sendbuf=[distributed_dem.data, MPI.FLOAT], \
            recvbuf=[dem_h.data, MPI.FLOAT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[distributed_dem.mask, MPI.BOOL], \
            recvbuf=[dem_h.mask, MPI.BOOL], op=MPI.LAND)
        distributed_dem = None
        # wait for all distributed processes to finish for beam
        comm.Barrier()

        # output interpolated DEM to HDF5
        dem_h.mask[np.abs(dem_h.data) >= 1e4] = True
        dem_h.data[dem_h.mask] = dem_h.fill_value
        IS2_atl06_dem[gtx]['land_ice_segments']['dem']['dem_h'] = dem_h
        IS2_atl06_fill[gtx]['land_ice_segments']['dem']['dem_h'] = dem_h.fill_value
        IS2_atl06_dims[gtx]['land_ice_segments']['dem']['dem_h'] = ['delta_time']
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h'] = {}
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['units'] = "meters"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['contentType'] = "referenceInformation"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['long_name'] = "DEM Height"
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['description'] = ("Height of "
            "the DEM, interpolated by bivariate-spline interpolation in the DEM coordinate system "
            "to the segment location.")
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['source'] = DEM_MODEL
        IS2_atl06_dem_attrs[gtx]['land_ice_segments']['dem']['dem_h']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

    # parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and bool(valid_tiles):
        # output HDF5 files with output masks
        fargs = (PRD,DEM_MODEL,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX)
        file_format = '{0}_{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
        # print file information
        logging.info(f'\t{output_file}')
        # write to output HDF5 file
        HDF5_ATL06_dem_write(IS2_atl06_dem, IS2_atl06_dem_attrs,
            CLOBBER=True, INPUT=os.path.basename(args.file),
            FILL_VALUE=IS2_atl06_fill, DIMENSIONS=IS2_atl06_dims,
            FILENAME=output_file)
        # change the permissions mode
        os.chmod(output_file, args.mode)
    # close the input file
    fileID.close()

# PURPOSE: outputting the interpolated DEM data for ICESat-2 data to HDF5
def HDF5_ATL06_dem_write(IS2_atl06_dem, IS2_atl06_attrs, INPUT=None,
    FILENAME='', FILL_VALUE=None, DIMENSIONS=None, CLOBBER=True):
    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)

    # create HDF5 records
    h5 = {}

    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k,v in IS2_atl06_dem['ancillary_data'].items():
        # Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # write each output beam
    beams = [k for k in IS2_atl06_dem.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        # add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl06_attrs[gtx][att_name]
        # create land_ice_segments group
        fileID[gtx].create_group('land_ice_segments')
        h5[gtx] = dict(land_ice_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl06_attrs[gtx]['land_ice_segments'][att_name]
            fileID[gtx]['land_ice_segments'].attrs[att_name] = att_val

        # segment_id, geolocation, time and height variables
        for k in ['delta_time','latitude','longitude','segment_id']:
            # values and attributes
            v = IS2_atl06_dem[gtx]['land_ice_segments'][k]
            attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][k]
            # Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments',k)
            if fillvalue:
                h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                    compression='gzip')
            else:
                h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                    np.shape(v), data=v, dtype=v.dtype, compression='gzip')
            # create or attach dimensions for HDF5 variable
            if DIMENSIONS[gtx]['land_ice_segments'][k]:
                # attach dimensions
                for i,dim in enumerate(DIMENSIONS[gtx]['land_ice_segments'][k]):
                    h5[gtx]['land_ice_segments'][k].dims[i].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
            else:
                # make dimension
                h5[gtx]['land_ice_segments'][k].make_scale(k)
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][k].attrs[att_name] = att_val

        # add to output variables
        for key in ['subsetting','dem']:
            fileID[gtx]['land_ice_segments'].create_group(key)
            h5[gtx]['land_ice_segments'][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name]
                fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
            for k,v in IS2_atl06_dem[gtx]['land_ice_segments'][key].items():
                # attributes
                attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][key][k]
                fillvalue = FILL_VALUE[gtx]['land_ice_segments'][key][k]
                # Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}/{3}'.format(gtx,'land_ice_segments',key,k)
                if fillvalue:
                    h5[gtx]['land_ice_segments'][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
                else:
                    h5[gtx]['land_ice_segments'][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                        dtype=v.dtype, compression='gzip')
                # attach dimensions
                for i,dim in enumerate(DIMENSIONS[gtx]['land_ice_segments'][key][k]):
                    h5[gtx]['land_ice_segments'][key][k].dims[i].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
                # add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx]['land_ice_segments'][key][k].attrs[att_name] = att_val

    # HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Land Ice Height'
    fileID.attrs['summary'] = ('Subsetting masks and geophysical parameters '
        'for land ice segments needed to interpret and assess the quality '
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
    # add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    # add attributes for input ATL06 file
    fileID.attrs['input_files'] = os.path.basename(INPUT)
    # find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl06_dem[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl06_dem[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl06_dem[gtx]['land_ice_segments']['delta_time']
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
    # convert start and end time from ATLAS SDP seconds into UTC time
    time_utc = is2tk.convert_delta_time(np.array([tmn,tmx]))
    # convert to calendar date
    YY,MM,DD,HH,MN,SS = is2tk.time.convert_julian(time_utc['julian'],
        format='tuple')
    # add attributes with measurement date start, end and duration
    tcs = datetime.datetime(int(YY[0]), int(MM[0]), int(DD[0]),
        int(HH[0]), int(MN[0]), int(SS[0]), int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(int(YY[1]), int(MM[1]), int(DD[1]),
        int(HH[1]), int(MN[1]), int(SS[1]), int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = f'{tmx-tmn:0.0f}'
    # add software information
    fileID.attrs['software_reference'] = is2tk.version.project_name
    fileID.attrs['software_version'] = is2tk.version.full_version
    # Closing the HDF5 file
    fileID.close()

# run main program
if __name__ == '__main__':
    main()
