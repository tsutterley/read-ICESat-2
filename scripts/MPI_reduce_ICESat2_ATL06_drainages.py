#!/usr/bin/env python
u"""
MPI_reduce_ICESat2_ATL06_drainages.py
Written by Tyler Sutterley (12/2022)

Create masks for reducing ICESat-2 data into IMBIE-2 drainage regions

COMMAND LINE OPTIONS:
    -D X, --directory X: Working Data Directory
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
    mpi4py: MPI for Python
        http://pythonhosted.org/mpi4py/
        http://mpi4py.readthedocs.org/en/stable/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
        http://docs.h5py.org/en/stable/mpi.html
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyshp: Python read/write support for ESRI Shapefile format
        https://github.com/GeospatialPython/pyshp
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of altimetry tools
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 11/2021: use delta time as output dimension for HDF5 variables
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 02/2021: use size of array to add to any valid check
        replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 12/2020: H5py deprecation warning change to use make_scale
    Updated 10/2020: using argparse to set parameters.  update pyproj transforms
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: using pyproj for coordinate conversions
    Updated 05/2019: check if beam exists in a try except else clause
    Updated 04/2019: check if subsetted beam contains land ice data
    Updated 02/2019: shapely updates for python3 compatibility
    Written 01/2018
"""
from __future__ import print_function

import sys
import os
import re
import pyproj
import logging
import datetime
import argparse
import warnings
import numpy as np
import icesat2_toolkit as is2tk

# attempt imports
try:
    import h5py
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("h5py not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("mpi4py not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import shapefile
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("shapefile not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    from shapely.geometry import MultiPoint, Polygon
except ModuleNotFoundError:
    warnings.filterwarnings("always")
    warnings.warn("shapely not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# IMBIE-2 Drainage basins
IMBIE_basin_file = {}
IMBIE_basin_file['N'] = ['GRE_Basins_IMBIE2_v1.3','GRE_Basins_IMBIE2_v1.3.shp']
IMBIE_basin_file['S'] = ['ANT_Basins_IMBIE2_v1.6','ANT_Basins_IMBIE2_v1.6.shp']
# basin titles within shapefile to extract
IMBIE_title = {}
IMBIE_title['N']=('CW','NE','NO','NW','SE','SW')
IMBIE_title['S']=('A-Ap','Ap-B','B-C','C-Cp','Cp-D','D-Dp','Dp-E','E-Ep','Ep-F',
    'F-G','G-H','H-Hp','Hp-I','I-Ipp','Ipp-J','J-Jpp','Jpp-K','K-A')

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
        description="""Create masks for reducing ICESat-2 data into IMBIE-2
            drainage regions
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL06 file to run')
    # working data directory for drainage basins shapefiles
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory for mask files')
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

# PURPOSE: set the hemisphere of interest based on the granule
def set_hemisphere(GRANULE):
    if GRANULE in ('10','11','12'):
        projection_flag = 'S'
    elif GRANULE in ('03','04','05'):
        projection_flag = 'N'
    return projection_flag

# PURPOSE: load Greenland or Antarctic drainage basins from IMBIE-2 (Mouginot)
def load_IMBIE2_basins(basin_dir, HEM, EPSG):
    # read drainage basin polylines from shapefile (using splat operator)
    basin_shapefile = os.path.join(basin_dir,*IMBIE_basin_file[HEM])
    shape_input = shapefile.Reader(basin_shapefile)
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    # projections for converting lat/lon to polar stereographic
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = pyproj.CRS.from_epsg(EPSG)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    # python dictionary with shapely polygon objects
    poly_dict = {}
    # for each region
    for REGION in IMBIE_title[HEM]:
        # find record index for region by iterating through shape attributes
        if (HEM == 'S'):
            i,=[i for i,a in enumerate(shape_attributes) if (a[1] == REGION)]
            # extract Polar-Stereographic coordinates for record
            points = np.array(shape_entities[i].points)
            # shapely polygon object for region outline
            poly_obj = Polygon(np.c_[points[:,0],points[:,1]])
        elif (HEM == 'N'):
            # no glaciers or ice caps
            i,=[i for i,a in enumerate(shape_attributes) if (a[0] == REGION)]
            # extract Polar-Stereographic coordinates for record
            points = np.array(shape_entities[i].points)
            # Greenland IMBIE-2 basins can have multiple parts
            parts = shape_entities[i].parts
            parts.append(len(points))
            # list object for x,y coordinates (exterior and holes)
            poly_list = []
            for p1,p2 in zip(parts[:-1],parts[1:]):
                # converting basin lat/lon into Polar-Stereographic
                X,Y = transformer.transform(points[p1:p2,0],points[p1:p2,1])
                poly_list.append(np.c_[X,Y])
            # convert poly_list into Polygon object with holes
            poly_obj = Polygon(poly_list[0],poly_list[1:])
        # check if polygon object is valid
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        # add to total polygon dictionary object
        poly_dict[REGION] = poly_obj
    # return the polygon object and the input file name
    return poly_dict, [basin_shapefile]

# PURPOSE: read ICESat-2 land ice height data (ATL06)
# reduce to IMBIE-2 (Rignot) drainage basins
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

    # Open the HDF5 file for reading
    fileID = h5py.File(args.file, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(args.file)
    # extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(args.file).pop()
    # set the hemisphere flag based on ICESat-2 granule
    HEM = set_hemisphere(GRN)
    # pyproj transformer for converting lat/lon to polar stereographic
    EPSG = dict(N=3413,S=3031)
    crs1 = pyproj.CRS.from_epsg(4326)
    crs2 = pyproj.CRS.from_epsg(EPSG[HEM])
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

    # read each basin and create shapely polygon objects
    # IMBIE-2 basin drainages
    BASIN_TITLE = 'IMBIE-2_BASIN_MASKS'
    DESCRIPTION = 'IMBIE-2_(Rignot_2016)'
    REFERENCE=dict(N='http://imbie.org/imbie-2016/drainage-basins/',
        S='http://imbie.org/imbie-2016/drainage-basins/')

    # read data on rank 0
    if (comm.rank == 0):
        # read each basin and create shapely polygon objects
        # load IMBIE-2 basin drainages
        BASIN,basin_files = load_IMBIE2_basins(args.directory,HEM,EPSG[HEM])
    else:
        # create empty object for list of shapely objects
        BASIN = None

    # Broadcast Shapely basin objects
    BASIN = comm.bcast(BASIN, root=0)
    # combined validity check for all beams
    valid_check = False

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
    IS2_atl06_mask = {}
    IS2_atl06_fill = {}
    IS2_atl06_dims = {}
    IS2_atl06_mask_attrs = {}
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    IS2_atl06_mask['ancillary_data'] = {}
    IS2_atl06_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl06_mask['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        IS2_atl06_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl06_mask_attrs['ancillary_data'][key][att_name] = att_val

    # for each input beam within the file
    for gtx in sorted(IS2_atl06_beams):
        # output data dictionaries for beam
        IS2_atl06_mask[gtx] = dict(land_ice_segments={})
        IS2_atl06_fill[gtx] = dict(land_ice_segments={})
        IS2_atl06_dims[gtx] = dict(land_ice_segments={})
        IS2_atl06_mask_attrs[gtx] = dict(land_ice_segments={})

        # number of segments
        segment_id = fileID[gtx]['land_ice_segments']['segment_id'][:]
        n_seg, = fileID[gtx]['land_ice_segments']['segment_id'].shape
        # invalid value for beam
        fv = fileID[gtx]['land_ice_segments']['h_li'].fillvalue
        # check if there are less segments than processes
        if (n_seg < comm.Get_size()):
            continue

        # define indices to run for specific process
        ind = np.arange(comm.Get_rank(), n_seg, comm.Get_size(), dtype=int)

        # extract delta time
        delta_time = fileID[gtx]['land_ice_segments']['delta_time'][:].copy()
        # extract lat/lon
        longitude = fileID[gtx]['land_ice_segments']['longitude'][:].copy()
        latitude = fileID[gtx]['land_ice_segments']['latitude'][:].copy()
        # convert lat/lon to polar stereographic
        X,Y = transformer.transform(longitude[ind], latitude[ind])
        # convert reduced x and y to shapely multipoint object
        xy_point = MultiPoint(np.c_[X, Y])

        # calculate mask for each drainage basin in the dictionary
        associated_map = {}
        for key,poly_obj in BASIN.items():
            # create distributed intersection map for calculation
            distributed_map = np.zeros((n_seg),dtype=bool)
            # create empty intersection map array for receiving
            associated_map[key] = np.zeros((n_seg),dtype=bool)
            # finds if points are encapsulated (within basin)
            int_test = poly_obj.intersects(xy_point)
            if int_test:
                # extract intersected points
                int_map = list(map(poly_obj.intersects,xy_point))
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

        # group attributes for beam
        IS2_atl06_mask_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl06_mask_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl06_mask_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl06_mask_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl06_mask_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl06_mask_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl06_mask_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        # group attributes for land_ice_segments
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        # geolocation, time and segment ID
        # delta time
        IS2_atl06_mask[gtx]['land_ice_segments']['delta_time'] = delta_time
        IS2_atl06_fill[gtx]['land_ice_segments']['delta_time'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['delta_time'] = None
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['standard_name'] = "time"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['calendar'] = "standard"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['delta_time']['coordinates'] = \
            "segment_id latitude longitude"
        # latitude
        IS2_atl06_mask[gtx]['land_ice_segments']['latitude'] = latitude
        IS2_atl06_fill[gtx]['land_ice_segments']['latitude'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['latitude'] = ['delta_time']
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['units'] = "degrees_north"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['valid_min'] = -90.0
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['valid_max'] = 90.0
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['latitude']['coordinates'] = \
            "segment_id delta_time longitude"
        # longitude
        IS2_atl06_mask[gtx]['land_ice_segments']['longitude'] = longitude
        IS2_atl06_fill[gtx]['land_ice_segments']['longitude'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['longitude'] = ['delta_time']
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['units'] = "degrees_east"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['valid_min'] = -180.0
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['valid_max'] = 180.0
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['longitude']['coordinates'] = \
            "segment_id delta_time latitude"
        # segment ID
        IS2_atl06_mask[gtx]['land_ice_segments']['segment_id'] = segment_id
        IS2_atl06_fill[gtx]['land_ice_segments']['segment_id'] = None
        IS2_atl06_dims[gtx]['land_ice_segments']['segment_id'] = ['delta_time']
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id']['units'] = "1"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id']['contentType'] = "referenceInformation"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id']['long_name'] = "Along-track segment ID number"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id']['description'] = ("A 7 digit number "
            "identifying the along-track geolocation segment number.  These are sequential, starting with "
            "1 for the first segment after an ascending equatorial crossing node. Equal to the segment_id for "
            "the second of the two 20m ATL03 segments included in the 40m ATL06 segment")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['segment_id']['coordinates'] = \
            "delta_time latitude longitude"

        # subsetting variables
        IS2_atl06_mask[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_dims[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce land ice segments to specific regions of interest.")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        # for each valid drainage
        valid_keys = np.array([k for k,v in associated_map.items() if v.any()])
        valid_check |= (np.size(valid_keys) > 0)
        for key in valid_keys:
            # output mask to HDF5
            IS2_atl06_mask[gtx]['land_ice_segments']['subsetting'][key] = associated_map[key]
            IS2_atl06_fill[gtx]['land_ice_segments']['subsetting'][key] = None
            IS2_atl06_dims[gtx]['land_ice_segments']['subsetting'][key] = ['delta_time']
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key] = {}
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key]['contentType'] = \
                "referenceInformation"
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key]['long_name'] = \
                f'{key} Mask'
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key]['description'] = \
                f'Mask calculated using the {key} drainage from {DESCRIPTION}.'
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key]['reference'] = REFERENCE[HEM]
            IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][key]['coordinates'] = \
                "../segment_id ../delta_time ../latitude ../longitude"
        # wait for all processes to finish calculation
        comm.Barrier()

    # parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and valid_check:
        # output HDF5 file with drainage basin masks
        fargs = (PRD,BASIN_TITLE,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX)
        file_format = '{0}_{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
        # print file information
        logging.into(f'\t{output_file}')
        # write to output HDF5 file
        HDF5_ATL06_mask_write(IS2_atl06_mask, IS2_atl06_mask_attrs,
            CLOBBER=True, INPUT=os.path.basename(args.file),
            FILL_VALUE=IS2_atl06_fill, DIMENSIONS=IS2_atl06_dims,
            FILENAME=output_file)
        # change the permissions mode
        os.chmod(output_file, args.mode)
    # close the input file
    fileID.close()

# PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL06_mask_write(IS2_atl06_mask, IS2_atl06_attrs, INPUT=None,
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
    for k,v in IS2_atl06_mask['ancillary_data'].items():
        # Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # write each output beam
    beams = [k for k in IS2_atl06_mask.keys() if bool(re.match(r'gt\d[lr]',k))]
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
        for k in ['segment_id','latitude','longitude','delta_time']:
            # values and attributes
            v = IS2_atl06_mask[gtx]['land_ice_segments'][k]
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

        # add to subsetting variables
        key = 'subsetting'
        fileID[gtx]['land_ice_segments'].create_group(key)
        h5[gtx]['land_ice_segments'][key] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name]
            fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
        for k,v in IS2_atl06_mask[gtx]['land_ice_segments'][key].items():
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
    fileID.attrs['summary'] = ('Subsetting masks for ice-sheets segments '
        'needed to interpret and assess the quality of the height estimates.')
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
        lon = IS2_atl06_mask[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl06_mask[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl06_mask[gtx]['land_ice_segments']['delta_time']
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
    fileID.attrs['software_revision'] = is2tk.utilities.get_git_revision_hash()
    # Closing the HDF5 file
    fileID.close()

# run main program
if __name__ == '__main__':
    main()
