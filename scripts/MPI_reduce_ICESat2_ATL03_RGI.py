#!/usr/bin/env python
u"""
MPI_reduce_ICESat2_ATL03_RGI.py
Written by Tyler Sutterley (12/2022)

Create masks for reducing ICESat-2 data to the Randolph Glacier Inventory
    https://www.glims.org/RGI/rgi60_dl.html

COMMAND LINE OPTIONS:
    -D X, --directory X: Working Data Directory
    -R X, --region X: region of Randolph Glacier Inventory to run
        1: Alaska
        2: Western Canada and USA
        3: Arctic Canada North
        4: Arctic Canada South
        5: Greenland Periphery
        6: Iceland
        7: Svalbard
        8: Scandinavia
        9: Russian Arctic
        10: North Asia
        11: Central Europe
        12: Caucasus, Middle East
        13: Central Asia
        14: South Asia West
        15: South Asia East
        16: Low Latitudes
        17: Southern Andes
        18: New Zealand
        19: Antarctic, Subantarctic
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

PROGRAM DEPENDENCIES:
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 12/2022: single implicit import of altimetry tools
    Updated 07/2022: place some imports behind try/except statements
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 12/2020: H5py deprecation warning change to use make_scale
    Updated 10/2020: using argparse to set parameters
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 06/2020: add additional beam check within heights groups
    Updated 10/2019: using delta_time as output HDF5 variable dimensions
        changing Y/N flags to True/False
    Updated 09/2019: using date functions paralleling public repository
    Updated 05/2019: read shapefiles from RGI provided zip files
        check if beam exists in a try except else clause
    Updated 04/2019: check if subsetted beam contains land ice data
        save RGIId for valid points to determine containing RGI polygons
    Written 04/2019
"""
from __future__ import print_function

import sys
import os
import re
import io
import logging
import zipfile
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
        description="""Create masks for reducing ICESat-2 photon event data to
            the Randolph Glacier Inventory (RGI)
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL03 file to run')
    # working data directory for location of RGI files
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory for mask files')
    # region of Randolph Glacier Inventory to run
    parser.add_argument('--region','-r',
        metavar='RGI', type=int, choices=range(1,20),
        help='region of Randolph Glacier Inventory to run')
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

# PURPOSE: load zip file containing Randolph Glacier Inventory shapefiles
def load_glacier_inventory(RGI_DIRECTORY,RGI_REGION):
    # list of Randolph Glacier Inventory files
    RGI_files = []
    RGI_files.append('01_rgi60_Alaska')
    RGI_files.append('02_rgi60_WesternCanadaUS')
    RGI_files.append('03_rgi60_ArcticCanadaNorth')
    RGI_files.append('04_rgi60_ArcticCanadaSouth')
    RGI_files.append('05_rgi60_GreenlandPeriphery')
    RGI_files.append('06_rgi60_Iceland')
    RGI_files.append('07_rgi60_Svalbard')
    RGI_files.append('08_rgi60_Scandinavia')
    RGI_files.append('09_rgi60_RussianArctic')
    RGI_files.append('10_rgi60_NorthAsia')
    RGI_files.append('11_rgi60_CentralEurope')
    RGI_files.append('12_rgi60_CaucasusMiddleEast')
    RGI_files.append('13_rgi60_CentralAsia')
    RGI_files.append('14_rgi60_SouthAsiaWest')
    RGI_files.append('15_rgi60_SouthAsiaEast')
    RGI_files.append('16_rgi60_LowLatitudes')
    RGI_files.append('17_rgi60_SouthernAndes')
    RGI_files.append('18_rgi60_NewZealand')
    RGI_files.append('19_rgi60_AntarcticSubantarctic')
    # read input zipfile containing RGI shapefiles
    zs = zipfile.ZipFile(os.path.join(RGI_DIRECTORY,
        f'{RGI_files[RGI_REGION-1]}.zip'))
    dbf,prj,shp,shx = [io.BytesIO(zs.read(s)) for s in sorted(zs.namelist())
        if re.match(r'(.*?)\.(dbf|prj|shp|shx)$',s)]
    # read the shapefile and extract entities
    shape_input = shapefile.Reader(dbf=dbf, prj=prj, shp=shp, shx=shx,
        encodingErrors='ignore')
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    # extract the RGI entities
    poly_dict = {}
    for i,att in enumerate(shape_attributes):
        # extract latitude and longitude coordinates for entity
        points = np.array(shape_entities[i].points)
        # entities can have multiple parts
        parts = shape_entities[i].parts
        parts.append(len(points))
        # list object for coordinates (exterior and holes)
        poly_list = []
        # add each part to list
        for p1,p2 in zip(parts[:-1],parts[1:]):
            poly_list.append(list(zip(points[p1:p2,0],points[p1:p2,1])))
        # convert poly_list into Polygon object with holes
        poly_obj = Polygon(poly_list[0],poly_list[1:])
        # Valid Polygon may not possess overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        # add to dictionary based on RGI identifier
        poly_dict[att[0]] = poly_obj
    # close the zipfile
    zs.close()
    # return the dictionary of polygon objects and the input file
    return (poly_dict, RGI_files[RGI_REGION-1])

# PURPOSE: read ICESat-2 geolocated photon height data (ATL03)
# reduce to the Randolph Glacier Inventory
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

    # read data on rank 0
    if (comm.rank == 0):
        # read RGI for region and create shapely polygon objects
        poly_dict,RGI_file = load_glacier_inventory(args.directory,args.region)
    else:
        # create empty object for list of shapely objects
        poly_dict = None
        RGI_file = None

    # Broadcast Shapely polygon objects
    poly_dict = comm.bcast(poly_dict, root=0)
    RGI_file = comm.bcast(RGI_file, root=0)
    # RGI version and name
    RGI_VERSION,RGI_NAME = re.findall(r'\d_rgi(\d+)_(.*?)$',RGI_file).pop()
    # combined validity check for all beams
    valid_check = False

    # read each input beam within the file
    IS2_atl03_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains data
        # check in both the geolocation and heights groups
        try:
            fileID[gtx]['geolocation']['segment_id']
            fileID[gtx]['heights']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl03_beams.append(gtx)

    # copy variables for outputting to HDF5 file
    IS2_atl03_mask = {}
    IS2_atl03_fill = {}
    IS2_atl03_dims = {}
    IS2_atl03_mask_attrs = {}
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    IS2_atl03_mask['ancillary_data'] = {}
    IS2_atl03_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        # get each HDF5 variable
        IS2_atl03_mask['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        IS2_atl03_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl03_mask_attrs['ancillary_data'][key][att_name] = att_val

    # for each input beam within the file
    for gtx in sorted(IS2_atl03_beams):
        # output data dictionaries for beam
        IS2_atl03_mask[gtx] = dict(heights={},subsetting={})
        IS2_atl03_fill[gtx] = dict(heights={},subsetting={})
        IS2_atl03_dims[gtx] = dict(heights={},subsetting={})
        IS2_atl03_mask_attrs[gtx] = dict(heights={},subsetting={})

        # number of photon events
        n_pe, = fileID[gtx]['heights']['h_ph'].shape
        # check if there are less photon events than processes
        if (n_pe < comm.Get_size()):
            continue
        # define indices to run for specific process
        ind = np.arange(comm.Get_rank(), n_pe, comm.Get_size(), dtype=int)

        # extract delta time
        delta_time = fileID[gtx]['heights']['delta_time'][:]
        # extract lat/lon
        longitude = fileID[gtx]['heights']['lon_ph'][:]
        latitude = fileID[gtx]['heights']['lat_ph'][:]

        # convert reduced lat/lon to shapely multipoint object
        xy_point = MultiPoint(list(zip(longitude[ind],latitude[ind])))

        # create distributed intersection map for calculation
        distributed_map = np.zeros((n_pe),dtype=bool)
        distributed_RGIId = np.zeros((n_pe),dtype='|S14')
        # create empty intersection map array for receiving
        associated_map = np.zeros((n_pe),dtype=bool)
        associated_RGIId = np.zeros((n_pe),dtype='|S14')
        for key,poly_obj in poly_dict.items():
            # finds if points are encapsulated (within RGI polygon)
            int_test = poly_obj.intersects(xy_point)
            if int_test:
                # extract intersected points
                int_map = list(map(poly_obj.intersects,xy_point))
                int_indices, = np.nonzero(int_map)
                # set distributed_map indices to True for intersected points
                distributed_map[ind[int_indices]] = True
                distributed_RGIId[ind[int_indices]] = key
        # communicate output MPI matrices between ranks
        # operation is a logical "or" across the elements.
        comm.Allreduce(sendbuf=[distributed_map, MPI.BOOL], \
            recvbuf=[associated_map, MPI.BOOL], op=MPI.LOR)
        # operation is a element summation.
        comm.Allreduce(sendbuf=[distributed_RGIId, MPI.CHAR], \
            recvbuf=[associated_RGIId, MPI.CHAR], op=MPI.SUM)
        distributed_map = None
        distributed_RGIId = None
        # wait for all processes to finish calculation
        comm.Barrier()
        # add to validity check
        valid_check |= np.any(associated_map)

        # group attributes for beam
        IS2_atl03_mask_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl03_mask_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl03_mask_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl03_mask_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl03_mask_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl03_mask_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl03_mask_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        # group attributes for heights
        IS2_atl03_mask_attrs[gtx]['heights']['Description'] = ("Contains arrays of the "
            "parameters for each received photon.")
        IS2_atl03_mask_attrs[gtx]['heights']['data_rate'] = ("Data are stored at the "
            "photon detection rate.")

        # geolocation, time and segment ID
        # delta time
        IS2_atl03_mask[gtx]['heights']['delta_time'] = delta_time
        IS2_atl03_fill[gtx]['heights']['delta_time'] = None
        IS2_atl03_dims[gtx]['heights']['delta_time'] = None
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time'] = {}
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['standard_name'] = "time"
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['calendar'] = "standard"
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl03_mask_attrs[gtx]['heights']['delta_time']['coordinates'] = \
            "lat_ph lon_ph"
        # latitude
        IS2_atl03_mask[gtx]['heights']['latitude'] = latitude
        IS2_atl03_fill[gtx]['heights']['latitude'] = None
        IS2_atl03_dims[gtx]['heights']['latitude'] = ['delta_time']
        IS2_atl03_mask_attrs[gtx]['heights']['latitude'] = {}
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['units'] = "degrees_north"
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['long_name'] = "Latitude"
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['standard_name'] = "latitude"
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['description'] = ("Latitude of each "
            "received photon. Computed from the ECF Cartesian coordinates of the bounce point.")
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['valid_min'] = -90.0
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['valid_max'] = 90.0
        IS2_atl03_mask_attrs[gtx]['heights']['latitude']['coordinates'] = \
            "delta_time lon_ph"
        # longitude
        IS2_atl03_mask[gtx]['heights']['longitude'] = longitude
        IS2_atl03_fill[gtx]['heights']['longitude'] = None
        IS2_atl03_dims[gtx]['heights']['longitude'] = ['delta_time']
        IS2_atl03_mask_attrs[gtx]['heights']['longitude'] = {}
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['units'] = "degrees_east"
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['long_name'] = "Longitude"
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['standard_name'] = "longitude"
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['description'] = ("Longitude of each "
            "received photon. Computed from the ECF Cartesian coordinates of the bounce point.")
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['valid_min'] = -180.0
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['valid_max'] = 180.0
        IS2_atl03_mask_attrs[gtx]['heights']['longitude']['coordinates'] = \
            "delta_time lat_ph"

        # subsetting variables
        IS2_atl03_mask_attrs[gtx]['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce photon events to specific regions of interest.")
        IS2_atl03_mask_attrs[gtx]['subsetting']['data_rate'] = ("Data are stored at the photon "
            "detection rate.")

        # output mask to HDF5
        key = RGI_NAME.replace('_',' ')
        IS2_atl03_mask[gtx]['subsetting'][RGI_NAME] = associated_map
        IS2_atl03_fill[gtx]['subsetting'][RGI_NAME] = None
        IS2_atl03_dims[gtx]['subsetting'][RGI_NAME] = ['delta_time']
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME] = {}
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['contentType'] = "referenceInformation"
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['long_name'] = f'{key} Mask'
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['description'] = ('Mask calculated '
            f'using the {key} region from the Randolph Glacier Inventory.')
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['source'] = f'RGIv{RGI_VERSION}'
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl03_mask_attrs[gtx]['subsetting'][RGI_NAME]['coordinates'] = \
            "../heights/delta_time ../heights/lat_ph ../heights/lon_ph"
        # output RGI identifier
        IS2_atl03_mask[gtx]['subsetting']['RGIId'] = associated_RGIId
        IS2_atl03_fill[gtx]['subsetting']['RGIId'] = None
        IS2_atl03_dims[gtx]['subsetting']['RGIId'] = ['delta_time']
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId'] = {}
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['contentType'] = "referenceInformation"
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['long_name'] = "RGI Identifier"
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['description'] = ('Identification '
            f'code within version {RGI_VERSION} of the Randolph Glacier Inventory (RGI).')
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['source'] = f'RGIv{RGI_VERSION}'
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl03_mask_attrs[gtx]['subsetting']['RGIId']['coordinates'] = \
            "../heights/delta_time ../heights/lat_ph ../heights/lon_ph"
        # wait for all processes to finish calculation
        comm.Barrier()

    # parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and valid_check:
        # output HDF5 file with RGI masks
        fargs=(PRD,RGI_VERSION,RGI_NAME,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX)
        file_format='{0}_RGI{1}_{2}_{3}{4}{5}{6}{7}{8}_{9}{10}{11}_{12}_{13}{14}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
        # print file information
        logging.info(f'\t{output_file}')
        # write to output HDF5 file
        HDF5_ATL03_mask_write(IS2_atl03_mask, IS2_atl03_mask_attrs,
            CLOBBER=True, INPUT=os.path.basename(args.file),
            FILL_VALUE=IS2_atl03_fill, DIMENSIONS=IS2_atl03_dims,
            FILENAME=output_file)
        # change the permissions mode
        os.chmod(output_file, args.mode)
    # close the input file
    fileID.close()

# PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL03_mask_write(IS2_atl03_mask, IS2_atl03_attrs, INPUT=None,
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
    for k,v in IS2_atl03_mask['ancillary_data'].items():
        # Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl03_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # write each output beam
    beams = [k for k in IS2_atl03_mask.keys() if bool(re.match(r'gt\d[lr]',k))]
    for gtx in beams:
        fileID.create_group(gtx)
        # add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl03_attrs[gtx][att_name]

        # for each output data group
        for key in ['heights','subsetting']:
            # create group
            fileID[gtx].create_group(key)
            h5[gtx][key] = {}
            for att_name in ['Description','data_rate']:
                att_val = IS2_atl03_attrs[gtx][key][att_name]
                fileID[gtx][key].attrs[att_name] = att_val

            # all variables for group
            groupkeys=set(IS2_atl03_mask[gtx][key].keys())-set(['delta_time'])
            for k in ['delta_time',*sorted(groupkeys)]:
                # values and attributes
                v = IS2_atl03_mask[gtx][key][k]
                attrs = IS2_atl03_attrs[gtx][key][k]
                fillvalue = FILL_VALUE[gtx][key][k]
                # Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}'.format(gtx,key,k)
                if fillvalue:
                    h5[gtx][key][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype,
                        fillvalue=fillvalue, compression='gzip')
                else:
                    h5[gtx][key][k] = fileID.create_dataset(val,
                        np.shape(v), data=v, dtype=v.dtype,
                        compression='gzip')
                # create or attach dimensions for HDF5 variable
                if DIMENSIONS[gtx][key][k]:
                    # attach dimensions
                    for i,dim in enumerate(DIMENSIONS[gtx][key][k]):
                        h5[gtx][key][k].dims[i].attach_scale(
                            h5[gtx][key][dim])
                else:
                    # make dimension
                    h5[gtx][key][k].make_scale(k)
                # add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx][key][k].attrs[att_name] = att_val

    # HDF5 file title
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
    # add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    # add attributes for input ATL03 and ATL09 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    # find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl03_mask[gtx]['heights']['longitude']
        lat = IS2_atl03_mask[gtx]['heights']['latitude']
        delta_time = IS2_atl03_mask[gtx]['heights']['delta_time']
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
