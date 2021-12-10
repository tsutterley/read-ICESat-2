#!/usr/bin/env python
u"""
MPI_reduce_ICESat2_ATL11_RGI.py
Written by Tyler Sutterley (10/2021)

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
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 05/2021: print full path of output filename
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Written 12/2020
"""
from __future__ import print_function

import sys
import os
import re
import io
import h5py
import logging
import zipfile
import datetime
import argparse
import shapefile
import numpy as np
import collections
from mpi4py import MPI
from shapely.geometry import MultiPoint, Polygon
from icesat2_toolkit.convert_delta_time import convert_delta_time
import icesat2_toolkit.time
import icesat2_toolkit.utilities

#-- PURPOSE: keep track of MPI threads
def info(rank, size):
    logging.info('Rank {0:d} of {1:d}'.format(rank+1,size))
    logging.info('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        logging.info('parent process: {0:d}'.format(os.getppid()))
    logging.info('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: load zip file containing Randolph Glacier Inventory shapefiles
def load_glacier_inventory(RGI_DIRECTORY,RGI_REGION):
    #-- list of Randolph Glacier Inventory files
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
    #-- read input zipfile containing RGI shapefiles
    zs = zipfile.ZipFile(os.path.join(RGI_DIRECTORY,
        '{0}.zip'.format(RGI_files[RGI_REGION-1])))
    dbf,prj,shp,shx = [io.BytesIO(zs.read(s)) for s in sorted(zs.namelist())
        if re.match(r'(.*?)\.(dbf|prj|shp|shx)$',s)]
    #-- read the shapefile and extract entities
    shape_input = shapefile.Reader(dbf=dbf, prj=prj, shp=shp, shx=shx,
        encodingErrors='ignore')
    shape_entities = shape_input.shapes()
    shape_attributes = shape_input.records()
    #-- extract the RGI entities
    poly_dict = {}
    for i,att in enumerate(shape_attributes):
        #-- extract latitude and longitude coordinates for entity
        points = np.array(shape_entities[i].points)
        #-- entities can have multiple parts
        parts = shape_entities[i].parts
        parts.append(len(points))
        #-- list object for coordinates (exterior and holes)
        poly_list = []
        #-- add each part to list
        for p1,p2 in zip(parts[:-1],parts[1:]):
            poly_list.append(list(zip(points[p1:p2,0],points[p1:p2,1])))
        #-- convert poly_list into Polygon object with holes
        poly_obj = Polygon(poly_list[0],poly_list[1:])
        #-- Valid Polygon may not possess overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        #-- add to dictionary based on RGI identifier
        poly_dict[att[0]] = poly_obj
    #-- close the zipfile
    zs.close()
    #-- return the dictionary of polygon objects and the input file
    return (poly_dict, RGI_files[RGI_REGION-1])

#-- PURPOSE: read ICESat-2 annual land ice height data (ATL11) from NSIDC
#-- reduce to the Randolph Glacier Inventory
def main():
    #-- start MPI communicator
    comm = MPI.COMM_WORLD

    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Create masks for reducing ICESat-2 ATL11 annual land
            ice height data to the Randolph Glacier Inventory (RGI)
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = \
        icesat2_toolkit.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL11 file to run')
    #-- working data directory for location of RGI files
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- region of Randolph Glacier Inventory to run
    parser.add_argument('--region','-r',
        metavar='RGI', type=int, choices=range(1,20),
        help='region of Randolph Glacier Inventory to run')
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

    #-- create logger
    loglevel = logging.INFO if args.verbose else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- output module information for process
    info(comm.rank,comm.size)
    if (comm.rank == 0):
        logging.info('{0} -->'.format(args.file))

    #-- Open the HDF5 file for reading
    fileID = h5py.File(args.file, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(args.file)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed_)?(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_'
        r'(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(args.file).pop()

    #-- read data on rank 0
    if (comm.rank == 0):
        #-- read RGI for region and create shapely polygon objects
        poly_dict,RGI_file = load_glacier_inventory(args.directory,args.region)
    else:
        #-- create empty object for list of shapely objects
        poly_dict = None
        RGI_file = None

    #-- Broadcast Shapely polygon objects
    poly_dict = comm.bcast(poly_dict, root=0)
    RGI_file = comm.bcast(RGI_file, root=0)
    #-- RGI version and name
    RGI_VERSION,RGI_NAME = re.findall(r'\d_rgi(\d+)_(.*?)$',RGI_file).pop()
    #-- combined validity check for all beam pairs
    valid_check = False

    #-- read each input beam pair within the file
    IS2_atl11_pairs = []
    for ptx in [k for k in fileID.keys() if bool(re.match(r'pt\d',k))]:
        #-- check if subsetted beam contains reference points
        try:
            fileID[ptx]['ref_pt']
        except KeyError:
            pass
        else:
            IS2_atl11_pairs.append(ptx)

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
        IS2_atl11_mask['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        IS2_atl11_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl11_mask_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam pair within the file
    for ptx in sorted(IS2_atl11_pairs):
        #-- output data dictionaries for beam pair
        IS2_atl11_mask[ptx] = dict(subsetting=collections.OrderedDict())
        IS2_atl11_fill[ptx] = dict(subsetting={})
        IS2_atl11_dims[ptx] = dict(subsetting={})
        IS2_atl11_mask_attrs[ptx] = dict(subsetting={})

        #-- number of average segments and number of included cycles
        delta_time = fileID[ptx]['delta_time'][:].copy()
        n_points,n_cycles = np.shape(delta_time)
        #-- check if there are less segments than processes
        if (n_points < comm.Get_size()):
            continue

        #-- define indices to run for specific process
        ind = np.arange(comm.Get_rank(),n_points,comm.Get_size(),dtype=int)

        #-- convert reduced lat/lon to shapely multipoint object
        longitude = fileID[ptx]['longitude'][:].copy()
        latitude = fileID[ptx]['latitude'][:].copy()
        xy_point = MultiPoint(list(zip(longitude[ind],latitude[ind])))

        #-- create distributed intersection map for calculation
        distributed_map = np.zeros((n_points),dtype=bool)
        distributed_RGIId = np.zeros((n_points),dtype='|S14')
        #-- create empty intersection map array for receiving
        associated_map = np.zeros((n_points),dtype=bool)
        associated_RGIId = np.zeros((n_points),dtype='|S14')
        for key,poly_obj in poly_dict.items():
            #-- finds if points are encapsulated (within RGI polygon)
            int_test = poly_obj.intersects(xy_point)
            if int_test:
                #-- extract intersected points
                int_map = list(map(poly_obj.intersects,xy_point))
                int_indices, = np.nonzero(int_map)
                #-- set distributed_map indices to True for intersected points
                distributed_map[ind[int_indices]] = True
                distributed_RGIId[ind[int_indices]] = key
        #-- communicate output MPI matrices between ranks
        #-- operation is a logical "or" across the elements.
        comm.Allreduce(sendbuf=[distributed_map, MPI.BOOL], \
            recvbuf=[associated_map, MPI.BOOL], op=MPI.LOR)
        #-- operation is a element summation.
        comm.Allreduce(sendbuf=[distributed_RGIId, MPI.CHAR], \
            recvbuf=[associated_RGIId, MPI.CHAR], op=MPI.SUM)
        distributed_map = None
        distributed_RGIId = None
        #-- wait for all processes to finish calculation
        comm.Barrier()
        #-- add to validity check
        valid_check |= np.any(associated_map)

        #-- group attributes for beam pair
        IS2_atl11_mask_attrs[ptx]['description'] = ('Contains the primary science parameters for this '
            'data set')
        IS2_atl11_mask_attrs[ptx]['beam_pair'] = fileID[ptx].attrs['beam_pair']
        IS2_atl11_mask_attrs[ptx]['ReferenceGroundTrack'] = fileID[ptx].attrs['ReferenceGroundTrack']
        IS2_atl11_mask_attrs[ptx]['first_cycle'] = fileID[ptx].attrs['first_cycle']
        IS2_atl11_mask_attrs[ptx]['last_cycle'] = fileID[ptx].attrs['last_cycle']
        IS2_atl11_mask_attrs[ptx]['equatorial_radius'] = fileID[ptx].attrs['equatorial_radius']
        IS2_atl11_mask_attrs[ptx]['polar_radius'] = fileID[ptx].attrs['polar_radius']

        #-- geolocation, time and reference point
        #-- reference point
        IS2_atl11_mask[ptx]['ref_pt'] = fileID[ptx]['ref_pt'][:].copy()
        IS2_atl11_fill[ptx]['ref_pt'] = None
        IS2_atl11_dims[ptx]['ref_pt'] = None
        IS2_atl11_mask_attrs[ptx]['ref_pt'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['ref_pt']['units'] = "1"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['contentType'] = "referenceInformation"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['long_name'] = "Reference point number"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['source'] = "ATL06"
        IS2_atl11_mask_attrs[ptx]['ref_pt']['description'] = ("The reference point is the 7 "
            "digit segment_id number corresponding to the center of the ATL06 data used for "
            "each ATL11 point.  These are sequential, starting with 1 for the first segment "
            "after an ascending equatorial crossing node.")
        IS2_atl11_mask_attrs[ptx]['ref_pt']['coordinates'] = \
            "delta_time latitude longitude"

        #-- cycle_number
        IS2_atl11_mask[ptx]['cycle_number'] = fileID[ptx]['cycle_number'][:].copy()
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
        IS2_atl11_mask[ptx]['delta_time'] = fileID[ptx]['delta_time'][:].copy()
        IS2_atl11_fill[ptx]['delta_time'] = fileID[ptx]['delta_time'].attrs['_FillValue']
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
        IS2_atl11_mask[ptx]['latitude'] = fileID[ptx]['latitude'][:].copy()
        IS2_atl11_fill[ptx]['latitude'] = fileID[ptx]['latitude'].attrs['_FillValue']
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
        IS2_atl11_mask[ptx]['longitude'] = fileID[ptx]['longitude'][:].copy()
        IS2_atl11_fill[ptx]['longitude'] = fileID[ptx]['longitude'].attrs['_FillValue']
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

        #-- output mask to HDF5
        key = RGI_NAME.replace('_',' ')
        IS2_atl11_mask[ptx]['subsetting'][RGI_NAME] = associated_map
        IS2_atl11_fill[ptx]['subsetting'][RGI_NAME] = None
        IS2_atl11_dims[ptx]['subsetting'][RGI_NAME] = ['ref_pt']
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['contentType'] = "referenceInformation"
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['long_name'] = '{0} Mask'.format(key)
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['description'] = ('Mask calculated '
            'using the {0} region from the Randolph Glacier Inventory.').format(key)
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['source'] = \
            'RGIv{0}'.format(RGI_VERSION)
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl11_mask_attrs[ptx]['subsetting'][RGI_NAME]['coordinates'] = \
            "../ref_pt ../delta_time ../latitude ../longitude"

        #-- output RGI identifier
        IS2_atl11_mask[ptx]['subsetting']['RGIId'] = associated_RGIId
        IS2_atl11_fill[ptx]['subsetting']['RGIId'] = None
        IS2_atl11_dims[ptx]['subsetting']['RGIId'] = ['ref_pt']
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId'] = collections.OrderedDict()
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['contentType'] = "referenceInformation"
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['long_name'] = "RGI Identifier"
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['description'] = ('Identification '
            'code within version {0} of the Randolph Glacier Inventory (RGI).').format(RGI_VERSION)
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['source'] = \
            'RGIv{0}'.format(RGI_VERSION)
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl11_mask_attrs[ptx]['subsetting']['RGIId']['coordinates'] = \
            "../ref_pt ../delta_time ../latitude ../longitude"

        #-- wait for all processes to finish calculation
        comm.Barrier()

    #-- parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and valid_check:
        #-- output HDF5 file with RGI masks
        fargs = (PRD,RGI_VERSION,RGI_NAME,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX)
        file_format = '{0}_RGI{1}_{2}_{3}{4}_{5}{6}_{7}_{8}{9}.h5'
        output_file = os.path.join(DIRECTORY,file_format.format(*fargs))
        #-- print file information
        logging.info('\t{0}'.format(output_file))
        #-- write to output HDF5 file
        HDF5_ATL11_mask_write(IS2_atl11_mask, IS2_atl11_mask_attrs,
            CLOBBER=True, INPUT=os.path.basename(args.file),
            FILL_VALUE=IS2_atl11_fill, DIMENSIONS=IS2_atl11_dims,
            FILENAME=output_file)
        #-- change the permissions mode
        os.chmod(output_file, args.mode)
    #-- close the input file
    fileID.close()

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

#-- run main program
if __name__ == '__main__':
    main()
