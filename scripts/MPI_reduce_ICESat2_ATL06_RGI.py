#!/usr/bin/env python
u"""
MPI_reduce_ICESat2_ATL06_RGI.py
Written by Tyler Sutterley (10/2020)

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
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    convert_calendar_decimal.py: converts from calendar date to decimal year
    time.py: Utilities for calculating time operations
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2020: using argparse to set parameters
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 10/2019: changing Y/N flags to True/False
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
import h5py
import zipfile
import datetime
import argparse
import shapefile
import numpy as np
from mpi4py import MPI
from shapely.geometry import MultiPoint, Polygon
from icesat2_toolkit.convert_julian import convert_julian
from icesat2_toolkit.convert_delta_time import convert_delta_time

#-- PURPOSE: keep track of MPI threads
def info(rank, size):
    print('Rank {0:d} of {1:d}'.format(rank+1,size))
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

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

#-- PURPOSE: read ICESat-2 data from NSIDC or MPI_ICESat2_ATL03.py
#-- reduce to the Randolph Glacier Inventory
def main():
    #-- start MPI communicator
    comm = MPI.COMM_WORLD

    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Create masks for reducing ICESat-2 land ice data to
            the Randolph Glacier Inventory (RGI)
            """
    )
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='ICESat-2 ATL06 file to run')
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
    args = parser.parse_args()

    #-- output module information for process
    if args.verbose:
        info(comm.rank,comm.size)
    if args.verbose and (comm.rank==0):
        print('{0} -->'.format(args.file))

    #-- Open the HDF5 file for reading
    fileID = h5py.File(args.file, 'r', driver='mpio', comm=comm)
    DIRECTORY = os.path.dirname(args.file)
    #-- extract parameters from ICESat-2 ATLAS HDF5 file name
    rx = re.compile(r'(processed)?(ATL\d{2})_(\d{4})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(args.file).pop()

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

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = fileID['ancillary_data']['atlas_sdp_gps_epoch'][:]

    #-- copy variables for outputting to HDF5 file
    IS2_atl06_mask = {}
    IS2_atl06_fill = {}
    IS2_atl06_mask_attrs = {}
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    IS2_atl06_mask['ancillary_data'] = {}
    IS2_atl06_mask_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl06_mask['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        IS2_atl06_mask_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl06_mask_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl06_beams):
        #-- number of segments
        segment_id = fileID[gtx]['land_ice_segments']['segment_id'][:]
        n_seg, = fileID[gtx]['land_ice_segments']['segment_id'].shape
        #-- invalid value for beam
        fv = fileID[gtx]['land_ice_segments']['h_li'].fillvalue
        #-- check if there are less segments than processes
        if (n_seg < comm.Get_size()):
            continue

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

        #-- convert reduced lat/lon to shapely multipoint object
        xy_point = MultiPoint(list(zip(longitude[ind],latitude[ind])))

        #-- create distributed intersection map for calculation
        distributed_map = np.zeros((n_seg),dtype=np.bool)
        distributed_RGIId = np.zeros((n_seg),dtype='|S14')
        #-- create empty intersection map array for receiving
        associated_map = np.zeros((n_seg),dtype=np.bool)
        associated_RGIId = np.zeros((n_seg),dtype='|S14')
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

        #-- output data dictionaries for beam
        IS2_atl06_mask[gtx] = dict(land_ice_segments={})
        IS2_atl06_fill[gtx] = dict(land_ice_segments={})
        IS2_atl06_mask_attrs[gtx] = dict(land_ice_segments={})

        #-- group attributes for beam
        IS2_atl06_mask_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl06_mask_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl06_mask_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl06_mask_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl06_mask_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl06_mask_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl06_mask_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        #-- group attributes for land_ice_segments
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        #-- geolocation, time and segment ID
        #-- delta time
        IS2_atl06_mask[gtx]['land_ice_segments']['delta_time'] = delta_time
        IS2_atl06_fill[gtx]['land_ice_segments']['delta_time'] = delta_time.fill_value
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
        #-- latitude
        IS2_atl06_mask[gtx]['land_ice_segments']['latitude'] = latitude
        IS2_atl06_fill[gtx]['land_ice_segments']['latitude'] = latitude.fill_value
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
        #-- longitude
        IS2_atl06_mask[gtx]['land_ice_segments']['longitude'] = longitude
        IS2_atl06_fill[gtx]['land_ice_segments']['longitude'] = longitude.fill_value
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
        #-- segment ID
        IS2_atl06_mask[gtx]['land_ice_segments']['segment_id'] = segment_id
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

        #-- subsetting variables
        IS2_atl06_mask[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['Description'] = ("The subsetting group "
            "contains parameters used to reduce land ice segments to specific regions of interest.")
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        #-- output mask to HDF5
        key = RGI_NAME.replace('_',' ')
        IS2_atl06_mask[gtx]['land_ice_segments']['subsetting'][RGI_NAME] = associated_map
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['contentType'] = "referenceInformation"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['long_name'] = '{0} Mask'.format(key)
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['description'] = ('Mask calculated '
            'using the {0} region from the Randolph Glacier Inventory.').format(key)
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['source'] = \
            'RGIv{0}'.format(RGI_VERSION)
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting'][RGI_NAME]['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- output RGI identifier
        IS2_atl06_mask[gtx]['land_ice_segments']['subsetting']['RGIId'] = associated_RGIId
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId'] = {}
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['contentType'] = "referenceInformation"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['long_name'] = "RGI Identifier"
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['description'] = ('Identification '
            'code within version {0} of the Randolph Glacier Inventory (RGI).').format(RGI_VERSION)
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['source'] = \
            'RGIv{0}'.format(RGI_VERSION)
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['reference'] = \
            'https://www.glims.org/RGI/'
        IS2_atl06_mask_attrs[gtx]['land_ice_segments']['subsetting']['RGIId']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- wait for all processes to finish calculation
        comm.Barrier()

    #-- parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0) and associated_map.any():
        #-- output HDF5 file with RGI masks
        arg=(PRD,RGI_VERSION,RGI_NAME,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX)
        file_format='{0}_RGI{1}_{2}_{3}{4}{5}{6}{7}{8}_{9}{10}{11}_{12}_{13}{14}.h5'
        #-- print file information
        if args.verbose:
            print('\t{0}'.format(file_format.format(*arg)))
        #-- write to output HDF5 file
        HDF5_ATL06_mask_write(IS2_atl06_mask, IS2_atl06_mask_attrs, CLOBBER=True,
            INPUT=os.path.basename(args.file), FILL_VALUE=IS2_atl06_fill,
            FILENAME=os.path.join(DIRECTORY,file_format.format(*arg)))
        #-- change the permissions mode
        os.chmod(os.path.join(DIRECTORY,file_format.format(*arg)), args.mode)
    #-- close the input file
    fileID.close()

#-- PURPOSE: outputting the masks for ICESat-2 data to HDF5
def HDF5_ATL06_mask_write(IS2_atl06_mask, IS2_atl06_attrs, INPUT=None,
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
    for k,v in IS2_atl06_mask['ancillary_data'].items():
        #-- Defining the HDF5 dataset variables
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype, compression='gzip')
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl06_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- write each output beam
    beams = [k for k in IS2_atl06_mask.keys() if bool(re.match(r'gt\d[lr]',k))]
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
        v = IS2_atl06_mask[gtx]['land_ice_segments']['segment_id']
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
            v = IS2_atl06_mask[gtx]['land_ice_segments'][k]
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

        #-- add to subsetting variables
        key = 'subsetting'
        fileID[gtx]['land_ice_segments'].create_group(key)
        h5[gtx]['land_ice_segments'][key] = {}
        for att_name in ['Description','data_rate']:
            att_val=IS2_atl06_attrs[gtx]['land_ice_segments'][key][att_name]
            fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
        for k,v in IS2_atl06_mask[gtx]['land_ice_segments'][key].items():
            #-- attributes
            attrs = IS2_atl06_attrs[gtx]['land_ice_segments'][key][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}/{3}'.format(gtx,'land_ice_segments',key,k)
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
    fileID.attrs['summary'] = ('Subsetting masks for ice-sheets segments '
        'needed to interpret and assess the quality of land height estimates.')
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
    #-- add attributes for input ATL06 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in beams:
        lon = IS2_atl06_mask[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl06_mask[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl06_mask[gtx]['land_ice_segments']['delta_time']
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
