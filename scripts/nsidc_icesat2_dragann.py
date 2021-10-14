#!/usr/bin/env python
u"""
nsidc_icesat2_dragann.py
Written by Tyler Sutterley (10/2021)

Acquires the ATL03 geolocated photon height product and appends the
    ATL08 DRAGANN classifications from NSIDC

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_dragann.py --user=<username> --release=003
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for alternative authentication
    -D X, --directory X: working data directory
    -Y X, --year X: years to sync
    -S X, --subdirectory X: specific subdirectories to sync
    -r X, --release X: ICESat-2 data release to sync
    -v X, --version X: ICESat-2 data version to sync
    -t X, --track X: ICESat-2 reference ground tracks to sync
    -g X, --granule X: ICESat-2 granule regions to sync
    -F, --flatten: Do not create subdirectories
    -T X, --timeout X: Timeout in seconds for blocking operations
    -R X, --retry X: Connection retry attempts
    -l, --log: output log of files downloaded
    -L, --list: print files to be transferred, but do not execute transfer
    -C, --clobber: Overwrite existing data in transfer
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
        http://docs.h5py.org/en/stable/index.html
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files
    read_ICESat2_ATL03.py: reads ICESat-2 global geolocated photon data files

UPDATE HISTORY:
    Updated 10/2021: using python logging for handling verbose output
    Updated 05/2021: added options for connection timeout and retry attempts
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Written 02/2021
"""
from __future__ import print_function

import sys
import os
import io
import re
import h5py
import netrc
import shutil
import getpass
import logging
import argparse
import builtins
import posixpath
import traceback
import lxml.etree
import numpy as np
import calendar, time
import multiprocessing as mp
import icesat2_toolkit.utilities
from icesat2_toolkit.read_ICESat2_ATL03 import find_HDF5_ATL03_beams

#-- PURPOSE: sync ATL03 geolocated photon height products and appends the
#-- ATL08 DRAGANN classifications from NSIDC
def nsidc_icesat2_dragann(DIRECTORY, RELEASE, VERSIONS, GRANULES, TRACKS,
    YEARS=None, SUBDIRECTORY=None, FORMAT=None, CHUNKS=None, FLATTEN=False,
    TIMEOUT=None, RETRY=1, LOG=False, LIST=False, CLOBBER=False, MODE=None):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- output of synchronized files
    if LOG:
        #-- format: NSIDC_ICESat-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NSIDC_ICESat-2_sync_{0}.log'.format(today)
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('ICESat-2 DRAGANN Sync Log ({0})'.format(today))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    #-- regular expression operator for finding files of a particular granule
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    regex_suffix = '(h5)'
    remote_regex_pattern=(r'({0})_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})(\d{{2}})({2})_({3})_({4})(.*?).{5}$')
    #-- regular expression operator for finding subdirectories
    if SUBDIRECTORY:
        #-- Sync particular subdirectories for product
        R2 = re.compile(r'('+'|'.join(SUBDIRECTORY)+')', re.VERBOSE)
    elif YEARS:
        #-- Sync particular years for product
        regex_pattern = '|'.join('{0:d}'.format(y) for y in YEARS)
        R2 = re.compile(r'({0}).(\d+).(\d+)'.format(regex_pattern), re.VERBOSE)
    else:
        #-- Sync all available subdirectories for product
        R2 = re.compile(r'(\d+).(\d+).(\d+)', re.VERBOSE)

    #-- get directories from remote directory
    atl03_directory = '{0}.{1}'.format('ATL03',RELEASE)
    atl08_directory = '{0}.{1}'.format('ATL08',RELEASE)
    PATH = [HOST,'ATLAS',atl08_directory]
    #-- compile regular expression operator
    args=('ATL08',regex_track,regex_granule,RELEASE,regex_version,regex_suffix)
    R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)
    #-- read and parse request for subdirectories (find column names)
    remote_sub,_,error = icesat2_toolkit.utilities.nsidc_list(PATH,
        build=False,
        timeout=TIMEOUT,
        parser=parser,
        pattern=R2,
        sort=True)
    #-- for each remote subdirectory
    for sd in remote_sub:
        #-- local directory for product and subdirectory
        if FLATTEN:
            local_dir = os.path.expanduser(DIRECTORY)
        else:
            local_dir = os.path.join(DIRECTORY,atl03_directory,sd)
        #-- check if data directory exists and recursively create if not
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- find ICESat-2 data files
        PATH = [HOST,'ATLAS',atl08_directory,sd]
        #-- find matching files (for granule, release, version, track)
        atl08s,lastmod,error = icesat2_toolkit.utilities.nsidc_list(PATH,
            build=False,timeout=TIMEOUT,parser=parser,pattern=R1,sort=True)
        #-- print if file was not found
        if not atl08s:
            logging.critical(error)
            continue
        #-- build lists of each ICESat-2 data file
        for i,atl08 in enumerate(atl08s):
            #-- extract base parameters
            PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX,SFX = \
                R1.findall(atl08).pop()
            #-- compile regular expression operator
            args = ('ATL03',TRK,GRAN,RL,regex_version,regex_suffix)
            R3 = re.compile(remote_regex_pattern.format(*args))
            PATH = [HOST,'ATLAS',atl03_directory,sd]
            #-- find associated ATL03 files
            atl03s,lastmod,error=icesat2_toolkit.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=R3,
                sort=True)
            #-- remote and local versions of the file
            for atl03,remote_mtime in zip(atl03s,lastmod):
                #-- sync ICESat-2 ATL03 files with NSIDC server
                remote_dir = posixpath.join(HOST,'ATLAS',atl03_directory,sd)
                remote_file = posixpath.join(remote_dir,atl03)
                local_file = os.path.join(local_dir,atl03)
                #-- download ATL03 file
                args = (remote_file, remote_mtime, local_file)
                kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY,
                    LIST=LIST, CLOBBER=CLOBBER)
                out = http_pull_file(*args, **kwds)
                logging.info(out) if out else None
                #-- append ATL08 dragann classifications
                PATH = [HOST,'ATLAS',atl08_directory,sd,atl08]
                logging.info(posixpath.join(*PATH))
                remote_buffer,_ = icesat2_toolkit.utilities.from_nsidc(PATH,
                    build=False,
                    timeout=TIMEOUT)
                #-- for each beam in the ATL03 file
                for gtx in find_HDF5_ATL03_beams(local_file):
                    #-- open ATL03 file in append mode
                    fileID = h5py.File(local_file, 'a')
                    #-- check if DRAGANN variables are already appended
                    if 'd_flag' in fileID[gtx]['heights'].keys():
                        #-- close the file and continue
                        fileID.close()
                        continue
                    #-- ATL03 20 meter segment id
                    segment_id=fileID[gtx]['geolocation']['segment_id'][:]
                    #-- first photon in the segment (0-based indexing)
                    ph_index_beg=fileID[gtx]['geolocation']['ph_index_beg'][:]-1
                    #-- photon event level delta time
                    delta_time = fileID[gtx]['heights']['delta_time']
                    #-- number of photon events in the beam group
                    n_pe, = delta_time.shape
                    #-- extract dragann classifiers for beam
                    mds,attrs = extract_dragann_classification(remote_buffer,
                        gtx,segment_id,ph_index_beg,n_pe)
                    for k,v in mds.items():
                        #-- create HDF5 variable
                        val = posixpath.join(gtx,'heights',k)
                        h5 = fileID.create_dataset(val, np.shape(v),
                            data=v, dtype=v.dtype, compression='gzip')
                        h5.dims[0].attach_scale(delta_time)
                        #-- add HDF5 variable attributes
                        for att_name,att_val in attrs[k].items():
                            h5.attrs[att_name] = att_val
                    #-- close the ATL03 file
                    fileID.close()
                #-- keep remote modification time of file and local access time
                os.utime(local_file,(os.stat(local_file).st_atime,remote_mtime))
                #-- change the permissions mode
                os.chmod(local_file, MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
    TIMEOUT=None, RETRY=1, LIST=False, CLOBBER=False):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if os.access(local_file, os.F_OK):
        #-- check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- output string for printing files transferred
        output = '{0} -->\n\t{1}{2}\n'.format(remote_file,local_file,OVERWRITE)
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- attempt to download up to the number of retries
            retry_counter = 0
            while (retry_counter < RETRY):
                #-- attempt to retrieve file from https server
                try:
                    #-- Create and submit request.
                    #-- There are a range of exceptions that can be thrown here
                    #-- including HTTPError and URLError.
                    request=icesat2_toolkit.utilities.urllib2.Request(remote_file)
                    response=icesat2_toolkit.utilities.urllib2.urlopen(request,
                        timeout=TIMEOUT)
                    #-- copy contents to file using chunked transfer encoding
                    #-- transfer should work with ascii and binary data formats
                    with open(local_file, 'wb') as f:
                        shutil.copyfileobj(response, f, CHUNK)
                except:
                    pass
                else:
                    break
                #-- add to retry counter
                retry_counter += 1
            #-- check if maximum number of retries were reached
            if (retry_counter == RETRY):
                raise TimeoutError('Maximum number of retries reached')
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        #-- return the output string
        return output

#-- PURPOSE: Extract photon event classification from ATL08
def extract_dragann_classification(buffer,gtx,segment_id,ph_index_beg,n_pe):
    #-- Open the HDF5 file for reading
    if isinstance(buffer, io.IOBase):
        fileID = h5py.File(buffer, 'r')
    else:
        fileID = h5py.File(os.path.expanduser(buffer), 'r')
    #-- allocate for output classifications
    output = {}
    #-- default is noise classification for all photon events
    output['d_flag'] = np.zeros((n_pe),dtype=np.int8)
    output['classed_pc_flag'] = np.zeros((n_pe),dtype=np.int8)
    attrs = dict(d_flag={},classed_pc_flag={})
    #-- attributes for the output remapped variables
    #-- DRAGANN classification
    attrs['d_flag']['units'] = 1
    attrs['d_flag']['long_name'] = "DRAGANN flag"
    attrs['d_flag']['description'] = ("Flag indicating the labeling of "
        "DRAGANN noise filtering for a given photon.")
    attrs['d_flag']['flag_values'] = [0,1]
    attrs['d_flag']['valid_min'] = 0
    attrs['d_flag']['valid_max'] = 1
    attrs['d_flag']['source'] = "Land ATBD section 2.3.5"
    attrs['d_flag']['contentType'] = "qualityInformation"
    attrs['d_flag']['coordinates'] = "delta_time lat_ph lon_ph"
    #-- ATL08 photon event classification
    attrs['classed_pc_flag']['units'] = 1
    attrs['classed_pc_flag']['long_name'] = "Photon Land Classification"
    attrs['classed_pc_flag']['description'] = ("Land Vegetation ATBD "
        "classification flag for each photon as either noise, ground, canopy, "
        "and top of canopy. 0 = noise,  1 = ground, 2 = canopy, or "
        "3 = top of canopy.")
    attrs['classed_pc_flag']['flag_values'] = [0,1,2,3]
    attrs['classed_pc_flag']['valid_min'] = 0
    attrs['classed_pc_flag']['valid_max'] = 3
    attrs['classed_pc_flag']['source'] = "Land ATBD section 4.10"
    attrs['classed_pc_flag']['contentType'] = "qualityInformation"
    attrs['classed_pc_flag']['coordinates'] = "delta_time lat_ph lon_ph"
    #-- check if subsetted beam contains photon event classifications
    try:
        fileID[gtx]['signal_photons']['delta_time']
    except KeyError:
        pass
    else:
        #-- segment id of the ATL08 classified photons
        ph_segment_id = fileID[gtx]['signal_photons']['ph_segment_id'][:]
        #-- index of the ATL08 classified photons within the ATL03 20m segments
        classed_pc_indx = fileID[gtx]['signal_photons']['classed_pc_indx'][:]-1
        #-- dragann flag and land classifications
        d_flag = fileID[gtx]['signal_photons']['d_flag']
        classed_pc_flag = fileID[gtx]['signal_photons']['classed_pc_flag']
        #-- for each ATL03 20 meter segment
        for s1,i1 in zip(segment_id,ph_index_beg):
            #-- indices within the ATL08 signal photons group
            i2, = np.nonzero(ph_segment_id == s1)
            #-- full index within the ATL03 file
            indx = i1 + classed_pc_indx[i2]
            #-- remapped dragann and land ice classifications
            output['d_flag'][indx] = d_flag[i2]
            output['classed_pc_flag'][indx] = classed_pc_flag[i2]
    #-- close the ATL08 HDF5 file
    fileID.close()
    #-- return the output variables and attributes
    return (output,attrs)

#-- Main program that calls nsidc_icesat2_dragann()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Acquires the ATL03 geolocated photon height product
            and appends the ATL08 DRAGANN classifications from NSIDC
            """
    )
    #-- command line parameters
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years of data to run
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to run')
    #-- subdirectories of data to run
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to run')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='004',
        help='ICESat-2 Data Release')
    #-- ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 Data Version')
    #-- ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
    #-- ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs)')
    #-- output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
    #-- connection timeout and number of retry attempts
    parser.add_argument('--timeout','-T',
        type=int, default=120,
        help='Timeout in seconds for blocking operations')
    parser.add_argument('--retry','-R',
        type=int, default=5,
        help='Connection retry attempts')
    #-- Output log file in form
    #-- NSIDC_IceSat-2_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    #-- clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    #-- permissions mode of the output files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- NASA Earthdata hostname
    HOST = 'urs.earthdata.nasa.gov'
    #-- get authentication
    try:
        args.user,_,PASSWORD = netrc.netrc(args.netrc).authenticators(HOST)
    except:
        #-- check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = 'Username for {0}: '.format(HOST)
            args.user = builtins.input(prompt)
        #-- enter password securely from command-line
        prompt = 'Password for {0}@{1}: '.format(args.user,HOST)
        PASSWORD = getpass.getpass(prompt)

    #-- build a urllib opener for NSIDC
    #-- Add the username and password for NASA Earthdata Login system
    icesat2_toolkit.utilities.build_opener(args.user,PASSWORD)

    #-- check internet connection before attempting to run program
    #-- check NASA earthdata credentials before attempting to run program
    if icesat2_toolkit.utilities.check_credentials():
        nsidc_icesat2_dragann(args.directory, args.release, args.version,
            args.granule, args.track, YEARS=args.year,
            SUBDIRECTORY=args.subdirectory, FLATTEN=args.flatten,
            TIMEOUT=args.timeout, RETRY=args.retry, LOG=args.log,
            LIST=args.list, CLOBBER=args.clobber, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
