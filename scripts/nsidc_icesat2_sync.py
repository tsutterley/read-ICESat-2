#!/usr/bin/env python
u"""
nsidc_icesat2_sync.py
Written by Tyler Sutterley (05/2022)

Acquires ICESat-2 datafiles from the National Snow and Ice Data Center (NSIDC)

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_sync.py --user=<username> --release=001 ATL06
    where <username> is your NASA Earthdata username

INPUTS:
    ATL03: Global Geolocated Photon Data
    ATL04: Normalized Relative Backscatter
    ATL06: Land Ice Height
    ATL07: Sea Ice Height
    ATL08: Land and Vegetation Height
    ATL09: Atmospheric Layer Characteristics
    ATL10: Sea Ice Freeboard
    ATL12: Ocean Surface Height
    ATL13: Inland Water Surface Height

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -W X, --password X: password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for alternative authentication
    -D X, --directory X: working data directory
    -Y X, --year X: years to sync
    -S X, --subdirectory X: specific subdirectories to sync
    -r X, --release X: ICESat-2 data release to sync
    -v X, --version X: ICESat-2 data version to sync
    -t X, --track X: ICESat-2 reference ground tracks to sync
    -g X, --granule X: ICESat-2 granule regions to sync
    -c X, --cycle=X: ICESat-2 cycles to sync
    -n X, --region X: ICESat-2 Named Region (ATL14/ATL15)
    -a X, --auxiliary: Sync ICESat-2 auxiliary files for each HDF5 file
    -I X, --index X: Input index of ICESat-2 files to sync
    -F, --flatten: Do not create subdirectories
    -P X, --np X: Number of processes to use in file downloads
    -T X, --timeout X: Timeout in seconds for blocking operations
    -R X, --retry X: Connection retry attempts
    -l, --log: output log of files downloaded
    -L, --list: print files to be transferred, but do not execute transfer
    -C, --clobber: Overwrite existing data in transfer
    --checksum: compare hashes to check if overwriting existing data
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 03/2022: use attempt login function to check credentials
    Updated 02/2022: added option to sync specific orbital cycles
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: set context for multiprocessing to fork child processes
        added option to compare checksums in order to overwrite data
        added a file length check to validate downloaded files
    Updated 05/2021: added options for connection timeout and retry attempts
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
        use regex backslash for comment special characters
    Updated 02/2021: added regular expression patterns for ATL11/14/15
    Updated 11/2020: nsidc_list will output a string for errors
    Updated 10/2020: using argparse to set parameters
    Updated 09/2020: use urllib imported in utilities
    Updated 08/2020: moved urllib opener to utilities. add credential check
        moved urllib directory listing to utilities
    Updated 07/2020: added option index to use a list of files to sync
    Updated 06/2020: added multiprocessing option for parallel download
    Updated 05/2020: added option netrc to use alternative authentication
        adjust regular expression to allow syncing of ATL07 sea ice products
        adjust regular expression for auxiliary products
    Updated 03/2020: added option flatten to not create subdirectories
    Updated 09/2019: added ssl context to urlopen headers
    Updated 07/2019: added options to sync specific granules, tracks and version
    Updated 06/2019: use strptime to extract last modified time of remote files
    Written 01/2019
"""
from __future__ import print_function

import sys
import os
import re
import io
import time
import shutil
import logging
import argparse
import posixpath
import traceback
import lxml.etree
import multiprocessing as mp
import icesat2_toolkit.utilities

#-- PURPOSE: sync the ICESat-2 elevation data from NSIDC
def nsidc_icesat2_sync(DIRECTORY, PRODUCTS, RELEASE, VERSIONS, GRANULES,
    TRACKS, YEARS=None, SUBDIRECTORY=None, CYCLES=None, REGION=None,
    AUXILIARY=False, INDEX=None, FLATTEN=False, TIMEOUT=None, RETRY=1,
    LOG=False, LIST=False, PROCESSES=0, CLOBBER=False, CHECKSUM=False,
    MODE=0o775):

    #-- check if directory exists and recursively create if not
    os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None

    #-- output of synchronized files
    if LOG:
        #-- format: NSIDC_ICESat-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NSIDC_ICESat-2_sync_{0}.log'.format(today)
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('ICESat-2 Data Sync Log ({0})'.format(today))

    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    #-- regular expression operator for finding files of a particular granule
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    if TRACKS:
        regex_track = r'|'.join(['{0:04d}'.format(T) for T in TRACKS])
    else:
        regex_track = r'\d{4}'
    if CYCLES:
        regex_cycle = r'|'.join(['{0:02d}'.format(C) for C in CYCLES])
    else:
        regex_cycle = r'\d{2}'
    regex_granule = r'|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = r'|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    regex_suffix = r'(.*?)' if AUXILIARY else r'(h5|nc)'
    default_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).{6}$')
    ATL11_pattern = r'({0})_({1})({2})_(\d{{2}})(\d{{2}})_({3})_({4})(.*?).{5}$'
    ATL1415_pattern = r'({0})_({1})_(\d{{2}})(\d{{2}})_({3})_({4})(.*?).{5}$'

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

    #-- build list of remote files, remote modification times and local files
    remote_files = []
    remote_mtimes = []
    local_files = []
    #-- build lists of files or use existing index file
    if INDEX:
        #-- read the index file, split at lines and remove all commented lines
        with open(os.path.expanduser(INDEX),'r') as f:
            files = [i for i in f.read().splitlines() if re.match(r'^(?!\#)',i)]
        #-- regular expression operator for extracting information from files
        rx = re.compile(r'(ATL\d{2})(-\d{2})?_(\d{4})(\d{2})(\d{2})(\d{2})'
            r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
        #-- for each line in the index
        for f in files:
            #-- extract parameters from ICESat-2 ATLAS HDF5 file
            PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(f).pop()
            #-- get directories from remote directory
            product_directory = '{0}.{1}'.format(PRD,RL)
            sd = '{0}.{1}.{2}'.format(YY,MM,DD)
            PATH = [HOST,'ATLAS',product_directory,sd]
            remote_dir = posixpath.join(HOST,'ATLAS',product_directory,sd)
            #-- local directory for product and subdirectory
            if FLATTEN:
                local_dir = os.path.expanduser(DIRECTORY)
            else:
                local_dir = os.path.join(DIRECTORY,product_directory,sd)
            #-- find ICESat-2 data file to get last modified time
            #-- find matching files (for granule, release, version, track)
            names,lastmod,error = icesat2_toolkit.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=f.strip())
            #-- print if file was not found
            if not names:
                logging.critical(error)
                continue
            #-- add to lists
            for colname,remote_mtime in zip(names,lastmod):
                #-- remote and local versions of the file
                remote_files.append(posixpath.join(remote_dir,colname))
                local_files.append(os.path.join(local_dir,colname))
                remote_mtimes.append(remote_mtime)
    else:
        #-- for each ICESat-2 product listed
        for p in PRODUCTS:
            logging.info('PRODUCT={0}'.format(p))
            #-- get directories from remote directory
            product_directory = '{0}.{1}'.format(p,RELEASE)
            PATH = [HOST,'ATLAS',product_directory]
            #-- compile regular expression operator
            if p in ('ATL11',):
                R1 = re.compile(ATL11_pattern.format(p,regex_track,
                    regex_granule,RELEASE,regex_version,regex_suffix))
            elif p in ('ATL14','ATL15'):
                regex_region = '|'.join(REGION)
                R1 = re.compile(ATL1415_pattern.format(p,regex_region,
                    RELEASE,regex_version,regex_suffix))
            else:
                R1 = re.compile(default_pattern.format(p,regex_track,
                    regex_cycle,regex_granule,RELEASE,regex_version,
                    regex_suffix))
            #-- read and parse request for subdirectories (find column names)
            remote_sub,_,error = icesat2_toolkit.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=R2,
                sort=True)
            #-- print if subdirectory was not found
            if not remote_sub:
                logging.critical(error)
                continue
            #-- for each remote subdirectory
            for sd in remote_sub:
                #-- local directory for product and subdirectory
                if FLATTEN:
                    local_dir = os.path.expanduser(DIRECTORY)
                else:
                    local_dir = os.path.join(DIRECTORY,product_directory,sd)
                logging.info("Building file list: {0}".format(sd))
                #-- find ICESat-2 data files
                PATH = [HOST,'ATLAS',product_directory,sd]
                remote_dir = posixpath.join(HOST,'ATLAS',product_directory,sd)
                #-- find matching files (for granule, release, version, track)
                names,lastmod,error = icesat2_toolkit.utilities.nsidc_list(PATH,
                    build=False,
                    timeout=TIMEOUT,
                    parser=parser,
                    pattern=R1,
                    sort=True)
                #-- print if file was not found
                if not names:
                    logging.critical(error)
                    continue
                #-- build lists of each ICESat-2 data file
                for colname,remote_mtime in zip(names,lastmod):
                    #-- remote and local versions of the file
                    remote_files.append(posixpath.join(remote_dir,colname))
                    local_files.append(os.path.join(local_dir,colname))
                    remote_mtimes.append(remote_mtime)

    #-- sync in series if PROCESSES = 0
    if (PROCESSES == 0):
        #-- sync each ICESat-2 data file
        for i,remote_file in enumerate(remote_files):
            #-- sync ICESat-2 files with NSIDC server
            args = (remote_file,remote_mtimes[i],local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY, LIST=LIST,
                CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)
            output = http_pull_file(*args, **kwds)
            #-- print the output string
            logging.info(output) if output else None
    else:
        #-- set multiprocessing start method
        ctx = mp.get_context("fork")
        #-- sync in parallel with multiprocessing Pool
        pool = ctx.Pool(processes=PROCESSES)
        #-- sync each ICESat-2 data file
        out = []
        for i,remote_file in enumerate(remote_files):
            #-- sync ICESat-2 files with NSIDC server
            args = (remote_file,remote_mtimes[i],local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY, LIST=LIST,
                CLOBBER=CLOBBER, CHECKSUM=CHECKSUM, MODE=MODE)
            out.append(pool.apply_async(multiprocess_sync,
                args=args,kwds=kwds))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()
        #-- print the output string
        for output in out:
            temp = output.get()
            logging.info(temp) if temp else None

    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: wrapper for running the sync program in multiprocessing mode
def multiprocess_sync(*args, **kwds):
    try:
        output = http_pull_file(*args, **kwds)
    except Exception as e:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        logging.critical('process id {0:d} failed'.format(os.getpid()))
        logging.error(traceback.format_exc())
    else:
        return output

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
#-- or if the checksums do not match between the files
def http_pull_file(remote_file, remote_mtime, local_file, TIMEOUT=None,
    RETRY=1, LIST=False, CLOBBER=False, CHECKSUM=False, MODE=0o775):
    #-- check if data directory exists and recursively create if not
    local_dir = os.path.dirname(local_file)
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
    #-- chunked transfer encoding size
    CHUNK = 16 * 1024
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
        #-- check if local version of file exists
    if CHECKSUM and os.access(local_file, os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        local_hash = icesat2_toolkit.utilities.get_hash(local_file)
        #-- generate checksum hash for remote file
        kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY, CHUNK=CHUNK)
        remote_buffer = retry_download(remote_file, **kwds)
        remote_hash = icesat2_toolkit.utilities.get_hash(remote_buffer)
        #-- compare checksums
        if (local_hash != remote_hash):
            TEST = True
            OVERWRITE = ' (checksums: {0} {1})'.format(local_hash,remote_hash)
    elif os.access(local_file, os.F_OK):
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
            #-- copy bytes or transfer file
            if CHECKSUM and os.access(local_file, os.F_OK):
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(local_file, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
            else:
                retry_download(remote_file, LOCAL=local_file,
                    TIMEOUT=TIMEOUT, RETRY=RETRY, CHUNK=CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)
        #-- return the output string
        return output

#-- PURPOSE: Try downloading a file up to a set number of times
def retry_download(remote_file, LOCAL=None, TIMEOUT=None, RETRY=1, CHUNK=0):
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
            #-- get the length of the remote file
            remote_length = int(response.headers['content-length'])
            #-- if copying to a local file
            if LOCAL:
                #-- copy contents to file using chunked transfer encoding
                #-- transfer should work with ascii and binary data formats
                with open(LOCAL, 'wb') as f:
                    shutil.copyfileobj(response, f, CHUNK)
                local_length = os.path.getsize(LOCAL)
            else:
                #-- copy remote file contents to bytesIO object
                remote_buffer = io.BytesIO()
                shutil.copyfileobj(response, remote_buffer, CHUNK)
                local_length = remote_buffer.getbuffer().nbytes
        except:
            pass
        else:
            #-- check that downloaded file matches original length
            if (local_length == remote_length):
                break
        #-- add to retry counter
        retry_counter += 1
    #-- check if maximum number of retries were reached
    if (retry_counter == RETRY):
        raise TimeoutError('Maximum number of retries reached')
    #-- return the bytesIO object
    if not LOCAL:
        #-- rewind bytesIO object to start
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Acquires ICESat-2 datafiles from the National Snow and
            Ice Data Center (NSIDC)
            """
    )
    #-- command line parameters
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('products',
        metavar='PRODUCTS', type=str, nargs='*', default=[],
        help='ICESat-2 products to sync')
    #-- NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years of data to sync
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to sync')
    #-- subdirectories of data to sync
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to sync')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='004',
        help='ICESat-2 Data Release')
    #-- ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 Data Version')
    #-- ICESat-2 granule region
    region = parser.add_mutually_exclusive_group(required=False)
    region.add_argument('--granule','-g',
        metavar='GRANULE', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
    #-- ICESat-2 orbital cycle
    parser.add_argument('--cycle','-c',
        type=int, nargs='+', default=None,
        help='ICESat-2 orbital cycles to sync')
    #-- ICESat-2 ATL14 and 15 named regions
    ATL1415_regions = ['AA','AK','CN','CS','GL','IC','SV','RU']
    region.add_argument('--region','-n',
        metavar='REGION', type=str, nargs='+',
        choices=ATL1415_regions, default=['AA','GL'],
        help='ICESat-2 Named Region (ATL14/ATL15)')
    #-- ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs)')
    #-- sync auxiliary files
    parser.add_argument('--auxiliary','-a',
        default=False, action='store_true',
        help='Sync ICESat-2 auxiliary files for each HDF5 file')
    #-- sync using files from an index
    group.add_argument('--index','-i',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Input index of ICESat-2 files to sync')
    #-- output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
    #-- run sync in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file downloads')
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
    parser.add_argument('--checksum',
        default=False, action='store_true',
        help='Compare hashes to check for overwriting existing data')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- NASA Earthdata hostname
    HOST = 'urs.earthdata.nasa.gov'
    #-- build a urllib opener for NASA Earthdata
    #-- check internet connection before attempting to run program
    opener = icesat2_toolkit.utilities.attempt_login(HOST,
        username=args.user, password=args.password,
        netrc=args.netrc)
    #-- check NASA earthdata credentials before attempting to run program
    nsidc_icesat2_sync(args.directory, args.products, args.release,
        args.version, args.granule, args.track, YEARS=args.year,
        SUBDIRECTORY=args.subdirectory, CYCLES=args.cycle,
        REGION=args.region, AUXILIARY=args.auxiliary, INDEX=args.index,
        FLATTEN=args.flatten, PROCESSES=args.np, TIMEOUT=args.timeout,
        RETRY=args.retry, LOG=args.log, LIST=args.list,
        CLOBBER=args.clobber, CHECKSUM=args.checksum, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
