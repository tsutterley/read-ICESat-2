#!/usr/bin/env python
u"""
nsidc_icesat2_sync.py
Written by Tyler Sutterley (08/2020)

Program to acquire ICESat-2 datafiles from NSIDC server:
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
    -U X, --user=X: username for NASA Earthdata Login
    -N X, --netrc=X: path to .netrc file for alternative authentication
    -D X, --directory: working data directory
    -Y X, --year=X: years to sync separated by commas
    -S X, --subdirectory=X: subdirectories to sync separated by commas
    --release=X: ICESat-2 data release to sync
    --version=X: ICESat-2 data version to sync
    --track=X: ICESat-2 reference ground tracks to sync
    --granule=X: ICESat-2 granule regions to sync
    --auxiliary: Sync ICESat-2 auxiliary files for each HDF5 file
    -I X, --index=X: Input index of ICESat-2 files to sync
    -F, --flatten: Do not create subdirectories
    -P X, --np=X: Number of processes to use in file downloads
    -M X, --mode=X: Local permissions mode of the directories and files synced
    -l, --log: output log of files downloaded
    -L, --list: print files to be transferred, but do not execute transfer
    -C, --clobber: Overwrite existing data in transfer

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

PROGRAM DEPENDENCIES:
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
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
import time
import netrc
import getopt
import shutil
import getpass
import builtins
import posixpath
import traceback
import lxml.etree
import numpy as np
import multiprocessing as mp
import icesat2_toolkit.utilities
if sys.version_info[0] == 2:
    import urllib2
else:
    import urllib.request as urllib2

#-- PURPOSE: sync the ICESat-2 elevation data from NSIDC
def nsidc_icesat2_sync(ddir, PRODUCTS, RELEASE, VERSIONS, GRANULES, TRACKS,
    YEARS=None, SUBDIRECTORY=None, AUXILIARY=False, INDEX=None, FLATTEN=False,
    LOG=False, LIST=False, PROCESSES=0, MODE=None, CLOBBER=False):

    #-- check if directory exists and recursively create if not
    os.makedirs(ddir,MODE) if not os.path.exists(ddir) else None

    #-- output of synchronized files
    if LOG:
        #-- format: NSIDC_IceBridge_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'NSIDC_IceSat-2_sync_{0}.log'.format(today)
        fid = open(os.path.join(ddir,LOGFILE),'w')
        print('ICESat-2 Data Sync Log ({0})'.format(today), file=fid)
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    #-- remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    #-- regular expression operator for finding files of a particular granule
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    regex_suffix = '(.*?)' if AUXILIARY else '(h5)'
    remote_regex_pattern=(r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
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

    #-- build list of remote files, remote modification times and local files
    remote_files = []
    remote_mtimes = []
    local_files = []
    #-- build lists of files or use existing index file
    if INDEX:
        #-- read the index file, split at lines and remove all commented lines
        with open(os.path.expanduser(INDEX),'r') as f:
            files = [i for i in f.read().splitlines() if re.match('^(?!#)',i)]
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
                local_dir = os.path.expanduser(ddir)
            else:
                local_dir = os.path.join(ddir,product_directory,sd)
            #-- check if data directory exists and recursively create if not
            os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
            #-- find ICESat-2 data file to get last modified time
            #-- find matching files (for granule, release, version, track)
            colnames,collastmod = icesat2_toolkit.utilities.nsidc_list(PATH,
                build=False,timeout=120,parser=parser,pattern=f.strip())
            #-- print if file was not found
            if not colnames:
                print('{0} not found on {1}'.format(f,remote_dir),file=fid)
            #-- add to lists
            for colname,remote_mtime in zip(colnames,collastmod):
                #-- remote and local versions of the file
                remote_files.append(posixpath.join(remote_dir,colname))
                local_files.append(os.path.join(local_dir,colname))
                remote_mtimes.append(remote_mtime)
    else:
        #-- for each ICESat-2 product listed
        for p in PRODUCTS:
            print('PRODUCT={0}'.format(p), file=fid)
            #-- get directories from remote directory
            product_directory = '{0}.{1}'.format(p,RELEASE)
            PATH = [HOST,'ATLAS',product_directory]
            #-- compile regular expression operator
            args=(p,regex_track,regex_granule,RELEASE,regex_version,regex_suffix)
            R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)
            #-- read and parse request for subdirectories (find column names)
            remote_sub,collastmod = icesat2_toolkit.utilities.nsidc_list(PATH,
                build=False,timeout=120,parser=parser,pattern=R2,sort=True)
            #-- for each remote subdirectory
            for sd in remote_sub:
                #-- local directory for product and subdirectory
                if FLATTEN:
                    local_dir = os.path.expanduser(ddir)
                else:
                    local_dir = os.path.join(ddir,product_directory,sd)
                #-- check if data directory exists and recursively create if not
                os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
                #-- find ICESat-2 data files
                PATH = [HOST,'ATLAS',product_directory,sd]
                remote_dir = posixpath.join(HOST,'ATLAS',product_directory,sd)
                #-- find matching files (for granule, release, version, track)
                colnames,collastmod = icesat2_toolkit.utilities.nsidc_list(PATH,
                    build=False,timeout=120,parser=parser,pattern=R1,sort=True)
                #-- build lists of each ICESat-2 data file
                for colname,remote_mtime in zip(colnames,collastmod):
                    #-- remote and local versions of the file
                    remote_files.append(posixpath.join(remote_dir,colname))
                    local_files.append(os.path.join(local_dir,colname))
                    remote_mtimes.append(remote_mtime)
            #-- close request
            req = None

    #-- sync in series if PROCESSES = 0
    if (PROCESSES == 0):
        #-- sync each ICESat-2 data file
        for i,remote_file in enumerate(remote_files):
            #-- sync ICESat-2 files with NSIDC server
            output = http_pull_file(remote_file, remote_mtimes[i],
                local_files[i], LIST, CLOBBER, MODE)
            #-- print the output string
            print(output, file=fid)
    else:
        #-- sync in parallel with multiprocessing Pool
        pool = mp.Pool(processes=PROCESSES)
        #-- sync each ICESat-2 data file
        output = []
        for i,remote_file in enumerate(remote_files):
            #-- sync ICESat-2 files with NSIDC server
            args = (remote_file,remote_mtimes[i],local_files[i])
            kwds = dict(LIST=LIST, CLOBBER=CLOBBER, MODE=MODE)
            output.append(pool.apply_async(multiprocess_sync,
                args=args,kwds=kwds))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()
        #-- print the output string
        for out in output:
            print(out.get(), file=fid)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(ddir,LOGFILE), MODE)

#-- PURPOSE: wrapper for running the sync program in multiprocessing mode
def multiprocess_sync(remote_file, remote_mtime, local_file,
    LIST=False, CLOBBER=False, MODE=0o775):
    try:
        output = http_pull_file(remote_file,remote_mtime,local_file,
            LIST,CLOBBER,MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
    else:
        return output

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(remote_file,remote_mtime,local_file,LIST,CLOBBER,MODE):
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
            #-- Create and submit request. There are a wide range of exceptions
            #-- that can be thrown here, including HTTPError and URLError.
            request = urllib2.Request(remote_file)
            response = urllib2.urlopen(request)
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- copy contents to local file using chunked transfer encoding
            #-- transfer should work properly with ascii and binary data formats
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)
        #-- return the output string
        return output

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -U X, --user=X\t\tUsername for NASA Earthdata Login')
    print(' -N X, --netrc=X\tPath to .netrc file for authentication')
    print(' -D X, --directory=X\tWorking data directory')
    print(' -Y X, --year=X\t\tYears to sync separated by commas')
    print(' -S X, --subdirectory=X\tSubdirectories to sync separated by commas')
    print(' --release=X\t\tICESat-2 data release to sync')
    print(' --version=X\t\tICESat-2 data version to sync')
    print(' --granule=X\t\tICESat-2 granule regions to sync')
    print(' --track=X\t\tICESat-2 reference ground tracks to sync')
    print(' --auxiliary\t\tSync ICESat-2 auxiliary files for each HDF5 file')
    print(' -I X, --index=X\t\tInput index of ICESat-2 files to sync')
    print(' -F, --flatten\t\tDo not create subdirectories')
    print(' -P X, --np=X\t\tNumber of processes to use in file downloads')
    print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
    print(' -L, --list\t\tOnly print files that are to be transferred')
    print(' -C, --clobber\t\tOverwrite existing data in transfer')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'NSIDC_IceSat-2_sync_{0}.log'.format(today)
    print('    Log file format: {0}\n'.format(LOGFILE))

#-- Main program that calls nsidc_icesat2_sync()
def main():
    #-- Read the system arguments listed after the program
    short_options = 'hU:N:D:Y:S:I:FP:LCM:l'
    long_options=['help','user=','netrc=','directory=','year=','subdirectory=',
        'release=','version=','granule=','track=','auxiliary','index=',
        'flatten','np=','list','log','mode=','clobber']
    optlist,arglist=getopt.getopt(sys.argv[1:],short_options,long_options)

    #-- command line parameters
    USER = ''
    NETRC = None
    #-- Working data directory
    DIRECTORY = os.getcwd()
    YEARS = None
    SUBDIRECTORY = None
    VERSIONS = np.arange(1,10)
    RELEASE = '003'
    GRANULES = np.arange(1,15)
    TRACKS = np.arange(1,1388)
    AUXILIARY = False
    INDEX = None
    FLATTEN = False
    #-- sync in series if processes is 0
    PROCESSES = 0
    LIST = False
    LOG = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    CLOBBER = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-Y","--year"):
            YEARS = np.array(arg.split(','), dtype=np.int)
        elif opt in ("-S","--subdirectory"):
            SUBDIRECTORY = arg.split(',')
        elif opt in ("-U","--user"):
            USER = arg
        elif opt in ("-N","--netrc"):
            NETRC = os.path.expanduser(arg)
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--release",):
            RELEASE = '{0:03d}'.format(int(arg))
        elif opt in ("--version",):
            VERSIONS = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--granule",):
            GRANULES = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--track",):
            TRACKS = np.sort(arg.split(',')).astype(np.int)
        elif opt in ("--auxiliary",):
            AUXILIARY = True
        elif opt in ("-I","--index"):
            INDEX = os.path.expanduser(arg)
        elif opt in ("-P","--np"):
            PROCESSES = int(arg)
        elif opt in ("-F","--flatten"):
            FLATTEN = True
        elif opt in ("-L","--list"):
            LIST = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-C","--clobber"):
            CLOBBER = True

    #-- Pre-ICESat-2 and IceBridge Products
    PROD = {}
    PROD['ATL03'] = 'Global Geolocated Photon Data'
    PROD['ATL04'] = 'Normalized Relative Backscatter'
    PROD['ATL06'] = 'Land Ice Height'
    PROD['ATL07'] = 'Sea Ice Height'
    PROD['ATL08'] = 'Land and Vegetation Height'
    PROD['ATL09'] = 'Atmospheric Layer Characteristics'
    PROD['ATL10'] = 'Sea Ice Freeboard'
    PROD['ATL12'] = 'Ocean Surface Height'
    PROD['ATL13'] = 'Inland Water Surface Height'

    #-- enter dataset to transfer as system argument
    if not INDEX and not arglist:
        for key,val in PROD.items():
            print('{0}: {1}'.format(key, val))
        raise Exception('No System Arguments Listed')

    #-- check that each data product entered was correctly typed
    keys = ','.join(sorted([key for key in PROD.keys()]))
    for p in arglist:
        if p not in PROD.keys():
            raise IOError('Incorrect Data Product Entered ({0})'.format(keys))

    #-- NASA Earthdata hostname
    HOST = 'urs.earthdata.nasa.gov'
    #-- get authentication
    if not USER and not NETRC:
        #-- check that NASA Earthdata credentials were entered
        USER = builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))
    elif NETRC:
        USER,LOGIN,PASSWORD = netrc.netrc(NETRC).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))

    #-- build a urllib opener for NSIDC
    #-- Add the username and password for NASA Earthdata Login system
    icesat2_toolkit.utilities.build_opener(USER,PASSWORD)

    #-- check internet connection before attempting to run program
    #-- check NASA earthdata credentials before attempting to run program
    if icesat2_toolkit.utilities.check_credentials():
        nsidc_icesat2_sync(DIRECTORY, arglist, RELEASE, VERSIONS, GRANULES,
            TRACKS, YEARS=YEARS, SUBDIRECTORY=SUBDIRECTORY, AUXILIARY=AUXILIARY,
            INDEX=INDEX, FLATTEN=FLATTEN, PROCESSES=PROCESSES, LOG=LOG,
            LIST=LIST, MODE=MODE, CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
    main()
