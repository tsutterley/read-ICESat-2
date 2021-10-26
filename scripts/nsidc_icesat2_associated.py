#!/usr/bin/env python
u"""
nsidc_icesat2_associated.py
Written by Tyler Sutterley (10/2021)

Acquires ICESat-2 datafiles from the National Snow and Ice Data Center (NSIDC)
    server that is associated with an input file

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_associated.py --user <username> --release 003
        --product ATL03 list_of_ATL06_files
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for alternative authentication
    -D X, --directory: working data directory
    -p X, --product: associated product to download
        ATL03: Global Geolocated Photon Data
        ATL04: Normalized Relative Backscatter
        ATL06: Land Ice Height
        ATL07: Sea Ice Height
        ATL08: Land and Vegetation Height
        ATL09: Atmospheric Layer Characteristics
        ATL10: Sea Ice Freeboard
        ATL12: Ocean Surface Height
        ATL13: Inland Water Surface Height
    -a, --auxiliary: Download ICESat-2 auxiliary files for each HDF5 file
    -F, --flatten: Do not create subdirectories
    -P X, --np X: Number of processes to use in file downloads
    -T X, --timeout X: Timeout in seconds for blocking operations
    -R X, --retry X: Connection retry attempts
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
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: set context for multiprocessing to fork child processes
        added a file length check to validate downloaded files
    Updated 05/2021: added options for connection timeout and retry attempts
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 11/2020: nsidc_list will output a string for errors
    Updated 10/2020: using argparse to set parameters
        added multiprocessing option for parallel download
    Updated 08/2020: moved urllib opener to utilities. add credential check
    Updated 05/2020: added option netrc to use alternative authentication
        adjust regular expression to allow syncing of ATL07 sea ice products
        adjust regular expression for auxiliary products
    Updated 09/2019: added ssl context to urlopen headers
    Written 08/2019
"""
from __future__ import print_function

import sys
import os
import re
import netrc
import shutil
import getpass
import logging
import builtins
import argparse
import traceback
import posixpath
import lxml.etree
import multiprocessing as mp
import icesat2_toolkit.utilities

#-- PURPOSE: download the ICESat-2 elevation data from NSIDC matching an file
def nsidc_icesat2_associated(file_list, PRODUCT, DIRECTORY=None,
    AUXILIARY=False, FLATTEN=False, PROCESSES=0, TIMEOUT=None,
    RETRY=1, MODE=0o775):

    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()
    #-- logging to standard output
    logging.basicConfig(level=logging.INFO)

    #-- remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    #-- regular expression operator for extracting information from files
    rx = re.compile(r'(processed_)?(ATL\d{2})(-\d{2})?_(\d{4})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    #-- regular expression pattern for finding specific files
    regex_suffix = '(.*?)$' if AUXILIARY else '(h5)$'
    remote_regex_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_(\d{{2}})(.*?).{5}')

    #-- build list of remote files, remote modification times and local files
    original_files = []
    remote_files = []
    remote_mtimes = []
    local_files = []
    #-- for each input file
    for input_file in file_list:
        #-- extract parameters from ICESat-2 ATLAS HDF5 file name
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX = \
            rx.findall(input_file).pop()
        #-- get directories from remote directory
        product_directory = '{0}.{1}'.format(PRODUCT,RL)
        sd = '{0}.{1}.{2}'.format(YY,MM,DD)
        PATH = [HOST,'ATLAS',product_directory,sd]
        #-- local and remote data directories
        remote_dir=posixpath.join(*PATH)
        temp=os.path.dirname(input_file) if (DIRECTORY is None) else DIRECTORY
        local_dir=os.path.expanduser(temp) if FLATTEN else os.path.join(temp,sd)
        #-- create output directory if not currently existing
        if not os.access(local_dir, os.F_OK):
            os.makedirs(local_dir, MODE)
        #-- compile regular expression operator for file parameters
        args = (PRODUCT,TRK,CYC,GRN,RL,regex_suffix)
        R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)
        #-- find associated ICESat-2 data file
        #-- find matching files (for granule, release, version, track)
        colnames,collastmod,colerror=icesat2_toolkit.utilities.nsidc_list(PATH,
            build=False,
            timeout=TIMEOUT,
            parser=parser,
            pattern=R1,
            sort=True)
        #-- print if file was not found
        if not colnames:
            logging.critical(colerror)
            continue
        #-- add to lists
        for colname,remote_mtime in zip(colnames,collastmod):
            #-- save original file to list (expands if getting auxiliary files)
            original_files.append(input_file)
            #-- remote and local versions of the file
            remote_files.append(posixpath.join(remote_dir,colname))
            local_files.append(os.path.join(local_dir,colname))
            remote_mtimes.append(remote_mtime)

    #-- download in series if PROCESSES = 0
    if (PROCESSES == 0):
        #-- download each associated ICESat-2 data file
        for i,input_file in enumerate(original_files):
            #-- download associated ICESat-2 files with NSIDC server
            args = (remote_files[i],remote_mtimes[i],local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT,RETRY=RETRY,MODE=MODE)
            out = http_pull_file(*args,**kwds)
            #-- print the output string
            logging.info('{0}\n{1}'.format(input_file,out))
    else:
        #-- set multiprocessing start method
        ctx = mp.get_context("fork")
        #-- download in parallel with multiprocessing Pool
        pool = ctx.Pool(processes=PROCESSES)
        #-- download each associated ICESat-2 data file
        output = []
        for i,input_file in enumerate(original_files):
            #-- download associated ICESat-2 files with NSIDC server
            args = (remote_files[i],remote_mtimes[i],local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT,RETRY=RETRY,MODE=MODE)
            out=pool.apply_async(multiprocess_sync,args=args,kwds=kwds)
            output.append('{0}\n{1}'.format(input_file,out))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()
        #-- print the output string
        for out in output:
            logging.info(out.get())

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
def http_pull_file(remote_file,remote_mtime,local_file,
    TIMEOUT=None, RETRY=1, MODE=0o775):
    #-- Printing files transferred
    output = '{0} -->\n\t{1}\n'.format(remote_file,local_file)
    #-- chunked transfer encoding size
    CHUNK = 16 * 1024
    #-- attempt to download a file up to a set number of times
    retry_download(remote_file, LOCAL=local_file, TIMEOUT=TIMEOUT,
        RETRY=RETRY, CHUNK=CHUNK)
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
            #-- copy contents to file using chunked transfer encoding
            #-- transfer should work with ascii and binary data formats
            with open(LOCAL, 'wb') as f:
                shutil.copyfileobj(response, f, CHUNK)
            local_length = os.path.getsize(LOCAL)
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

#-- Main program that calls nsidc_icesat2_associated()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Program to acquire ICESat-2 datafiles from the NSIDC
            server that is associated with an input file
            """
    )
    #-- ICESat-2 Products
    PRODUCTS = {}
    PRODUCTS['ATL03'] = 'Global Geolocated Photon Data'
    PRODUCTS['ATL04'] = 'Normalized Relative Backscatter'
    PRODUCTS['ATL06'] = 'Land Ice Height'
    PRODUCTS['ATL07'] = 'Sea Ice Height'
    PRODUCTS['ATL08'] = 'Land and Vegetation Height'
    PRODUCTS['ATL09'] = 'Atmospheric Layer Characteristics'
    PRODUCTS['ATL10'] = 'Sea Ice Freeboard'
    PRODUCTS['ATL12'] = 'Ocean Surface Height'
    PRODUCTS['ATL13'] = 'Inland Water Surface Height'
    #-- command line parameters
    parser.add_argument('file',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='ICESat-2 products to associate')
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
    #-- ICESat-2 parameters
    #-- ICESat-2 data product
    parser.add_argument('--product','-p',
        metavar='PRODUCTS', type=str,
        choices=PRODUCTS.keys(), default='ATL06',
        help='Associated ICESat-2 data product to download')
    #-- download auxiliary files
    parser.add_argument('--auxiliary','-a',
        default=False, action='store_true',
        help='Sync ICESat-2 auxiliary files for each HDF5 file')
    #-- output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
    #-- run download in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in downloading files')
    #-- connection timeout and number of retry attempts
    parser.add_argument('--timeout','-T',
        type=int, default=120,
        help='Timeout in seconds for blocking operations')
    parser.add_argument('--retry','-R',
        type=int, default=5,
        help='Connection retry attempts')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args,_ = parser.parse_known_args()

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
        nsidc_icesat2_associated(args.file, args.product,
            DIRECTORY=args.directory, AUXILIARY=args.auxiliary,
            FLATTEN=args.flatten, PROCESSES=args.np, TIMEOUT=args.timeout,
            RETRY=args.retry, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
