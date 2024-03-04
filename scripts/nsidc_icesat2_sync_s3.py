#!/usr/bin/env python
u"""
nsidc_icesat2_sync_s3.py
Written by Tyler Sutterley (03/2024)

Acquires ICESat-2 datafiles from the National Snow and Ice Data Center (NSIDC)
    and transfers to an AWS S3 bucket using a local machine as pass through

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_sync_s3.py --user=<username> --release=001 ATL06
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
    -W X, --password X: Password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for alternative authentication
    --aws-access-key-id X: AWS Access Key ID
    --aws-secret-access-key X: AWS Secret Key
    --aws-region-name X: AWS Region Name
    --s3-bucket-name X: AWS S3 Bucket Name
    --s3-bucket-path X: AWS S3 Bucket Path
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
    -C, --clobber: Overwrite existing data in transfer

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    boto3: Amazon Web Services (AWS) SDK for Python
        https://boto3.amazonaws.com/v1/documentation/api/latest/index.html

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 03/2024: use pathlib to define and operate on paths
    Updated 09/2023: generalized regular expressions for non-entered cases
    Updated 12/2022: single implicit import of altimetry tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 03/2022: use attempt login function to check credentials
    Updated 02/2022: added option to sync specific orbital cycles
    Updated 10/2021: using python logging for handling verbose output
    Forked 08/2021 from nsidc_icesat2_sync.py
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
import logging
import pathlib
import argparse
import warnings
import posixpath
import traceback
import lxml.etree
import multiprocessing as mp
import icesat2_toolkit as is2tk

# attempt imports
try:
    import boto3
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    warnings.warn("boto3 not available", ImportWarning)

# PURPOSE: sync the ICESat-2 elevation data from NSIDC
def nsidc_icesat2_sync_s3(aws_access_key_id, aws_secret_access_key,
    aws_region_name, s3_bucket_name, s3_bucket_path,
    PRODUCTS, RELEASE, VERSIONS, GRANULES, TRACKS,
    YEARS=None, SUBDIRECTORY=None, CYCLES=None, REGION=None,
    AUXILIARY=False, INDEX=None, FLATTEN=False, TIMEOUT=None,
    RETRY=1, PROCESSES=0, CLOBBER=False):

    # get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    # get s3 object and bucket object
    s3 = session.resource('s3')
    bucket = s3.Bucket(s3_bucket_name)

    # logging to standard output
    logging.basicConfig(level=logging.INFO)

    # compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()

    # remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    # regular expression operator for finding files of a particular granule
    # find ICESat-2 HDF5 files in the subdirectory for product and release
    if TRACKS:
        regex_track = r'|'.join([rf'{T:04d}' for T in TRACKS])
    else:
        regex_track = r'\d{4}'
    if CYCLES:
        regex_cycle = r'|'.join([rf'{C:02d}' for C in CYCLES])
    else:
        regex_cycle = r'\d{2}'
    if GRANULES:
        regex_granule = r'|'.join([rf'{G:02d}' for G in GRANULES])
    else:
        regex_granule = r'\d{2}'
    if VERSIONS:
        regex_version = r'|'.join([rf'{V:02d}' for V in VERSIONS])
    else:
        regex_version = r'\d{2}'
    regex_suffix = r'(.*?)' if AUXILIARY else r'(h5|nc)'
    default_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).{6}$')
    ATL11_pattern = r'({0})_({1})({2})_(\d{{2}})(\d{{2}})_({3})_({4})(.*?).{5}$'
    ATL1415_pattern = r'({0})_({1})_(\d{{2}})(\d{{2}})_({3})_({4})(.*?).{5}$'

    # regular expression operator for finding subdirectories
    if SUBDIRECTORY:
        # Sync particular subdirectories for product
        R2 = re.compile(r'('+r'|'.join(SUBDIRECTORY)+r')', re.VERBOSE)
    elif YEARS:
        # Sync particular years for product
        regex_pattern = r'|'.join(rf'{y:d}' for y in YEARS)
        R2 = re.compile(rf'({regex_pattern}).(\d+).(\d+)', re.VERBOSE)
    else:
        # Sync all available subdirectories for product
        R2 = re.compile(r'(\d+).(\d+).(\d+)', re.VERBOSE)

    # build list of remote files, remote modification times and AWS S3 files
    remote_files = []
    remote_mtimes = []
    s3_files = []
    # build lists of files or use existing index file
    if INDEX:
        # read the index file, split at lines and remove all commented lines
        INDEX = pathlib.Path(INDEX).expanduser().absolute()
        with INDEX.open(mode='r', encoding='utf8') as f:
            files = [i for i in f.read().splitlines() if re.match(r'^(?!\#)',i)]
        # regular expression operator for extracting information from files
        rx = re.compile(r'(ATL\d{2})(-\d{2})?_(\d{4})(\d{2})(\d{2})(\d{2})'
            r'(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
        # for each line in the index
        for f in files:
            # extract parameters from ICESat-2 ATLAS HDF5 file
            PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX = \
                rx.findall(f).pop()
            # get directories from remote directory
            product_directory = f'{PRD}.{RL}'
            sd = f'{YY}.{MM}.{DD}'
            PATH = [HOST,'ATLAS',product_directory,sd]
            remote_dir = posixpath.join(HOST,'ATLAS',product_directory,sd)
            # AWS S3 directory for product and subdirectory
            if FLATTEN:
                s3_path = posixpath.expanduser(s3_bucket_path)
            else:
                s3_path = posixpath.join(s3_bucket_path,product_directory,sd)
            # find ICESat-2 data file to get last modified time
            # find matching files (for granule, release, version, track)
            names,lastmod,error = is2tk.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=f.strip()
            )
            # print if file was not found
            if not names:
                logging.critical(error)
                continue
            # add to lists
            for colname,remote_mtime in zip(names,lastmod):
                # remote and AWS S3 versions of the file
                remote_files.append(posixpath.join(remote_dir,colname))
                s3_files.append(posixpath.join(s3_path,colname))
                remote_mtimes.append(remote_mtime)
    else:
        # for each ICESat-2 product listed
        for p in PRODUCTS:
            logging.info(f'PRODUCT={p}')
            # get directories from remote directory
            product_directory = f'{p}.{RELEASE}'
            PATH = [HOST,'ATLAS',product_directory]
            # compile regular expression operator
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
            # read and parse request for subdirectories (find column names)
            remote_sub,_,error = is2tk.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=R2,
                sort=True
            )
            # print if subdirectory was not found
            if not remote_sub:
                logging.critical(error)
                continue
            # for each remote subdirectory
            for sd in remote_sub:
                # AWS S3 directory for product and subdirectory
                if FLATTEN:
                    s3_path = posixpath.expanduser(s3_bucket_path)
                else:
                    s3_path = posixpath.join(s3_bucket_path,product_directory,sd)
                logging.info(f"Building file list: {sd}")
                # find ICESat-2 data files
                PATH = [HOST,'ATLAS',product_directory,sd]
                remote_dir = posixpath.join(HOST,'ATLAS',product_directory,sd)
                # find matching files (for granule, release, version, track)
                names,lastmod,error = is2tk.utilities.nsidc_list(PATH,
                    build=False,
                    timeout=TIMEOUT,
                    parser=parser,
                    pattern=R1,
                    sort=True
                )
                # print if file was not found
                if not names:
                    logging.critical(error)
                    continue
                # build lists of each ICESat-2 data file
                for colname,remote_mtime in zip(names,lastmod):
                    # remote and AWS S3 versions of the file
                    remote_files.append(posixpath.join(remote_dir,colname))
                    s3_files.append(posixpath.join(s3_path,colname))
                    remote_mtimes.append(remote_mtime)

    # sync in series if PROCESSES = 0
    if (PROCESSES == 0):
        # sync each ICESat-2 data file
        for i,remote_file in enumerate(remote_files):
            # sync ICESat-2 files with NSIDC server
            args = (bucket,remote_file,remote_mtimes[i],s3_files[i])
            kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY, CLOBBER=CLOBBER)
            output = http_pull_file(*args, **kwds)
            # print the output string
            logging.info(output) if output else None
    else:
        # set multiprocessing start method
        ctx = mp.get_context("fork")
        # sync in parallel with multiprocessing Pool
        pool = ctx.Pool(processes=PROCESSES)
        # sync each ICESat-2 data file
        out = []
        for i,remote_file in enumerate(remote_files):
            # sync ICESat-2 files with NSIDC server
            args = (bucket,remote_file,remote_mtimes[i],s3_files[i])
            kwds = dict(TIMEOUT=TIMEOUT, RETRY=RETRY, CLOBBER=CLOBBER)
            out.append(pool.apply_async(multiprocess_sync,
                args=args,kwds=kwds))
        # start multiprocessing jobs
        # close the pool
        # prevents more tasks from being submitted to the pool
        pool.close()
        # exit the completed processes
        pool.join()
        # print the output string
        for output in out:
            temp = output.get()
            logging.info(temp) if temp else None

# PURPOSE: wrapper for running the sync program in multiprocessing mode
def multiprocess_sync(*args, **kwds):
    try:
        output = http_pull_file(*args, **kwds)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
    else:
        return output

# PURPOSE: pull file from a remote host checking if file exists on S3
# and if the remote file is newer than the AWS S3 bucket file
def http_pull_file(bucket, remote_file, remote_mtime, s3_file,
    TIMEOUT=None, RETRY=1, CLOBBER=False):
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if s3 bucket version of file exists
    try:
        # check last modification time of s3 file
        obj = bucket.Object(key=s3_file)
        s3_mtime = obj.last_modified
    except:
        TEST = True
        OVERWRITE = ' (new)'
    else:
        # if remote file is newer: overwrite the s3 bucket file
        if (remote_mtime > s3_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    # if file does not exist on S3, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        # output string for printing files transferred
        output = f'{remote_file} -->\n\t{s3_file}{OVERWRITE}\n'
        # copy bytes to s3 object
        retry_download(remote_file, BUCKET=bucket, LOCAL=s3_file,
            TIMEOUT=TIMEOUT, RETRY=RETRY)
        # return the output string
        return output

# PURPOSE: Try downloading a file up to a set number of times
def retry_download(remote_file, BUCKET=None, LOCAL=None, TIMEOUT=None, RETRY=1):
    # attempt to download up to the number of retries
    retry_counter = 0
    while (retry_counter < RETRY):
        # attempt to retrieve file from https server
        try:
            # Create and submit request.
            # There are a range of exceptions that can be thrown here
            # including HTTPError and URLError.
            request = is2tk.utilities.urllib2.Request(remote_file)
            response = is2tk.utilities.urllib2.urlopen(request,
                timeout=TIMEOUT)
            # get the length of the remote file
            remote_length = int(response.headers['content-length'])
            # upload file response to s3 bucket
            BUCKET.upload_fileobj(response, LOCAL)
            obj = BUCKET.Object(key=LOCAL)
            local_length = obj.content_length
        except:
            pass
        else:
            # check that downloaded file matches original length
            if (local_length == remote_length):
                break
        # add to retry counter
        retry_counter += 1
    # check if maximum number of retries were reached
    if (retry_counter == RETRY):
        raise TimeoutError('Maximum number of retries reached')

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Acquires ICESat-2 datafiles from the National Snow
            and Ice Data Center (NSIDC) and transfers to an AWS S3 bucket
            using a local machine as pass through
            """
    )
    # command line parameters
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('products',
        metavar='PRODUCTS', type=str, nargs='*', default=[],
        help='ICESat-2 products to sync')
    # NASA Earthdata credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('EARTHDATA_USERNAME'),
        help='Username for NASA Earthdata Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('EARTHDATA_PASSWORD'),
        help='Password for NASA Earthdata Login')
    parser.add_argument('--netrc','-N',
        type=pathlib.Path,
        default=pathlib.Path.home().joinpath('.netrc'),
        help='Path to .netrc file for authentication')
    # AWS credentials, bucket and path
    parser.add_argument('--aws-access-key-id',
        type=str, required=True,
        help='AWS Access Key ID')
    parser.add_argument('--aws-secret-access-key',
        type=str, required=True,
        help='AWS Secret Key')
    parser.add_argument('--aws-region-name',
        type=str, default='us-west-2',
        help='AWS Region Name')
    parser.add_argument('--s3-bucket-name',
        type=str, required=True,
        help='AWS S3 Bucket Name')
    parser.add_argument('--s3-bucket-path',
        type=str, default='',
        help='AWS S3 Bucket Path')
    # years of data to sync
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to sync')
    # subdirectories of data to sync
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to sync')
    # ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='006',
        help='ICESat-2 Data Release')
    # ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+',
        help='ICESat-2 Data Version')
    # ICESat-2 granule region
    region = parser.add_mutually_exclusive_group(required=False)
    region.add_argument('--granule','-g',
        metavar='GRANULE', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
    # ICESat-2 orbital cycle
    parser.add_argument('--cycle','-c',
        type=int, nargs='+', default=None,
        help='ICESat-2 orbital cycles to sync')
    # ICESat-2 ATL14 and 15 named regions
    ATL1415_regions = ['AA','AK','CN','CS','GL','IC','SV','RU']
    region.add_argument('--region','-n',
        metavar='REGION', type=str, nargs='+',
        choices=ATL1415_regions, default=['AA','GL'],
        help='ICESat-2 Named Region (ATL14/ATL15)')
    # ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs)')
    # sync auxiliary files
    parser.add_argument('--auxiliary','-a',
        default=False, action='store_true',
        help='Sync ICESat-2 auxiliary files for each HDF5 file')
    # sync using files from an index
    group.add_argument('--index','-i',
        type=pathlib.Path,
        help='Input index of ICESat-2 files to sync')
    # output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
    # run sync in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file downloads')
    # connection timeout and number of retry attempts
    parser.add_argument('--timeout','-T',
        type=int, default=120,
        help='Timeout in seconds for blocking operations')
    parser.add_argument('--retry','-R',
        type=int, default=5,
        help='Connection retry attempts')
    # clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # NASA Earthdata hostname
    HOST = 'urs.earthdata.nasa.gov'
    # build a urllib opener for NASA Earthdata
    # check internet connection before attempting to run program
    opener = is2tk.utilities.attempt_login(HOST,
        username=args.user, password=args.password,
        netrc=args.netrc)

    # check NASA earthdata credentials before attempting to run program
    nsidc_icesat2_sync_s3(args.aws_access_key_id,
        args.aws_secret_access_key, args.aws_region_name,
        args.s3_bucket_name, args.s3_bucket_path,
        args.products, args.release, args.version, args.granule,
        args.track, YEARS=args.year, SUBDIRECTORY=args.subdirectory,
        CYCLES=args.cycle, REGION=args.region, AUXILIARY=args.auxiliary,
        INDEX=args.index, FLATTEN=args.flatten, PROCESSES=args.np,
        TIMEOUT=args.timeout, RETRY=args.retry, CLOBBER=args.clobber)

# run main program
if __name__ == '__main__':
    main()
