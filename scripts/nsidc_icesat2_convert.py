#!/usr/bin/env python
u"""
nsidc_icesat2_convert.py
Written by Tyler Sutterley (03/2024)

Acquires ICESat-2 datafiles from NSIDC and directly converts to
    zarr datafiles or rechunked HDF5 files

zarr files make large datasets easily accessible to distributed computing on
    both local filesystems and cloud-based object stores
    - arrays are divided into chunks and compressed
    - metadata are stored in lightweight .json files

rechunked HDF5 files can be more optimized for cloud-based object stores

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_convert.py --user=<username> --release=003 ATL06
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
    -D X, --directory X: working data directory
    -Y X, --year X: years to sync
    -S X, --subdirectory X: specific subdirectories to sync
    -r X, --release X: ICESat-2 data release to sync
    -v X, --version X: ICESat-2 data version to sync
    -t X, --track X: ICESat-2 reference ground tracks to sync
    -c X, --cycle=X: ICESat-2 cycles to sync
    -g X, --granule X: ICESat-2 granule regions to sync
    -f X, --format X: output file format (zarr, HDF5)
    --chunks X: Rechunk output files to size
    -a X, --auxiliary: Sync ICESat-2 auxiliary files for each HDF5 file
    -I X, --index X: Input index of ICESat-2 files to sync
    -F, --flatten: Do not create subdirectories
    -P X, --np X: Number of processes to use in file downloads
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
    zarr: Chunked, compressed, N-dimensional arrays in Python
        https://github.com/zarr-developers/zarr-python
        https://zarr.readthedocs.io/en/stable/index.html
    pandas: Python Data Analysis Library
        https://pandas.pydata.org/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

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
    Updated 07/2021: set context for multiprocessing to fork child processes
    Updated 05/2021: added options for connection timeout and retry attempts
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
        use regex backslash for comment special characters
    Updated 01/2021: generalized to output either zarr or rechunked HDF5
    Updated 11/2020: nsidc_list will output a string for errors
    Updated 10/2020: using argparse to set parameters
        added chunks keyword to rechunk output zarr files
        using convert module to convert from HDF5 to zarr
    Updated 09/2020: use urllib imported in utilities
    Updated 08/2020: moved urllib opener to utilities. add credential check
        moved urllib directory listing to utilities
    Updated 07/2020: added option index to use a list of files to sync
    Updated 06/2020: transfer to BytesIO obj in chunks to prevent overflow error
    Written 06/2020
"""
from __future__ import print_function

import sys
import os
import re
import shutil
import logging
import pathlib
import argparse
import posixpath
import traceback
import lxml.etree
import calendar, time
import multiprocessing as mp
import icesat2_toolkit as is2tk

# PURPOSE: sync the ICESat-2 elevation data from NSIDC and convert to
# either zarr or rechunked HDF5
def nsidc_icesat2_convert(DIRECTORY, PRODUCTS, RELEASE, VERSIONS, GRANULES, TRACKS,
    YEARS=None, SUBDIRECTORY=None, CYCLES=None, FORMAT=None, CHUNKS=None,
    AUXILIARY=False, INDEX=None, FLATTEN=False, TIMEOUT=None, RETRY=1,
    LOG=False, LIST=False, PROCESSES=0, CLOBBER=False, MODE=None):

    # check if directory exists and recursively create if not
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output of synchronized files
    if LOG:
        # format: NSIDC_ICESat-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = DIRECTORY.joinpath(f'NSIDC_ICESat-2_sync_{today}.log')
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'ICESat-2 Data Convert Log ({today})')
    else:
        # standard output (terminal output)
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
    regex_suffix = r'(.*?)' if AUXILIARY else r'(h5)'
    remote_regex_pattern=(r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).{6}$')

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

    # build list of remote files, remote modification times and local files
    remote_files = []
    remote_mtimes = []
    local_files = []
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
            PATH = [HOST, 'ATLAS', product_directory, sd]
            remote_dir = posixpath.join(*PATH)
            # local directory
            local_dir = pathlib.Path(DIRECTORY).expanduser().absolute()
            if not FLATTEN:
                local_dir = local_dir.joinpath(product_directory,sd)
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
                # remote and local versions of the file
                remote_files.append(posixpath.join(remote_dir,colname))
                local_files.append(local_dir.joinpath(colname))
                remote_mtimes.append(remote_mtime)
    else:
        # for each ICESat-2 product listed
        for p in PRODUCTS:
            logging.info(f'PRODUCT={p}')
            # get directories from remote directory
            product_directory = f'{p}.{RELEASE}'
            PATH = [HOST, 'ATLAS', product_directory]
            # compile regular expression operator
            args = (p, regex_track, regex_cycle, regex_granule, RELEASE,
                regex_version, regex_suffix)
            R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)
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
                # local directory
                local_dir = pathlib.Path(DIRECTORY).expanduser().absolute()
                if not FLATTEN:
                    local_dir = local_dir.joinpath(product_directory,sd)
                # find ICESat-2 data files
                PATH = [HOST, 'ATLAS', product_directory, sd]
                remote_dir = posixpath.join(*PATH)
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
                    # remote and local versions of the file
                    remote_files.append(posixpath.join(remote_dir,colname))
                    local_files.append(local_dir.joinpath(colname))
                    remote_mtimes.append(remote_mtime)

    # sync in series if PROCESSES = 0
    if (PROCESSES == 0):
        # sync each ICESat-2 data file
        for i,remote_file in enumerate(remote_files):
            # sync ICESat-2 files with NSIDC server
            args = (remote_file, remote_mtimes[i], local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT,
                RETRY=RETRY,
                FORMAT=FORMAT,
                CHUNKS=CHUNKS,
                LIST=LIST,
                CLOBBER=CLOBBER,
                MODE=MODE
            )
            output = http_pull_file(*args,**kwds)
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
            args = (remote_file,remote_mtimes[i],local_files[i])
            kwds = dict(TIMEOUT=TIMEOUT,
                FORMAT=FORMAT,
                CHUNKS=CHUNKS,
                RETRY=RETRY,
                LIST=LIST,
                CLOBBER=CLOBBER,
                MODE=MODE
            )
            out.append(pool.apply_async(multiprocess_sync,args=args,kwds=kwds))
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

    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

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

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
        TIMEOUT=None,
        RETRY=1,
        FORMAT=None,
        CHUNKS=None,
        LIST=False,
        CLOBBER=False,
        MODE=0o775
    ):
    # check if data directory exists and recursively create if not
    local_file = pathlib.Path(local_file).expanduser().absolute()
    local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
    # convert HDF5 file from server into a new file format
    if (local_file.suffix == '.h5') and (FORMAT == 'zarr'):
        output_file = local_file.with_suffix('.zarr')
    elif (local_file.suffix == '.h5') and (FORMAT == 'HDF5'):
        output_file = local_file.with_suffix('.h5')
    else:
        output_file = local_file
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if output_file.exists():
        # check last modification time of local file
        local_mtime = output_file.stat().st_mtime
        # if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    # if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        # output string for printing files transferred
        output = f'{remote_file} -->\n\t{str(output_file)}{OVERWRITE}\n'
        # if executing copy command (not only printing the files)
        if not LIST:
            # chunked transfer encoding size
            CHUNK = 16 * 1024
            # attempt to download up to the number of retries
            retry_counter = 0
            while (retry_counter < RETRY):
                # attempt to retrieve file from https server
                try:
                    # Create and submit request
                    # There are a range of exceptions that can be thrown
                    # including HTTPError and URLError.
                    fid = is2tk.utilities.from_http(remote_file,
                        timeout=TIMEOUT,
                        context=None,
                        chunk=CHUNK
                    )
                except:
                    pass
                else:
                    break
                # add to retry counter
                retry_counter += 1
            # check if maximum number of retries were reached
            if (retry_counter == RETRY):
                raise TimeoutError('Maximum number of retries reached')
            # rewind retrieved binary to start of file
            fid.seek(0)
            if (local_file.suffix == '.h5'):
                # copy local HDF5 filename to BytesIO object attribute
                fid.filename = local_file
                # copy everything from the HDF5 file to the output file
                conv = is2tk.convert(filename=fid, reformat=FORMAT)
                conv.file_converter(chunks=CHUNKS)
            else:
                # copy contents to file using chunked transfer encoding
                # transfer should work with ascii and binary data formats
                with output_file.open('wb') as f:
                    shutil.copyfileobj(fid, f, CHUNK)
            # keep remote modification time of file and local access time
            os.utime(output_file, (output_file.stat().st_atime, remote_mtime))
            output_file.chmod(mode=MODE)
        # return the output string
        return output

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Acquires ICESat-2 datafiles from NSIDC and directly
            converts to zarr datafiles or rechunked HDF5 files
            """
    )
    # command line parameters
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('products',
        metavar='PRODUCTS', type=str, nargs='*', default=[],
        help='ICESat-2 products to convert')
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
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Working data directory')
    # years of data to run
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to run')
    # subdirectories of data to run
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to run')
    # ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='006',
        help='ICESat-2 Data Release')
    # ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+',
        help='ICESat-2 Data Version')
    # ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
    # ICESat-2 orbital cycle
    parser.add_argument('--cycle','-c',
        type=int, nargs='+', default=None,
        help='ICESat-2 orbital cycles to sync')
    # ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs)')
    # output file format
    parser.add_argument('--format','-f',
        type=str, choices=('zarr','HDF5'), default='zarr',
        help='Output file format')
    # rechunk output data
    parser.add_argument('--chunks',
        type=int,
        help='Rechunk output files to size')
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
    # run conversion in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file conversion')
    # connection timeout and number of retry attempts
    parser.add_argument('--timeout','-T',
        type=int, default=120,
        help='Timeout in seconds for blocking operations')
    parser.add_argument('--retry','-R',
        type=int, default=5,
        help='Connection retry attempts')
    # Output log file in form
    # NSIDC_IceSat-2_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    # clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    # permissions mode of the converted files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
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
    nsidc_icesat2_convert(args.directory, args.products, args.release,
        args.version, args.granule, args.track, YEARS=args.year,
        SUBDIRECTORY=args.subdirectory, CYCLES=args.cycle,
        FORMAT=args.format, CHUNKS=args.chunks, AUXILIARY=args.auxiliary,
        INDEX=args.index, FLATTEN=args.flatten, PROCESSES=args.np,
        TIMEOUT=args.timeout, RETRY=args.retry, LOG=args.log,
        LIST=args.list, CLOBBER=args.clobber, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
