#!/usr/bin/env python
u"""
nsidc_icesat2_dragann.py
Written by Tyler Sutterley (09/2023)

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
    -W X, --password X: Password for NASA Earthdata Login
    -N X, --netrc X: path to .netrc file for alternative authentication
    -D X, --directory X: working data directory
    -Y X, --year X: years to sync
    -S X, --subdirectory X: specific subdirectories to sync
    -r X, --release X: ICESat-2 data release to sync
    -v X, --version X: ICESat-2 data version to sync
    -t X, --track X: ICESat-2 reference ground tracks to sync
    -g X, --granule X: ICESat-2 granule regions to sync
    -c X, --cycle=X: ICESat-2 cycles to sync
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

UPDATE HISTORY:
    Updated 09/2023: generalized regular expressions for non-entered cases
    Updated 12/2022: single implicit import of altimetry tools
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 03/2022: use attempt login function to check credentials
    Updated 02/2022: added option to sync specific orbital cycles
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
import shutil
import logging
import pathlib
import argparse
import warnings
import posixpath
import traceback
import lxml.etree
import numpy as np
import calendar, time
import multiprocessing as mp
import icesat2_toolkit as is2tk

# attempt imports
try:
    import h5py
except ModuleNotFoundError:
    warnings.warn("h5py not available", ImportWarning)

# PURPOSE: sync ATL03 geolocated photon height products and appends the
# ATL08 DRAGANN classifications from NSIDC
def nsidc_icesat2_dragann(DIRECTORY, RELEASE, VERSIONS, GRANULES, TRACKS,
    YEARS=None, SUBDIRECTORY=None, CYCLES=None, FLATTEN=False,
    TIMEOUT=None, RETRY=1, LOG=False, LIST=False, CLOBBER=False, MODE=None):

    # check if directory exists and recursively create if not
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    DIRECTORY.mkdir(mode=MODE, parents=True, exist_ok=True)

    # output of synchronized files
    if LOG:
        # format: NSIDC_ICESat-2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = DIRECTORY.joinpath(f'NSIDC_ICESat-2_sync_{today}.log')
        logging.basicConfig(filename=LOGFILE, level=logging.INFO)
        logging.info(f'ICESat-2 DRAGANN Sync Log ({today})')
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
    regex_suffix = r'(h5)'
    remote_regex_pattern=(r'({0})_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
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

    # get directories from remote directory
    atl03_directory = f'ATL03.{RELEASE}'
    atl08_directory = f'ATL08.{RELEASE}'
    PATH = [HOST, 'ATLAS', atl08_directory]
    # compile regular expression operator
    args = ('ATL08', regex_track, regex_cycle, regex_granule,
        RELEASE, regex_version, regex_suffix)
    R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)
    # read and parse request for subdirectories (find column names)
    remote_sub,_,error = is2tk.utilities.nsidc_list(PATH,
        build=False,
        timeout=TIMEOUT,
        parser=parser,
        pattern=R2,
        sort=True
    )
    # for each remote subdirectory
    for sd in remote_sub:
        # local directory
        local_dir = pathlib.Path(DIRECTORY).expanduser().absolute()
        if not FLATTEN:
            local_dir = local_dir.joinpath(atl03_directory,sd)
        # find ICESat-2 data files
        PATH = [HOST, 'ATLAS', atl08_directory, sd]
        # find matching files (for granule, release, version, track)
        atl08s,lastmod,error = is2tk.utilities.nsidc_list(PATH,
            build=False,
            timeout=TIMEOUT,
            parser=parser,
            pattern=R1,
            sort=True
        )
        # print if file was not found
        if not atl08s:
            logging.critical(error)
            continue
        # build lists of each ICESat-2 data file
        for i,atl08 in enumerate(atl08s):
            # extract base parameters
            PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX,SFX = \
                R1.findall(atl08).pop()
            # compile regular expression operator
            args = ('ATL03',TRK,CYCL,GRAN,RL,regex_version,regex_suffix)
            R3 = re.compile(remote_regex_pattern.format(*args))
            PATH = [HOST,'ATLAS',atl03_directory,sd]
            # find associated ATL03 files
            atl03s,lastmod,error = is2tk.utilities.nsidc_list(PATH,
                build=False,
                timeout=TIMEOUT,
                parser=parser,
                pattern=R3,
                sort=True)
            # remote and local versions of the file
            for atl03,remote_mtime in zip(atl03s,lastmod):
                # sync ICESat-2 ATL03 files with NSIDC server
                remote_dir = posixpath.join(HOST, 'ATLAS', atl03_directory, sd)
                remote_file = posixpath.join(remote_dir, atl03)
                local_file = local_dir.joinpath(atl03)
                # recursively create data directory if not existing
                local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
                # download ATL03 file
                args = (remote_file, remote_mtime, local_file)
                kwds = dict(TIMEOUT=TIMEOUT,
                    RETRY=RETRY,
                    LIST=LIST,
                    CLOBBER=CLOBBER,
                    MODE=MODE
                )
                out = http_pull_file(*args, **kwds)
                logging.info(out) if out else None
                # append ATL08 dragann classifications
                PATH = [HOST, 'ATLAS', atl08_directory, sd, atl08]
                logging.info(posixpath.join(*PATH))
                remote_buffer,_ = is2tk.utilities.from_nsidc(PATH,
                    build=False,
                    timeout=TIMEOUT
                )
                # for each beam in the ATL03 file
                for gtx in is2tk.io.ATL03.find_beams(local_file):
                    # open ATL03 file in append mode
                    fileID = h5py.File(local_file, 'a')
                    # check if DRAGANN variables are already appended
                    if 'd_flag' in fileID[gtx]['heights'].keys():
                        # close the file and continue
                        fileID.close()
                        continue
                    # ATL03 20 meter segment id
                    segment_id=fileID[gtx]['geolocation']['segment_id'][:]
                    # first photon in the segment (0-based indexing)
                    ph_index_beg=fileID[gtx]['geolocation']['ph_index_beg'][:]-1
                    # photon event level delta time
                    delta_time = fileID[gtx]['heights']['delta_time']
                    # number of photon events in the beam group
                    n_pe, = delta_time.shape
                    # extract dragann classifiers for beam
                    mds,attrs = extract_dragann_classification(remote_buffer,
                        gtx, segment_id, ph_index_beg, n_pe)
                    for k,v in mds.items():
                        # create HDF5 variable
                        val = posixpath.join(gtx, 'heights', k)
                        h5 = fileID.create_dataset(val, np.shape(v),
                            data=v, dtype=v.dtype, compression='gzip')
                        h5.dims[0].attach_scale(delta_time)
                        # add HDF5 variable attributes
                        for att_name,att_val in attrs[k].items():
                            h5.attrs[att_name] = att_val
                    # close the ATL03 file
                    fileID.close()
                # keep remote modification time of file and local access time
                os.utime(local_file, (local_file.stat().st_atime, remote_mtime))
                # change the permissions mode
                local_file.chmod(mode=MODE)

    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

# PURPOSE: pull file from a remote host checking if file exists locally
# and if the remote file is newer than the local file
def http_pull_file(remote_file, remote_mtime, local_file,
    TIMEOUT=None, RETRY=1, LIST=False, CLOBBER=False, MODE=0o775):
    # check if data directory exists and recursively create if not
    local_file = pathlib.Path(local_file).expanduser().absolute()
    local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
    # if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # check if local version of file exists
    if local_file.exists():
        # check last modification time of local file
        local_mtime = local_file.stat().st_mtime
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
        output = f'{remote_file} -->\n\t{str(local_file)}{OVERWRITE}\n'
        # if executing copy command (not only printing the files)
        if not LIST:
            # chunked transfer encoding size
            CHUNK = 16 * 1024
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
                    # copy contents to file using chunked transfer encoding
                    # transfer should work with ascii and binary data formats
                    with local_file.open('wb') as f:
                        shutil.copyfileobj(response, f, CHUNK)
                except:
                    pass
                else:
                    break
                # add to retry counter
                retry_counter += 1
            # check if maximum number of retries were reached
            if (retry_counter == RETRY):
                raise TimeoutError('Maximum number of retries reached')
            # keep remote modification time of file and local access time
            os.utime(local_file, (local_file.stat().st_atime, remote_mtime))
            # change the permissions mode
            local_file.chmod(mode=MODE)
        # return the output string
        return output

# PURPOSE: Extract photon event classification from ATL08
def extract_dragann_classification(buffer, gtx, segment_id, ph_index_beg, n_pe):
    # expand path to HDF5 file
    if isinstance(buffer, (str, pathlib.Path)):
        buffer = pathlib.Path(buffer).expanduser().absolute()
    # open the HDF5 file for reading
    fileID = h5py.File(buffer, 'r')
    # allocate for output classifications
    output = {}
    # default is noise classification for all photon events
    output['d_flag'] = np.zeros((n_pe),dtype=np.int8)
    output['classed_pc_flag'] = np.zeros((n_pe),dtype=np.int8)
    attrs = dict(d_flag={},classed_pc_flag={})
    # attributes for the output remapped variables
    # DRAGANN classification
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
    # ATL08 photon event classification
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
    # check if subsetted beam contains photon event classifications
    try:
        fileID[gtx]['signal_photons']['delta_time']
    except KeyError:
        pass
    else:
        # segment id of the ATL08 classified photons
        ph_segment_id = fileID[gtx]['signal_photons']['ph_segment_id'][:]
        # index of the ATL08 classified photons within the ATL03 20m segments
        classed_pc_indx = fileID[gtx]['signal_photons']['classed_pc_indx'][:]-1
        # dragann flag and land classifications
        d_flag = fileID[gtx]['signal_photons']['d_flag']
        classed_pc_flag = fileID[gtx]['signal_photons']['classed_pc_flag']
        # for each ATL03 20 meter segment
        for s1,i1 in zip(segment_id,ph_index_beg):
            # indices within the ATL08 signal photons group
            i2, = np.nonzero(ph_segment_id == s1)
            # full index within the ATL03 file
            indx = i1 + classed_pc_indx[i2]
            # remapped dragann and land ice classifications
            output['d_flag'][indx] = d_flag[i2]
            output['classed_pc_flag'][indx] = classed_pc_flag[i2]
    # close the ATL08 HDF5 file
    fileID.close()
    # return the output variables and attributes
    return (output, attrs)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Acquires the ATL03 geolocated photon height product
            and appends the ATL08 DRAGANN classifications from NSIDC
            """
    )
    # command line parameters
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
    # output subdirectories
    parser.add_argument('--flatten','-F',
        default=False, action='store_true',
        help='Do not create subdirectories')
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
    # permissions mode of the output files (number in octal)
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
    nsidc_icesat2_dragann(args.directory, args.release, args.version,
        args.granule, args.track, YEARS=args.year,
        SUBDIRECTORY=args.subdirectory, CYCLES=args.cycle,
        FLATTEN=args.flatten, TIMEOUT=args.timeout, RETRY=args.retry,
        LOG=args.log, LIST=args.list, CLOBBER=args.clobber,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
