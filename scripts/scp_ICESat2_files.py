#!/usr/bin/env python
u"""
scp_ICESat2_files.py
Written by Tyler Sutterley (05/2022)
Copies ICESat-2 HDF5 data from between a local host and a remote host
can switch between pushing and pulling to/from remote
    PUSH to remote: s.put(local_file, remote_file)
    PULL from remote: s.get(remote_file,local_path=local_file)

CALLING SEQUENCE:
    python scp_ICESat2_files.py --host <host> --user <username> \
        --product ATL06 --release 003 --granule 10 11 12 --cycle 1 2 \
        --remote <path_to_remote> --verbose --mode 0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    --host X: Remote server host
    --user X: Remote server username
    -D X, --directory X: Local working directory
    --remote X: Remote working directory
    --product X: ICESat-2 data product to copy
    --release X: ICESat-2 data release to copy
    --version X: ICESat-2 data version to copy
    --granule X: ICESat-2 granule regions to copy
    --cycle X: ICESat-2 cycle to copy
    --track X: ICESat-2 tracks to copy
    -C, --clobber: overwrite existing data in transfer
    -V, --verbose: output information about each synced file
    --push: Transfer files from local computer to remote server
    -L, --list: only list files to be transferred
    -M X, --mode X: permission mode of directories and files copied

PYTHON DEPENDENCIES:
    paramiko: Native Python SSHv2 protocol library
        http://www.paramiko.org/
        https://github.com/paramiko/paramiko
    scp: scp module for paramiko
        https://github.com/jbardin/scp.py

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 10/2020: using argparse to set parameters
    Updated 05/2020: adjust regular expression to run ATL07 sea ice products
    Updated 09/2019: sort subdirectories.
    Updated 07/2019: using Python3 compliant division.  regex for file versions
    Written 05/2019
"""
from __future__ import print_function, division

import sys
import os
import re
import io
import scp
import getpass
import logging
import argparse
import builtins
import paramiko
import posixpath
import numpy as np

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Copies ICESat-2 HDF5 data from between a local host and
            remote host
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
    #-- remote server credentials
    parser.add_argument('--host','-H',
        type=str, default='',
        help='Hostname of the remote server')
    parser.add_argument('--user','-U',
        type=str, default='',
        help='Remote server username')
    #-- working data directories
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Local working directory')
    parser.add_argument('--remote','-R',
        type=str, default='',
        help='Remote working directory')
    #-- ICESat-2 parameters
    #-- ICESat-2 data product
    parser.add_argument('--product','-p',
        metavar='PRODUCTS', type=str,
        choices=PRODUCTS.keys(), default='ATL06',
        help='ICESat-2 data product to copy')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='004',
        help='ICESat-2 data release to copy')
    #-- ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 data versions to copy')
    #-- ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 granule regions to copy')
    #-- ICESat-2 orbital cycle
    parser.add_argument('--cycle','-c',
        type=int, nargs='+',
        default=range(1,10),
        help='ICESat-2 orbital cycles to copy')
    #-- ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs) to copy')
    #-- sync options
    parser.add_argument('--push','-P',
        default=False, action='store_true',
        help='Transfer files from local computer to remote server')
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    #-- verbose will output information about each copied file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output directories and files')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- use entered host and username
    client_kwds = {}
    client_kwds.setdefault('hostname',args.host)
    client_kwds.setdefault('username',args.user)
    #-- use ssh configuration file to extract hostname, user and identityfile
    user_config_file = os.path.join(os.environ['HOME'],".ssh","config")
    if os.path.exists(user_config_file):
        #-- read ssh configuration file and parse with paramiko
        ssh_config = paramiko.SSHConfig()
        with open(user_config_file) as f:
            ssh_config.parse(f)
        #-- lookup hostname from list of hosts
        user_config = ssh_config.lookup(args.host)
        client_kwds['hostname'] = user_config['hostname']
        #-- get username if not entered from command-line
        if args.user is None and 'username' in user_config.keys():
            client_kwds['username'] = user_config['user']
        #-- use identityfile if in ssh configuration file
        if 'identityfile' in user_config.keys():
            client_kwds['key_filename'] = user_config['identityfile']

    #-- open HOST ssh client for USER (and use password if no IDENTITYFILE)
    client = attempt_login(**client_kwds)
    #-- open secure FTP client
    client_ftp = client.open_sftp()
    #-- verbosity settings
    if args.verbose or args.list:
        logging.getLogger("paramiko").setLevel(logging.INFO)
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.CRITICAL)
    #-- print username for remote client
    logging.info('{0}@{1}:\n'.format(client_kwds['username'],
        client_kwds['hostname']))

    #-- run program
    scp_ICESat2_files(client, client_ftp, args.directory, args.remote,
        args.product, args.release, args.version, args.granule, args.cycle,
        args.track, PUSH=args.push, LIST=args.list, CLOBBER=args.clobber,
        MODE=args.mode)

    #-- close the secure FTP server
    client_ftp.close()
    #-- close the ssh client
    client = None

#-- PURPOSE: try logging onto the server and catch authentication errors
def attempt_login(**client_kwds):
    #-- open HOST ssh client
    kwds = client_kwds.copy()
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    tryagain = True
    #-- add initial attempt
    attempts = 1
    #-- use identification file
    try:
        client.connect(**kwds)
    except paramiko.ssh_exception.AuthenticationException:
        pass
    else:
        return client
    #-- add attempt
    attempts += 1
    #-- phrase for entering password
    phrase = 'Password for {0}@{1}: '.format(kwds['username'],kwds['hostname'])
    #-- remove key_filename from keywords
    kwds.pop('key_filename') if 'key_filename' in kwds.keys() else None
    #-- enter password securely from command-line
    while tryagain:
        kwds['password'] = getpass.getpass(phrase)
        try:
            client.connect(*kwds)
        except paramiko.ssh_exception.AuthenticationException:
            pass
        else:
            kwds.pop('password')
            return client
        #-- retry with new password
        logging.critical('Authentication Failed (Attempt {0:d})'.format(attempts))
        tryagain = builtins.input('Try Different Password? (Y/N): ') in ('Y','y')
        #-- add attempt
        attempts += 1
    #-- exit program if not trying again
    sys.exit()

#-- PURPOSE: copies ICESat-2 HDF5 files between a remote host and a local host
def scp_ICESat2_files(client, client_ftp, DIRECTORY, REMOTE, PRODUCT,
    RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, CLOBBER=False,
    PUSH=False, LIST=False, MODE=0o775):
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    TRACKS = np.arange(1,1388) if not np.any(TRACKS) else TRACKS
    CYCLES = np.arange(1,3) if not np.any(CYCLES) else CYCLES
    GRANULES = np.arange(1,15) if not np.any(GRANULES) else GRANULES
    VERSIONS = np.arange(1,10) if not np.any(VERSIONS) else VERSIONS
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_cycle = '|'.join(['{0:02d}'.format(C) for C in CYCLES])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    #-- compile regular expression operator for finding subdirectories
    #-- and extracting date information from the subdirectory
    rx1 = re.compile(r'(\d+)\.(\d+)\.(\d+)',re.VERBOSE)
    #-- compile regular expression operator for extracting data from files
    args = (PRODUCT,regex_track,regex_cycle,regex_granule,RELEASE,regex_version)
    regex_pattern = (r'(processed_)?({0})(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).h5$')
    rx2 = re.compile(regex_pattern.format(*args,re.VERBOSE))
    #-- if pushing from local directory to remote directory
    if PUSH:
        #-- find all local subdirectories
        SUBDIRECTORY = [s for s in os.listdir(DIRECTORY) if rx1.match(s)]
        #-- for each subdirectory to run
        for sub in sorted(SUBDIRECTORY):
            #-- find files within local directory
            local_dir = os.path.join(DIRECTORY,sub)
            remote_path = os.path.join(DIRECTORY,sub)
            file_list = [f for f in os.listdir(local_dir) if rx2.match(f)]
            for fi in sorted(file_list):
                #-- check if data directory exists and recursively create if not
                remote_makedirs(client_ftp, remote_path, LIST=LIST, MODE=MODE)
                #-- push file from local to remote
                scp_push_file(client, client_ftp, fi, local_dir, remote_path,
                    CLOBBER=CLOBBER, LIST=LIST, MODE=MODE)
    else:
        #-- find all remote subdirectories
        SUBDIRECTORY = [s for s in client_ftp.listdir(REMOTE) if rx1.match(s)]
        #-- for each subdirectory to run
        for sub in sorted(SUBDIRECTORY):
            #-- local and remote directories
            local_dir = os.path.join(DIRECTORY,sub)
            remote_path = posixpath.join(REMOTE,sub)
            #-- find remote files for hemisphere
            file_list=[f for f in client_ftp.listdir(remote_path) if rx2.match(f)]
            for fi in sorted(file_list):
                #-- check if data directory exists and recursively create if not
                if not os.access(local_dir, os.F_OK) and not LIST:
                    os.makedirs(local_dir, MODE)
                #-- push file from local to remote
                scp_pull_file(client, client_ftp, fi, local_dir, remote_path,
                    CLOBBER=CLOBBER, LIST=LIST, MODE=MODE)

#-- PURPOSE: recursively create directories on remote server
def remote_makedirs(client_ftp, remote_dir, LIST=False, MODE=0o775):
    dirs = remote_dir.split(posixpath.sep)
    remote_path = dirs[0] if dirs[0] else posixpath.sep
    for s in dirs:
        if (s not in client_ftp.listdir(remote_path)) and not LIST:
            client_ftp.mkdir(posixpath.join(remote_path,s), MODE)
        remote_path = posixpath.join(remote_path,s)

#-- PURPOSE: push a local file to a remote host checking if file exists
#-- and if the local file is newer than the remote file (reprocessed)
#-- set the permissions mode of the remote transferred file to MODE
def scp_push_file(client, client_ftp, transfer_file, local_dir, remote_dir,
    CLOBBER=False, LIST=False, MODE=0o775):
    #-- local and remote versions of file
    local_file = os.path.join(local_dir,transfer_file)
    remote_file = posixpath.join(remote_dir,transfer_file)
    #-- check if local file is newer than the remote file
    TEST = False
    OVERWRITE = 'clobber'
    if (transfer_file in client_ftp.listdir(remote_dir)):
        local_mtime = os.stat(local_file).st_mtime
        remote_mtime = client_ftp.stat(remote_file).st_mtime
        #-- if local file is newer: overwrite the remote file
        if (even(local_mtime) > even(remote_mtime)):
            TEST = True
            OVERWRITE = 'overwrite'
    else:
        TEST = True
        OVERWRITE = 'new'
    #-- if file does not exist remotely, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        logging.info('{0} --> '.format(local_file))
        logging.info('\t{0} ({1})\n'.format(remote_file,OVERWRITE))
        #-- if not only listing files
        if not LIST:
            #-- copy local files to remote server
            with scp.SCPClient(client.get_transport(), socket_timeout=20) as s:
                s.put(local_file, remote_file, preserve_times=True)
            #-- change the permissions level of the transported file to MODE
            client_ftp.chmod(remote_file, MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file (reprocessed)
#-- set the permissions mode of the local transferred file to MODE
def scp_pull_file(client, client_ftp, transfer_file, local_dir, remote_dir,
    CLOBBER=False, LIST=False, MODE=0o775):
    #-- local and remote versions of file
    local_file = os.path.join(local_dir,transfer_file)
    remote_file = posixpath.join(remote_dir,transfer_file)
    #-- check if remote file is newer than the local file
    TEST = False
    OVERWRITE = 'clobber'
    if os.access(local_file, os.F_OK):
        local_mtime = os.stat(local_file).st_mtime
        remote_mtime = client_ftp.stat(remote_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (even(remote_mtime) > even(local_mtime)):
            TEST = True
            OVERWRITE = 'overwrite'
    else:
        TEST = True
        OVERWRITE = 'new'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        logging.info('{0} --> '.format(remote_file))
        logging.info('\t{0} ({1})\n'.format(local_file,OVERWRITE))
        #-- if not only listing files
        if not LIST:
            #-- copy local files from remote server
            with scp.SCPClient(client.get_transport(), socket_timeout=20) as s:
                s.get(remote_file, local_path=local_file, preserve_times=True)
            #-- change the permissions level of the transported file to MODE
            os.chmod(local_file, MODE)

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(i):
    return 2*int(i//2)

#-- run main program
if __name__ == '__main__':
    main()
