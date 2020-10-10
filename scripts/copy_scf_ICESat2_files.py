#!/usr/bin/env python
u"""
copy_scf_ICESat2_files.py
Written by Tyler Sutterley (10/2020)
Copies ICESat-2 HDF5 files from the SCF server

CALLING SEQUENCE:
    python copy_scf_ICESat2_files.py --scf_host <host> --scf_user <username> \
        --product ATL06 --release 003 --granule 10 11 12 --cycle 1 2 \
        --scf_outgoing <path_to_outgoing> --verbose --mode 0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    --scf_host X: hostname of the SCF server
    --scf_user X: SCF server username
    -D X, --directory X: local working directory for receiving data
    --product X: ICESat-2 data product to copy
    --release X: ICESat-2 data release to copy
    --version X: ICESat-2 data version to copy
    --granule X: ICESat-2 granule regions to copy
    --cycle X: ICESat-2 cycle to copy
    --track X: ICESat-2 tracks to copy
     --scf_incoming X: directory on the SCF where the rscf sends PANS
     --scf_outgoing X: directory on the SCF where the data resides
    -C, --clobber: overwrite existing data in transfer
    -V, --verbose: output information about each synced file
    -M X, --mode X: permission mode of directories and files synced
    -L, --list: only list files to be transferred

PYTHON DEPENDENCIES:
    paramiko: Native Python SSHv2 protocol library
        http://www.paramiko.org/
        https://github.com/paramiko/paramiko

UPDATE HISTORY:
    Updated 10/2020: using argparse to set parameters
    Updated 05/2020: adjust regular expression to run ATL07 sea ice products
    Updated 07/2019: using python3 compliant division
    Updated 05/2019: changed host and user to scf_host and scf_user
    Updated 04/2019: added parameters for selecting the ICESat-2 product,
        release, version, granule, cycle and track
    Written 04/2019
"""
from __future__ import print_function, division

import sys
import os
import re
import logging
import argparse
import paramiko
import posixpath
import numpy as np

#-- Main program that calls copy_scf_files()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Copies ICESat-2 HDF5 files from the SCF server
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
    #-- SCF credentials credentials
    parser.add_argument('--scf_host',
        type=str, default='',
        help='Hostname of the SCF server')
    parser.add_argument('--scf_user',
        type=str, default='',
        help='SCF server username')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- years of data to copy
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to copy')
    #-- subdirectories of data to copy
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to copy')
    #-- ICESat-2 parameters
    #-- ICESat-2 data product
    parser.add_argument('--product','-p',
        metavar='PRODUCTS', type=str,
        choices=PRODUCTS.keys(), default='ATL06',
        help='ICESat-2 data product to copy')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='003',
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
    #-- ICESat-2 Science Computing Facility (SCF) parameters
    parser.add_argument('--scf_incoming',
        type=str,
        help='Directory on the SCF where the rscf sends PANS')
    parser.add_argument('--scf_outgoing',
        type=str,
        help='Directory on the SCF where the data resides')
    #-- sync options
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
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- use entered host and username
    scf_kwds = {}
    scf_kwds.setdefault('hostname',args.scf_host)
    scf_kwds.setdefault('username',args.scf_user)
    #-- use ssh configuration file to extract hostname, user and identityfile
    user_config_file = os.path.join(os.environ['HOME'],".ssh","config")
    if os.path.exists(user_config_file):
        #-- read ssh configuration file and parse with paramiko
        ssh_config = paramiko.SSHConfig()
        with open(user_config_file) as f:
            ssh_config.parse(f)
        #-- lookup hostname from list of hosts
        scf_user_config = ssh_config.lookup(args.scf_host)
        scf_kwds['hostname'] = scf_user_config['hostname']
        #-- get username if not entered from command-line
        if args.scf_user is None and 'username' in scf_user_config.keys():
            scf_kwds['username'] = scf_user_config['user']
        #-- use identityfile if in ssh configuration file
        if 'identityfile' in scf_user_config.keys():
            scf_kwds['key_filename'] = scf_user_config['identityfile']

    #-- open HOST ssh client for USER
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.connect(**scf_kwds)

    #-- open secure FTP client
    client_ftp = client.open_sftp()
    #-- verbosity settings
    if args.verbose or args.list:
        logging.getLogger("paramiko").setLevel(logging.WARNING)
        print('{0}@{1}:\n'.format(scf_kwds['username'],scf_kwds['hostname']))

    #-- run SCF copy program
    copy_scf_files(client, client_ftp, args.directory, args.scf_incoming,
        args.scf_outgoing, args.product, args.release, args.version,
        args.granule, args.cycle, args.track, CLOBBER=args.clobber,
        VERBOSE=args.verbose, LIST=args.list, MODE=args.mode)

    #-- close the secure FTP server
    client_ftp.close()
    #-- close the ssh client
    client = None

#-- PURPOSE: copy ICESat-2 files to data directory with data subdirectories
def copy_scf_files(client, client_ftp, base_dir, scf_incoming, scf_outgoing,
    PRODUCT, RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, CLOBBER=False,
    VERBOSE=False, LIST=False, MODE=0o775):
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    TRACKS = np.arange(1,1388) if not np.any(TRACKS) else TRACKS
    CYCLES = np.arange(1,3) if not np.any(CYCLES) else CYCLES
    GRANULES = np.arange(1,15) if not np.any(GRANULES) else GRANULES
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_cycle = '|'.join(['{0:02d}'.format(C) for C in CYCLES])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    #-- compile regular expression operator for extracting data from files
    args = (PRODUCT,regex_track,regex_cycle,regex_granule,RELEASE,regex_version)
    regex_pattern = (r'(processed_)?({0})(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).h5$')
    rx = re.compile(regex_pattern.format(*args),re.VERBOSE)
    #-- find files within scf_outgoing
    file_transfers = [f for f in client_ftp.listdir(scf_outgoing) if rx.match(f)]
    for f in file_transfers:
        #-- extract parameters from file
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(f).pop()
        #-- local directory set by product
        #-- check if data directory exists and recursively create if not
        local_dir = os.path.join(base_dir,'{0}.{1}.{2}'.format(YY,MM,DD))
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- pull file from remote to local
        scp_pull_file(client, client_ftp, f, local_dir, scf_outgoing,
            CLOBBER=CLOBBER, VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
#-- set the permissions mode of the local transferred file to MODE
def scp_pull_file(client, client_ftp, transfer_file, local_dir, remote_dir,
    CLOBBER=False, VERBOSE=False, LIST=False, MODE=0o775):
    #-- local and remote versions of file
    local_file = os.path.join(local_dir,transfer_file)
    remote_file = posixpath.join(remote_dir,transfer_file)
    #-- get access and modification time of remote file
    remote_atime = client_ftp.stat(remote_file).st_atime
    remote_mtime = client_ftp.stat(remote_file).st_mtime
    #-- check if remote file is newer than the local file
    TEST = False
    OVERWRITE = 'clobber'
    if os.access(local_file, os.F_OK):
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (even(remote_mtime) > even(local_mtime)):
            TEST = True
            OVERWRITE = 'overwrite'
    else:
        TEST = True
        OVERWRITE = 'new'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        if VERBOSE or LIST:
            print('{0} --> '.format(remote_file))
            print('\t{0} ({1})\n'.format(local_file,OVERWRITE))
        #-- if not only listing files
        if not LIST:
            #-- copy local files from remote server
            client_ftp.get(remote_file, local_file)
            #-- keep modification and access time of input file
            os.utime(local_file, (remote_atime, remote_mtime))
            #-- change the permissions level of the transported file to MODE
            os.chmod(local_file, MODE)

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(i):
    return 2*int(i//2)

#-- run main program
if __name__ == '__main__':
    main()
