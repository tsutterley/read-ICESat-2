#!/usr/bin/env python
u"""
scp_scf_ICESat2_files.py
Written by Tyler Sutterley (05/2020)
Copies ICESat-2 HDF5 files from the SCF server to a remote host

CALLING SEQUENCE:
    python scp_scf_ICESat2_files.py --host=<host> --user=<username> \
        --scf_host=<scf_host> --scf_user=<scf_username> \
        --product=ATL06 --release=205 --granule=10,11,12 --cycle=1,2 \
        --remote=<path_to_remote> --scf_outgoing=<path_to_outgoing> \
        --verbose --mode=0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    --host=X: Remote server host
    --user=X: Remote server username
    --scf_host=X: hostname of the SCF server
    --scf_user=X: SCF server username
    --remote=X: Remote working directory for receiving data
    --product=X: ICESat-2 data product to copy
    --release=X: ICESat-2 data release to copy
    --version=X: ICESat-2 data version to copy
    --granule=X: ICESat-2 granule regions to copy
    --cycle=X: ICESat-2 cycle to copy
    --track=X: ICESat-2 tracks to copy
    --scf_incoming=X: directory on the SCF where the rscf sends PANS
    --scf_outgoing=X: directory on the SCF where the data resides
    -C, --clobber: overwrite existing data in transfer
    -V, --verbose: output information about each synced file
    -M X, --mode=X: permission mode of directories and files synced
    -L, --list: only list files to be transferred

PYTHON DEPENDENCIES:
    paramiko: Native Python SSHv2 protocol library
        http://www.paramiko.org/
        https://github.com/paramiko/paramiko
    scp: scp module for paramiko
        https://github.com/jbardin/scp.py

UPDATE HISTORY:
    Updated 05/2020: adjust regular expression to run ATL07 sea ice products
    Updated 07/2019: using Python3 compliant division.  regex for file versions
    Updated 05/2019: only create directories if not --list
    Written 05/2019
"""
from __future__ import print_function, division

import sys
import os
import re
import io
import scp
import getopt
import getpass
import logging
import paramiko
import posixpath
import numpy as np

#-- PURPOSE: help module to describe the optional input command-line parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' --host=X\t\tRemote server host')
    print(' --user=X\t\tRemote server user')
    print(' --scf_host=X\t\tHostname of the SCF server')
    print(' --scf_user=X\t\tSCF server username')
    print(' --remote=X\t\tRemote working directory for receiving data')
    print(' --product=X\t\tICESat-2 data product to copy')
    print(' --release=X\t\tICESat-2 data release to copy')
    print(' --version=X\t\tICESat-2 data version to copy')
    print(' --granule=X\t\tICESat-2 granule regions to copy')
    print(' --cycle=X\t\tICESat-2 cycles to copy')
    print(' --track=X\t\tICESat-2 tracks to copy')
    print(' --scf_incoming=X\tDirectory on the SCF where the rscf sends PANS')
    print(' --scf_outgoing=X\tDirectory on the SCF where the data resides')
    print(' -C, --clobber\t\tOverwrite existing data in transfer')
    print(' -V, --verbose\t\tOutput information about each created file')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' -L, --list\t\tOnly list files to be transferred\n')

#-- Main program that calls scp_scf_files()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','host=','user=','scf_host=','scf_user=','remote=',
        'product=','release=','version=','granule=','cycle=','track=',
        'scf_incoming=','scf_outgoing=','verbose','clobber','mode=','list']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hVCM:L',long_options)

    #-- command line parameters
    HOST = ''
    USER = None
    SCF_HOST = ''
    SCF_USER = None
    IDENTITYFILE = None
    SCF_IDENTITY = None
    #-- working data directories
    remote_dir = ''
    #-- ICESat-2 parameters
    PRODUCT = 'ATL06'
    RELEASE = '002'
    VERSIONS = None
    CYCLES = np.arange(1,6)
    GRANULES = None
    TRACKS = None
    scf_incoming = ''
    scf_outgoing = ''
    VERBOSE = False
    CLOBBER = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    LIST = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("--host"):
            HOST = arg
        elif opt in ("--user"):
            USER = arg
        elif opt in ("--scf_host"):
            SCF_HOST = arg
        elif opt in ("--scf_user"):
            SCF_USER = arg
        elif opt in ("--remote"):
            remote_dir = arg
        elif opt in ("--product"):
            PRODUCT = arg
        elif opt in ("--release"):
            RELEASE = arg
        elif opt in ("--version"):
            VERSIONS = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--granule"):
            GRANULES = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--cycle"):
            CYCLES = np.sort(arg.split(',')).astype(np.int)
        elif opt in ("--track"):
            TRACKS = np.sort(arg.split(',')).astype(np.int)
        elif opt in ("--scf_incoming"):
            scf_incoming = arg
        elif opt in ("--scf_outgoing"):
            scf_outgoing = arg
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-C","--clobber"):
            CLOBBER = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-L","--list"):
            LIST = True

    #-- use ssh configuration file to extract hostname, user and identityfile
    user_config_file = os.path.join(os.environ['HOME'],".ssh","config")
    if os.path.exists(user_config_file):
        #-- read ssh configuration file and parse with paramiko
        ssh_config = paramiko.SSHConfig()
        with open(user_config_file) as f:
            ssh_config.parse(f)
        #-- lookup hostname from list of hosts
        user_config = ssh_config.lookup(HOST)
        HOST = user_config['hostname']
        scf_user_config = ssh_config.lookup(SCF_HOST)
        SCF_HOST = scf_user_config['hostname']
        #-- get username if not entered from command-line
        if USER is None and 'username' in user_config.keys():
            USER = user_config['username']
        if SCF_USER is None and 'username' in scf_user_config.keys():
            SCF_USER = user_config['username']
        #-- use identityfile if in ssh configuration file
        if 'identityfile' in user_config.keys():
            IDENTITYFILE = user_config['identityfile']
        if 'identityfile' in scf_user_config.keys():
            SCF_IDENTITY = user_config['identityfile']

    #-- open HOST ssh client for USER
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.connect(HOST, username=USER, key_filename=IDENTITYFILE)
    #-- open SCF_HOST ssh client for SCF_USER
    scf_client = paramiko.SSHClient()
    scf_client.load_system_host_keys()
    scf_client.connect(SCF_HOST, username=SCF_USER, key_filename=SCF_IDENTITY)

    #-- open secure FTP client
    client_ftp = client.open_sftp()
    #-- open secure FTP scf_client
    scf_client_ftp = scf_client.open_sftp()
    #-- verbosity settings
    if VERBOSE or LIST:
        logging.getLogger("paramiko").setLevel(logging.WARNING)
        print('{0}@{1} --> {2}@{3}\n'.format(SCF_USER, SCF_HOST, USER, HOST))

    #-- run program
    scp_scf_files(client, client_ftp, scf_client, scf_client_ftp, remote_dir,
        scf_incoming, scf_outgoing, PRODUCT, RELEASE, VERSIONS, GRANULES, CYCLES,
        TRACKS, CLOBBER=CLOBBER, VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)

    #-- close the secure FTP server
    client_ftp.close()
    scf_client_ftp.close()
    #-- close the ssh clients
    client = None
    scf_client = None

#-- PURPOSE: copy ICESat-2 files to data directory with data subdirectories
def scp_scf_files(client, client_ftp, scf_client, scf_client_ftp, remote_dir,
    scf_incoming, scf_outgoing, PRODUCT, RELEASE, VERSIONS, GRANULES, CYCLES,
    TRACKS, CLOBBER=False, VERBOSE=False, LIST=False, MODE=0o775):
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    TRACKS = np.arange(1,1388) if not np.any(TRACKS) else TRACKS
    CYCLES = np.arange(1,3) if not np.any(CYCLES) else CYCLES
    GRANULES = np.arange(1,15) if not np.any(GRANULES) else GRANULES
    VERSIONS = np.arange(1,10) if not np.any(VERSIONS) else VERSIONS
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
    file_list = [f for f in scf_client_ftp.listdir(scf_outgoing) if rx.match(f)]
    for f in sorted(file_list):
        #-- extract parameters from file
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(f).pop()
        #-- put data in directories similar to NSIDC
        #-- check if data directory exists and recursively create if not
        remote_path = posixpath.join(remote_dir,'{0}.{1}.{2}'.format(YY,MM,DD))
        remote_makedirs(client_ftp, remote_path, LIST=LIST, MODE=MODE)
        #-- pull file from scf to remote
        scp_pull_file(client_ftp, scf_client_ftp, f, remote_path, scf_outgoing,
            CLOBBER=CLOBBER, VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)

#-- PURPOSE: recursively create directories on remote server
def remote_makedirs(client_ftp, remote_dir, LIST=False, MODE=0o775):
    dirs = remote_dir.split(posixpath.sep)
    remote_path = dirs[0] if dirs[0] else posixpath.sep
    for s in dirs:
        if (s not in client_ftp.listdir(remote_path)) and not LIST:
            client_ftp.mkdir(posixpath.join(remote_path,s), MODE)
        remote_path = posixpath.join(remote_path,s)

#-- PURPOSE: pull file from scf host checking if file exists locally
#-- and if the remote file is newer than the local file
#-- set the permissions mode of the local transferred file to MODE
def scp_pull_file(client_ftp, scf_client_ftp, transfer_file, remote_dir,
    scf_outgoing, CLOBBER=False, VERBOSE=False, LIST=False, MODE=0o775):
    #-- remote and scf outgoing versions of file
    remote_file = os.path.join(remote_dir,transfer_file)
    outgoing_file = posixpath.join(scf_outgoing,transfer_file)
    #-- get access and modification time of remote file
    outgoing_atime = scf_client_ftp.stat(outgoing_file).st_atime
    outgoing_mtime = scf_client_ftp.stat(outgoing_file).st_mtime
    #-- check if remote file is newer than the local file
    TEST = False
    OVERWRITE = 'clobber'
    if (transfer_file in client_ftp.listdir(remote_dir)):
        remote_mtime = client_ftp.stat(remote_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (even(outgoing_mtime) > even(remote_mtime)):
            TEST = True
            OVERWRITE = 'overwrite'
    else:
        TEST = True
        OVERWRITE = 'new'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        if VERBOSE or LIST:
            print('{0} --> '.format(outgoing_file))
            print('\t{0} ({1})\n'.format(remote_file,OVERWRITE))
        #-- if not only listing files
        if not LIST:
            #-- load scf file contents to BytesIO object
            fileobj = io.BytesIO()
            scf_client_ftp.getfo(outgoing_file, fileobj)
            #-- rewind retrieved binary to start of file
            fileobj.seek(0)
            #-- copy BytesIO object to remote server
            client_ftp.putfo(fileobj, remote_file)
            #-- change the permissions level of the transported file to MODE
            client_ftp.chmod(remote_file, MODE)
            #-- keep modification and access time of scf file
            client_ftp.utime(remote_file, (outgoing_atime, outgoing_mtime))

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(i):
    return 2*int(i//2)

#-- run main program
if __name__ == '__main__':
    main()
