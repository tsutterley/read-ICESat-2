#!/usr/bin/env python
u"""
scp_ICESat2_files.py
Written by Tyler Sutterley (09/2019)
Copies ICESat-2 HDF5 data from between a local host and a remote host
can switch between pushing and pulling to/from remote
    PUSH to remote: s.put(local_file, remote_file)
    PULL from remote: s.get(remote_file,local_path=local_file)

CALLING SEQUENCE:
    python scp_ICESat2_files.py --host=<host> --user=<username> \
        --product=ATL06 --release=205 --granule=10,11,12 --cycle=1,2 \
        --remote=<path_to_remote> --verbose --mode=0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    --host=X: Remote server host
    --user=X: Remote server username
    -D X, --directory=X: Local working directory
    --remote=X: Remote working directory
    --product=X: ICESat-2 data product to copy
    --release=X: ICESat-2 data release to copy
    --version=X: ICESat-2 data version to copy
    --granule=X: ICESat-2 granule regions to copy
    --cycle=X: ICESat-2 cycle to copy
    --track=X: ICESat-2 tracks to copy
    -C, --clobber: overwrite existing data in transfer
    -V, --verbose: output information about each synced file
    -M X, --mode=X: permission mode of directories and files synced
    --push: Transfer files from local computer to remote server
    -L, --list: only list files to be transferred

PYTHON DEPENDENCIES:
    paramiko: Native Python SSHv2 protocol library
        http://www.paramiko.org/
        https://github.com/paramiko/paramiko
    scp: scp module for paramiko
        https://github.com/jbardin/scp.py

UPDATE HISTORY:
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
    print(' -D X, --directory=X\tLocal working directory')
    print(' --remote=X\t\tRemote working directory')
    print(' --product=X\t\tICESat-2 data product to copy')
    print(' --release=X\t\tICESat-2 data release to copy')
    print(' --version=X\t\tICESat-2 data version to copy')
    print(' --granule=X\t\tICESat-2 granule regions to copy')
    print(' --cycle=X\t\tICESat-2 cycles to copy')
    print(' --track=X\t\tICESat-2 tracks to copy')
    print(' -C, --clobber\t\tOverwrite existing data in transfer')
    print(' -V, --verbose\t\tOutput information about each created file')
    print(' -M X, --mode=X\t\tPermission mode of directories and files created')
    print(' --push\t\t\tTransfer files from local computer to remote server')
    print(' -L, --list\t\tOnly list files to be transferred\n')

#-- Main program that calls scp_ICESat2_files()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','host=','user=','directory=','remote=','product=',
        'release=','version=','granule=','cycle=','track=','verbose','clobber',
        'mode=','push','list']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:VCM:L',long_options)

    #-- command line parameters
    HOST = ''
    USER = None
    IDENTITYFILE = None
    #-- working data directories
    DIRECTORY = os.getcwd()
    REMOTE = ''
    #-- ICESat-2 parameters
    PRODUCT = 'ATL06'
    RELEASE = '002'
    VERSIONS = None
    CYCLES = np.arange(1,6)
    GRANULES = None
    TRACKS = None
    VERBOSE = False
    CLOBBER = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    PUSH = False
    LIST = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("--host"):
            HOST = arg
        elif opt in ("--user"):
            USER = arg
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--remote"):
            REMOTE = arg
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
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-C","--clobber"):
            CLOBBER = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("--push"):
            PUSH = True
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
        #-- get username if not entered from command-line
        if USER is None and 'username' in user_config.keys():
            USER = user_config['username']
        #-- use identityfile if in ssh configuration file
        if 'identityfile' in user_config.keys():
            IDENTITYFILE = user_config['identityfile']

    #-- open HOST ssh client for USER (and use password if no IDENTITYFILE)
    client = attempt_login(HOST, USER, IDENTITYFILE=IDENTITYFILE)

    #-- open secure FTP client
    client_ftp = client.open_sftp()
    #-- verbosity settings
    if VERBOSE or LIST:
        logging.getLogger("paramiko").setLevel(logging.WARNING)
        print('{0}@{1}:\n'.format(USER, HOST))

    #-- run program
    scp_ICESat2_files(client, client_ftp, DIRECTORY, REMOTE, PRODUCT,
        RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, CLOBBER=CLOBBER,
        VERBOSE=VERBOSE, PUSH=PUSH, LIST=LIST, MODE=MODE)

    #-- close the secure FTP server
    client_ftp.close()
    #-- close the ssh client
    client = None

#-- PURPOSE: try logging onto the server and catch authentication errors
def attempt_login(HOST, USER, IDENTITYFILE=None):
    #-- open HOST ssh client
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    tryagain = True
    attempts = 1
    #-- use identification file
    if IDENTITYFILE:
        try:
            client.connect(HOST, username=USER, key_filename=IDENTITYFILE)
        except paramiko.ssh_exception.AuthenticationException:
            pass
        else:
            return client
        attempts += 1
    #-- enter password securely from command-line
    while tryagain:
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))
        try:
            client.connect(HOST, username=USER, password=PASSWORD)
        except paramiko.ssh_exception.AuthenticationException:
            pass
        else:
            del PASSWORD
            return client
        #-- retry with new password
        print('Authentication Failed (Attempt {0:d})'.format(attempts))
        tryagain = builtins.input('Try Different Password? (Y/N): ') in ('Y','y')
        attempts += 1
    #-- exit program if not trying again
    sys.exit()

#-- PURPOSE: copies ICESat-2 HDF5 files between a remote host and a local host
def scp_ICESat2_files(client, client_ftp, DIRECTORY, REMOTE, PRODUCT,
    RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, CLOBBER=False, VERBOSE=False,
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
    rx1 = re.compile('(\d+)\.(\d+)\.(\d+)',re.VERBOSE)
    #-- compile regular expression operator for extracting data from files
    args = (PRODUCT,regex_track,regex_cycle,regex_granule,RELEASE,regex_version)
    rx2 = re.compile(('({0})(.*?)_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})(\d{{2}})'
        '(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).h5$'.format(*args)))

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
                    CLOBBER=CLOBBER, VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)
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
                    CLOBBER=CLOBBER, VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)

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
    CLOBBER=False, VERBOSE=False, LIST=False, MODE=0o775):
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
        if VERBOSE or LIST:
            print('{0} --> '.format(local_file))
            print('\t{0} ({1})\n'.format(remote_file,OVERWRITE))
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
    CLOBBER=False, VERBOSE=False, LIST=False, MODE=0o775):
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
        if VERBOSE or LIST:
            print('{0} --> '.format(remote_file))
            print('\t{0} ({1})\n'.format(local_file,OVERWRITE))
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
