#!/usr/bin/env python
u"""
copy_scf_ICESat2_files.py
Written by Tyler Sutterley (05/2019)
Copies ICESat-2 HDF5 files from the SCF server

CALLING SEQUENCE:
	python copy_scf_ICESat2_files.py --scf_host=<host> --scf_user=<username> \
		--product=ATL06 --release=205 --granule=10,11,12 --cycle=1,2 \
		--scf_outgoing=<path_to_outgoing> --verbose --mode=0o775

COMMAND LINE OPTIONS:
	-h, --help: list the command line options
	--scf_host=X: hostname of the SCF server
	--scf_user=X: SCF server username
	-D X, --directory=X: local working directory for receiving data
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

UPDATE HISTORY:
	Updated 05/2019: changed host and user to scf_host and scf_user
	Updated 04/2019: added parameters for selecting the ICESat-2 product,
		release, version, granule, cycle and track
	Written 04/2019
"""
from __future__ import print_function

import sys
import os
import re
import getopt
import getpass
import logging
import paramiko
import posixpath
import numpy as np

#-- PURPOSE: help module to describe the optional input command-line parameters
def usage():
	print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
	print(' --scf_host=X\t\tHostname of the SCF server')
	print(' --scf_user=X\t\tSCF server username')
	print(' -D X, --directory=X\tLocal working directory for receiving data')
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

#-- Main program that calls copy_scf_files()
def main():
	#-- Read the system arguments listed after the program
	long_options = ['help','scf_host=','scf_user=','directory=','product=',
		'release=','version=','granule=','cycle=','track=','scf_incoming=',
		'scf_outgoing=','verbose','clobber','mode=','list']
	optlist,arglist = getopt.getopt(sys.argv[1:],'hD:VCM:L',long_options)

	#-- command line parameters
	HOST = ''
	USER = None
	IDENTITYFILE = None
	#-- working data directories
	base_dir = os.getcwd()
	#-- ICESat-2 parameters
	PRODUCT = 'ATL06'
	RELEASE = '203'
	VERSION = '01'
	CYCLES = np.arange(1,3)
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
		elif opt in ("--scf_host"):
			HOST = arg
		elif opt in ("--scf_user"):
			USER = arg
		elif opt in ("-D","--directory"):
			base_dir = os.path.expanduser(arg)
		elif opt in ("--product"):
			PRODUCT = arg
		elif opt in ("--release"):
			RELEASE = arg
		elif opt in ("--version"):
			VERSION = arg
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
		#-- get username if not entered from command-line
		if USER is None and 'username' in user_config.keys():
			USER = user_config['username']
		#-- use identityfile if in ssh configuration file
		if 'identityfile' in user_config.keys():
			IDENTITYFILE = user_config['identityfile']

	#-- open HOST ssh client for USER
	client = paramiko.SSHClient()
	client.load_system_host_keys()
	client.connect(HOST, username=USER, key_filename=IDENTITYFILE)

	#-- open secure FTP client
	client_ftp = client.open_sftp()
	#-- verbosity settings
	if VERBOSE or LIST:
		logging.getLogger("paramiko").setLevel(logging.WARNING)
		print('{0}@{1}:\n'.format(USER, HOST))

	#-- run program
	copy_scf_files(client, client_ftp, base_dir, scf_incoming, scf_outgoing,
		PRODUCT, RELEASE, VERSION, GRANULES, CYCLES, TRACKS, CLOBBER=CLOBBER,
		VERBOSE=VERBOSE, LIST=LIST, MODE=MODE)

	#-- close the secure FTP server
	client_ftp.close()
	#-- close the ssh client
	client = None

#-- PURPOSE: copy ICESat-2 files to data directory with data subdirectories
def copy_scf_files(client, client_ftp, base_dir, scf_incoming, scf_outgoing,
	PRODUCT, RELEASE, VERSION, GRANULES, CYCLES, TRACKS, CLOBBER=False,
	VERBOSE=False, LIST=False, MODE=0o775):
	#-- find ICESat-2 HDF5 files in the subdirectory for product and release
	TRACKS = np.arange(1,1388) if not np.any(TRACKS) else TRACKS
	CYCLES = np.arange(1,3) if not np.any(CYCLES) else CYCLES
	GRANULES = np.arange(1,15) if not np.any(GRANULES) else GRANULES
	regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
	regex_cycle = '|'.join(['{0:02d}'.format(C) for C in CYCLES])
	regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
	#-- compile regular expression operator for extracting data from files
	args = (PRODUCT,regex_track,regex_cycle,regex_granule,RELEASE,VERSION)
	rx = re.compile(('({0})_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})(\d{{2}})'
		'(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).h5$'.format(*args)))
	#-- find files within scf_outgoing
	file_transfers = [f for f in client_ftp.listdir(scf_outgoing) if rx.match(f)]
	for fi in file_transfers:
		#-- extract parameters from file
		PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX = rx.findall(fi).pop()
		#-- local directory set by product
		#-- check if data directory exists and recursively create if not
		local_dir = os.path.join(base_dir,'{0}.{1}.{2}'.format(YY,MM,DD))
		os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
		#-- pull file from remote to local
		scp_pull_file(client, client_ftp, fi, local_dir, scf_outgoing,
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
	return 2*int(i/2)

#-- run main program
if __name__ == '__main__':
	main()
