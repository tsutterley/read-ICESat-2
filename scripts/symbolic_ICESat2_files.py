#!/usr/bin/env python
u"""
symbolic_ICESat2_files.py
Written by Tyler Sutterley (07/2019)
Creates symbolic links for ICESat-2 HDF5 files organized by date

CALLING SEQUENCE:
    python symbolic_ICESat2_files.py --product=ATL06 --release=001 \
        --granule=10,11,12 --cycle=1,2 --directory=<path_to_directory>
        --scf_outgoing=<path_to_outgoing> --verbose --mode=0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    -D X, --directory=X: local working directory for creating symbolic links
    --product=X: ICESat-2 data product to create symbolic links
    --release=X: ICESat-2 data release to create symbolic links
    --version=X: ICESat-2 data version to create symbolic links
    --granule=X: ICESat-2 granule regions to create symbolic links
    --cycle=X: ICESat-2 cycle to create symbolic links
    --track=X: ICESat-2 tracks to create symbolic links
    --scf_incoming=X: directory on the SCF where the rscf sends PANS
    --scf_outgoing=X: directory on the SCF where the data resides
    -V, --verbose: output information about each symbolic link
    -M X, --mode=X: permission mode of directories

UPDATE HISTORY:
    Updated 05/2020: adjust regular expression to run ATL07 sea ice products
    Written 07/2019
"""
from __future__ import print_function

import sys
import os
import re
import getopt
import numpy as np

#-- PURPOSE: help module to describe the optional input command-line parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tLocal working directory for symbolic links')
    print(' --product=X\t\tICESat-2 data product to create symbolic links')
    print(' --release=X\t\tICESat-2 data release to create symbolic links')
    print(' --version=X\t\tICESat-2 data version to create symbolic links')
    print(' --granule=X\t\tICESat-2 granule regions to create symbolic links')
    print(' --cycle=X\t\tICESat-2 cycles to create symbolic links')
    print(' --track=X\t\tICESat-2 tracks to create symbolic links')
    print(' --scf_incoming=X\tDirectory on the SCF where the rscf sends PANS')
    print(' --scf_outgoing=X\tDirectory on the SCF where the data resides')
    print(' -V, --verbose\t\tOutput information about each symbolic link')
    print(' -M X, --mode=X\t\tPermission mode of directories\n')

#-- Main program that calls symbolic_ICESat2_files()
def main():
    #-- Read the system arguments listed after the program
    long_opt = ['help','directory=','product=','release=','version=','granule=',
        'cycle=','track=','scf_incoming=','scf_outgoing=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:VM:',long_opt)

    #-- command line parameters
    #-- working data directories
    base_dir = os.getcwd()
    #-- ICESat-2 parameters
    PRODUCT = 'ATL06'
    RELEASE = '002'
    VERSIONS = None
    CYCLES = None
    GRANULES = None
    TRACKS = None
    scf_incoming = ''
    scf_outgoing = ''
    VERBOSE = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            base_dir = os.path.expanduser(arg)
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
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- run program
    symbolic_ICESat2_files(base_dir, scf_incoming, scf_outgoing, PRODUCT,
        RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, VERBOSE=VERBOSE, MODE=MODE)

#-- PURPOSE: copy ICESat-2 files to data directory with data subdirectories
def symbolic_ICESat2_files(base_dir, scf_incoming, scf_outgoing, PRODUCT,
    RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, VERBOSE=False, MODE=0o775):
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    TRACKS = np.arange(1,1388) if not np.any(TRACKS) else TRACKS
    CYCLES = np.arange(1,10) if not np.any(CYCLES) else CYCLES
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
    file_transfers = [f for f in os.listdir(scf_outgoing) if rx.match(f)]
    for f in sorted(file_transfers):
        #-- extract parameters from file
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(f).pop()
        #-- put symlinks in directories similar to NSIDC
        #-- check if data directory exists and recursively create if not
        local_dir = os.path.join(base_dir,'{0}.{1}.{2}'.format(YY,MM,DD))
        os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
        #-- print original and symbolic link of file
        if VERBOSE:
            print('{0} -->'.format(os.path.join(scf_outgoing,f)))
            print('\t{0}'.format(os.path.join(local_dir,f)))
        #-- create symbolic link of file from scf_outgoing to local
        os.symlink(os.path.join(scf_outgoing,f), os.path.join(local_dir,f))

#-- run main program
if __name__ == '__main__':
    main()
