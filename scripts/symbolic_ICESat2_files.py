#!/usr/bin/env python
u"""
symbolic_ICESat2_files.py
Written by Tyler Sutterley (10/2021)
Creates symbolic links for ICESat-2 HDF5 files organized by date

CALLING SEQUENCE:
    python symbolic_ICESat2_files.py --product ATL06 --release 003 \
        --granule 10 11 12 --cycle 1 2 --directory <path_to_directory>
        --scf_outgoing <path_to_outgoing> --verbose --mode 0o775

COMMAND LINE OPTIONS:
    -h, --help: list the command line options
    -D X, --directory X: local working directory for creating symbolic links
    --product X: ICESat-2 data product to create symbolic links
    --release X: ICESat-2 data release to create symbolic links
    --version X: ICESat-2 data version to create symbolic links
    --granule X: ICESat-2 granule regions to create symbolic links
    --cycle X: ICESat-2 cycle to create symbolic links
    --track X: ICESat-2 tracks to create symbolic links
    --scf_incoming X: directory on the SCF where the rscf sends PANS
    --scf_outgoing X: directory on the SCF where the data resides
    -V, --verbose: output information about each symbolic link
    -M X, --mode X: permission mode of directories

UPDATE HISTORY:
    Updated 10/2021: using python logging for handling verbose output
    Updated 11/2020: add exception for FileExistsError to skip files
    Updated 10/2020: using argparse to set parameters
    Updated 05/2020: adjust regular expression to run ATL07 sea ice products
    Written 07/2019
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import numpy as np

#-- Main program that calls symbolic_ICESat2_files()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Creates symbolic links for ICESat-2 HDF5 files from the
            local scf directory to a separate directory organized by date
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
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory for symbolic link')
    #-- ICESat-2 parameters
    #-- ICESat-2 data product
    parser.add_argument('--product','-p',
        metavar='PRODUCTS', type=str,
        choices=PRODUCTS.keys(), default='ATL06',
        help='ICESat-2 data product to create symbolic links')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='004',
        help='ICESat-2 data release to create symbolic links')
    #-- ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 data versions to create symbolic links')
    #-- ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 granule regions to create symbolic links')
    #-- ICESat-2 orbital cycle
    parser.add_argument('--cycle','-c',
        type=int, nargs='+',
        default=range(1,10),
        help='ICESat-2 orbital cycles to create symbolic links')
    #-- ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs) to create symbolic links')
    #-- ICESat-2 Science Computing Facility (SCF) parameters
    parser.add_argument('--scf_incoming',
        type=str,
        help='Directory on the SCF where the rscf sends PANS')
    parser.add_argument('--scf_outgoing',
        type=str,
        help='Directory on the SCF where the data resides')
    #-- verbose will output information about each symbolic link
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Output information about each symbolic link')
    #-- permissions mode of the local directories (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output directories')
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevel = logging.INFO if args.verbose else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- run program
    symbolic_ICESat2_files(args.directory, args.scf_incoming, args.scf_outgoing,
        args.product, args.release, args.version, args.granule, args.cycle,
        args.track, MODE=args.mode)

#-- PURPOSE: copy ICESat-2 files to data directory with data subdirectories
def symbolic_ICESat2_files(base_dir, scf_incoming, scf_outgoing, PRODUCT,
    RELEASE, VERSIONS, GRANULES, CYCLES, TRACKS, MODE=0o775):
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
        #-- attempt to create the symbolic link else continue
        try:
            #-- create symbolic link of file from scf_outgoing to local
            os.symlink(os.path.join(scf_outgoing,f), os.path.join(local_dir,f))
        except FileExistsError:
            continue
        else:
            #-- print original and symbolic link of file
            args = (os.path.join(scf_outgoing,f),os.path.join(local_dir,f))
            logging.info('{0} -->\n\t{1}'.format(*args))

#-- run main program
if __name__ == '__main__':
    main()
