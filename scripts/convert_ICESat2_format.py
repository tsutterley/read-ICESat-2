#!/usr/bin/env python
u"""
convert_ICESat2_format.py
Written by Tyler Sutterley (12/2022)

Converts ICESat-2 HDF5 datafiles to zarr or rechunked HDF5 datafiles

zarr files make large datasets easily accessible to distributed computing on
    both local filesystems and cloud-based object stores
    - arrays are divided into chunks and compressed
    - metadata are stored in lightweight .json files

rechunked HDF5 files can be more optimized for cloud-based object stores

CALLING SEQUENCE:
    python convert_ICESat2_format.py --release=003 ATL06

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
    -D X, --directory X: working data directory
    -Y X, --year X: years to run separated by commas
    -S X, --subdirectory X: subdirectories to run separated by commas
    -r X, --release X: ICESat-2 data release to run
    -v X, --version X: ICESat-2 data version to run
    -t X, --track X: ICESat-2 reference ground tracks to run
    -g X, --granule X: ICESat-2 granule regions to run
    -f X, --format X: output file format (zarr, HDF5)
    -c X, --chunks X: Rechunk output files to size
    -P X, --np X: Number of processes to use in file conversion
    -C, --clobber: Overwrite existing files
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Local permissions mode of the converted files

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

UPDATE HISTORY:
    Updated 12/2022: single implicit import of altimetry tools
    Updated 06/2022: use explicit import of convert functions
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
        added parsing for converting file lines to arguments
    Updated 07/2021: set context for multiprocessing to fork child processes
    Updated 01/2021: generalized to output either zarr or rechunked HDF5
    Updated 10/2020: using argparse to set parameters. added verbose keyword
        added chunks keyword to rechunk output zarr files
        using convert module to convert from HDF5 to zarr
    Written 06/2020
"""
from __future__ import print_function

import sys
import os
import re
import logging
import argparse
import traceback
import multiprocessing as mp
import icesat2_toolkit as is2tk

# PURPOSE: convert the ICESat-2 elevation data from HDF5 to zarr
# or rechunked HDF5 formats
def convert_ICESat2_format(DIRECTORY, PRODUCTS, RELEASE, VERSIONS, GRANULES,
    TRACKS, YEARS=None, SUBDIRECTORY=None, FORMAT=None, CHUNKS=None,
    PROCESSES=0, CLOBBER=False, VERBOSE=False, MODE=0o775):

    # create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # regular expression operator for finding files of a particular granule
    # find ICESat-2 HDF5 files in the subdirectory for product and release
    regex_track = '|'.join([rf'{T:04d}' for T in TRACKS])
    regex_granule = '|'.join([rf'{G:02d}' for G in GRANULES])
    regex_version = '|'.join([rf'{V:02d}' for V in VERSIONS])
    file_regex_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})(\d{{2}})({2})_({3})_({4})(.*?).(h5)$')

    # regular expression operator for finding subdirectories
    if SUBDIRECTORY:
        # convert particular subdirectories for product
        R2 = re.compile(r'('+r'|'.join(SUBDIRECTORY)+r')', re.VERBOSE)
    elif YEARS:
        # convert particular years for product
        regex_pattern = '|'.join(rf'{y:d}' for y in YEARS)
        R2 = re.compile(rf'({regex_pattern}).(\d+).(\d+)', re.VERBOSE)
    else:
        # convert all available subdirectories for product
        R2 = re.compile(r'(\d+).(\d+).(\d+)', re.VERBOSE)

    # build list of HDF5 files
    hdf5_file_list = []
    # for each ICESat-2 product listed
    for p in PRODUCTS:
        logging.info(f'PRODUCT={p}')
        # local file directory
        ddir = os.path.join(DIRECTORY,f'{p}.{RELEASE}')
        subdirectories = [sd for sd in os.listdir(ddir) if R2.match(sd)]
        # compile regular expression operator for product, release and version
        args = (p,regex_track,regex_granule,RELEASE,regex_version)
        R1 = re.compile(file_regex_pattern.format(*args), re.VERBOSE)
        # for each subdirectory
        for sd in subdirectories:
            # find matching files (for granule, release, version, track)
            # add to list of HDF5 files
            hdf5_file_list.extend([os.path.join(ddir,sd,f) for f in
                os.listdir(os.path.join(ddir,sd)) if R1.match(f)])

    # convert in series if PROCESSES = 0
    if (PROCESSES == 0):
        # convert each ICESat-2 data file
        for f in hdf5_file_list:
            # convert ICESat-2 file to output format
            output = convert_HDF5(f,FORMAT=FORMAT,CHUNKS=CHUNKS,
                CLOBBER=CLOBBER,MODE=MODE)
            # print the output string
            logging.info(output)
    else:
        # set multiprocessing start method
        ctx = mp.get_context("fork")
        # convert in parallel with multiprocessing Pool
        pool = ctx.Pool(processes=PROCESSES)
        # convert each ICESat-2 data file
        output = []
        for hdf5_file in hdf5_file_list:
            # convert ICESat-2 file to output format
            kwds = dict(FORMAT=FORMAT, CHUNKS=CHUNKS, CLOBBER=CLOBBER,
                MODE=MODE)
            output.append(pool.apply_async(multiprocess_convert,
                args=(hdf5_file,),kwds=kwds))
        # start multiprocessing jobs
        # close the pool
        # prevents more tasks from being submitted to the pool
        pool.close()
        # exit the completed processes
        pool.join()
        # print the output string
        for out in output:
            logging.info(out.get())

# PURPOSE: wrapper for running conversion program in multiprocessing mode
def multiprocess_convert(hdf5_file, FORMAT=None, CHUNKS=None, CLOBBER=False,
    MODE=0o775):
    try:
        output = convert_HDF5(hdf5_file,FORMAT=FORMAT,CHUNKS=CHUNKS,
            CLOBBER=CLOBBER,MODE=MODE)
    except Exception as exc:
        # if there has been an error exception
        # print the type, value, and stack trace of the
        # current exception being handled
        logging.critical(f'process id {os.getpid():d} failed')
        logging.error(traceback.format_exc())
    else:
        return output

# PURPOSE: convert the HDF5 file and change permissions
def convert_HDF5(hdf5_file,FORMAT=None,CHUNKS=None,CLOBBER=False,MODE=0o775):
    # split extension from input HDF5 file
    fileBasename,fileExtension = os.path.splitext(hdf5_file)
    # convert HDF5 file into output format
    if (FORMAT == 'zarr'):
        output_file = f'{fileBasename}.zarr'
    elif (FORMAT == 'HDF5'):
        output_file = f'{fileBasename}.h5'
    # if output file exists in file system: check if HDF5 file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    # last modification time of HDF5 file
    hdf5_mtime = os.stat(hdf5_file).st_mtime
    # check if local version of file exists
    if os.access(output_file, os.F_OK):
        # check last modification time of output file
        zarr_mtime = os.stat(output_file).st_mtime
        # if HDF5 file is newer: overwrite the output file
        if (hdf5_mtime > zarr_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'

    # if file does not exist, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        # output string for printing files transferred
        output = '{0} -->\n\t{1}{2}\n'.format(hdf5_file,output_file,OVERWRITE)
        # copy everything from the HDF5 file to the output file
        conv = is2tk.convert(filename=hdf5_file, reformat=FORMAT)
        conv.file_converter(chunks=CHUNKS)
        # keep remote modification time of file and local access time
        os.utime(output_file, (os.stat(output_file).st_atime, hdf5_mtime))
        os.chmod(output_file, MODE)
        # return the output string
        return output

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts ICESat-2 HDF5 datafiles to zarr or
            rechunked HDF5 datafiles
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # ICESat-2 Products
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
    # command line parameters
    parser.add_argument('products',
        metavar='PRODUCTS', type=str, nargs='+',
        choices=PRODUCTS.keys(),
        help='ICESat-2 products to convert')
    # working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
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
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 Data Version')
    # ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
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
    parser.add_argument('--chunks','-c',
        type=int,
        help='Rechunk output files to size')
    # run conversion in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file conversion')
    # clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
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

    # convert HDF5 files for each data product
    convert_ICESat2_format(args.directory, args.products, args.release,
        args.version, args.granule, args.track, YEARS=args.year,
        SUBDIRECTORY=args.subdirectory, FORMAT=args.format,
        CHUNKS=args.chunks, PROCESSES=args.np, CLOBBER=args.clobber,
        VERBOSE=args.verbose, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
