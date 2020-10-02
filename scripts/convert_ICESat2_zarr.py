#!/usr/bin/env python
u"""
convert_ICESat2_zarr.py
Written by Tyler Sutterley (10/2020)

Converts ICESat-2 HDF5 datafiles to zarr datafiles

zarr files make large datasets easily accessible to distributed computing on
    both local filesystems and cloud-based object stores
    - arrays are divided into chunks and compressed
    - metadata are stored in lightweight .json files

CALLING SEQUENCE:
    python convert_ICESat2_zarr.py --release=001 ATL06

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
    --chunks X: Rechunk zarr files to size
    -P X, --np X: Number of processes to use in file conversion
    -C, --clobber: Overwrite existing zarr files
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

UPDATE HISTORY:
    Updated 10/2020: using argparse to set parameters. added verbose keyword
        added chunks keyword to rechunk output zarr files
    Written 06/2020
"""
from __future__ import print_function

import sys
import os
import re
import zarr
import h5py
import argparse
import posixpath
import traceback
import itertools
import numpy as np
import calendar, time
import multiprocessing as mp

#-- PURPOSE: convert the ICESat-2 elevation data from HDF5 to zarr
def convert_ICESat2_zarr(DIRECTORY, PRODUCTS, RELEASE, VERSIONS, GRANULES,
    TRACKS, YEARS=None, SUBDIRECTORY=None, CHUNKS=None, PROCESSES=0,
    CLOBBER=False, VERBOSE=False, MODE=0o775):

    #-- regular expression operator for finding files of a particular granule
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    file_regex_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})_({1})(\d{{2}})({2})_({3})_({4})(.*?).(h5)$')

    #-- regular expression operator for finding subdirectories
    if SUBDIRECTORY:
        #-- convert particular subdirectories for product
        R2 = re.compile(r'('+'|'.join(SUBDIRECTORY)+')', re.VERBOSE)
    elif YEARS:
        #-- convert particular years for product
        regex_pattern = '|'.join('{0:d}'.format(y) for y in YEARS)
        R2 = re.compile(r'({0}).(\d+).(\d+)'.format(regex_pattern), re.VERBOSE)
    else:
        #-- convert all available subdirectories for product
        R2 = re.compile(r'(\d+).(\d+).(\d+)', re.VERBOSE)

    #-- build list of HDF5 files
    hdf5_file_list = []
    #-- for each ICESat-2 product listed
    for p in PRODUCTS:
        print('PRODUCT={0}'.format(p))
        #-- local file directory
        ddir = os.path.join(DIRECTORY,'{0}.{1}'.format(p,RELEASE))
        subdirectories = [sd for sd in os.listdir(ddir) if R2.match(sd)]
        #-- compile regular expression operator for product, release and version
        args = (p,regex_track,regex_granule,RELEASE,regex_version)
        R1 = re.compile(file_regex_pattern.format(*args), re.VERBOSE)
        #-- for each subdirectory
        for sd in subdirectories:
            #-- find matching files (for granule, release, version, track)
            #-- add to list of HDF5 files
            hdf5_file_list.extend([os.path.join(ddir,sd,f) for f in
                os.listdir(os.path.join(ddir,sd)) if R1.match(f)])

    #-- convert in series if PROCESSES = 0
    if (PROCESSES == 0):
        #-- convert each ICESat-2 data file
        for f in hdf5_file_list:
            #-- convert ICESat-2 file to zarr
            output = HDF5_to_zarr(f,CHUNKS=CHUNKS,CLOBBER=CLOBBER,MODE=MODE)
            #-- print the output string
            print(output) if VERBOSE else None
    else:
        #-- convert in parallel with multiprocessing Pool
        pool = mp.Pool(processes=PROCESSES)
        #-- convert each ICESat-2 data file
        output = []
        for hdf5_file in hdf5_file_list:
            #-- convert ICESat-2 file to zarr
            kwds = dict(CHUNKS=CHUNKS, CLOBBER=CLOBBER, MODE=MODE)
            output.append(pool.apply_async(multiprocess_convert,
                args=(hdf5_file,),kwds=kwds))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()
        #-- print the output string
        for out in output:
            print(out.get()) if VERBOSE else None

#-- PURPOSE: wrapper for running conversion program in multiprocessing mode
def multiprocess_convert(hdf5_file, CHUNKS=None, CLOBBER=False, MODE=0o775):
    try:
        output = HDF5_to_zarr(hdf5_file,CHUNKS=CHUNKS,CLOBBER=CLOBBER,MODE=MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
    else:
        return output

#-- PURPOSE: convert the HDF5 file to zarr and change permissions
def HDF5_to_zarr(hdf5_file,CHUNKS=None,CLOBBER=False,MODE=0o775):
    #-- split extension from input HDF5 file
    fileBasename,fileExtension = os.path.splitext(hdf5_file)
    #-- convert HDF5 file into zarr file
    zarr_file = '{0}.zarr'.format(fileBasename)
    #-- if zarr file exists in file system: check if HDF5 file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- last modification time of HDF5 file
    hdf5_mtime = os.stat(hdf5_file).st_mtime
    #-- check if local version of file exists
    if os.access(zarr_file, os.F_OK):
        #-- check last modification time of zarr file
        zarr_mtime = os.stat(zarr_file).st_mtime
        #-- if HDF5 file is newer: overwrite the zarr file
        if (hdf5_mtime > zarr_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'

    #-- if zarr file does not exist, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- output string for printing files transferred
        output = '{0} -->\n\t{1}{2}\n'.format(hdf5_file,zarr_file,OVERWRITE)
        #-- copy everything from the HDF5 file to the zarr file
        with h5py.File(hdf5_file,mode='r') as source:
            dest = zarr.open_group(zarr_file,mode='w')
            #-- value checks on output zarr
            if not hasattr(dest, 'create_dataset'):
                raise ValueError('dest must be a group, got {!r}'.format(dest))
            #-- for each key in the root of the hdf5 file structure
            for k in source.keys():
                copy_from_HDF5(source[k], dest, name=k, chunks=CHUNKS)
        #-- keep remote modification time of file and local access time
        os.utime(zarr_file, (os.stat(zarr_file).st_atime, hdf5_mtime))
        os.chmod(zarr_file, MODE)
        #-- return the output string
        return output

#-- PURPOSE: Copy a named variable from the HDF5 file to the zarr file
def copy_from_HDF5(source, dest, name, **create_kws):
    """Copy a named variable from the `source` HDF5 into the `dest` zarr"""
    if hasattr(source, 'shape'):
        #-- copy a dataset/array
        if dest is not None and name in dest:
            raise Exception('an object {!r} already exists in destination '
                '{!r}'.format(name, dest.name))
        #-- setup creation keyword arguments
        kws = create_kws.copy()
        #-- setup chunks option, preserve by default
        kws.setdefault('chunks', source.chunks)
        #-- setup compression options
        #-- from h5py to zarr: use zarr default compression options
        kws.setdefault('fill_value', source.fillvalue)
        #-- create new dataset in destination
        ds=dest.create_dataset(name,shape=source.shape,dtype=source.dtype,**kws)
        #-- copy data going chunk by chunk to avoid loading in entirety
        shape = ds.shape
        chunks = ds.chunks
        chunk_offsets = [range(0, s, c) for s, c in zip(shape, chunks)]
        for offset in itertools.product(*chunk_offsets):
            sel = tuple(slice(o, min(s, o + c)) for o, s, c in
                zip(offset, shape, chunks))
            ds[sel] = source[sel]
        #-- copy attributes
        attrs = {key:attributes_encoder(source.attrs[key]) for key in
            source.attrs.keys() if attributes_encoder(source.attrs[key])}
        ds.attrs.update(attrs)
    else:
        #-- copy a group
        if (dest is not None and name in dest and hasattr(dest[name], 'shape')):
            raise Exception('an array {!r} already exists in destination '
                '{!r}'.format(name, dest.name))
        #-- require group in destination
        grp = dest.require_group(name)
        #-- copy attributes
        attrs = {key:attributes_encoder(source.attrs[key]) for key in
            source.attrs.keys() if attributes_encoder(source.attrs[key])}
        grp.attrs.update(attrs)
        #-- recursively copy from source
        for k in source.keys():
            copy_from_HDF5(source[k], grp, name=k)

#-- PURPOSE: encoder for copying the file attributes
def attributes_encoder(attr):
    """Custom encoder for copying file attributes in Python 3"""
    if isinstance(attr, (bytes, bytearray)):
        return attr.decode('utf-8')
    if isinstance(attr, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32,
        np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
        return int(attr)
    elif isinstance(attr, (np.float_, np.float16, np.float32, np.float64)):
        return float(attr)
    elif isinstance(attr, (np.ndarray)):
        if not isinstance(attr[0], (object)):
            return attr.tolist()
    elif isinstance(attr, (np.bool_)):
        return bool(attr)
    elif isinstance(attr, (np.void)):
        return None
    else:
        return attr

#-- Main program that calls convert_ICESat2_zarr()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Converts ICESat-2 HDF5 datafiles to zarr datafiles
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
    parser.add_argument('products',
        metavar='PRODUCTS', type=str, nargs='+',
        choices=PRODUCTS.keys(),
        help='ICESat-2 products to convert to zarr')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=os.path.expanduser, default=os.getcwd(),
        help='Working data directory')
    #-- years of data to run
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to run')
    #-- subdirectories of data to run
    parser.add_argument('--subdirectory','-S',
        type=str, nargs='+',
        help='subdirectories of data to run')
    #-- ICESat-2 data release
    parser.add_argument('--release','-r',
        type=str, default='003',
        help='ICESat-2 Data Release')
    #-- ICESat-2 data version
    parser.add_argument('--version','-v',
        type=int, nargs='+', default=range(1,10),
        help='ICESat-2 Data Version')
    #-- ICESat-2 granule region
    parser.add_argument('--granule','-g',
        metavar='REGION', type=int, nargs='+',
        choices=range(1,15), default=range(1,15),
        help='ICESat-2 Granule Region')
    #-- ICESat-2 reference ground tracks
    parser.add_argument('--track','-t',
        metavar='RGT', type=int, nargs='+',
        choices=range(1,1388), default=range(1,1388),
        help='ICESat-2 Reference Ground Tracks (RGTs)')
    #-- rechunk output zarr data
    parser.add_argument('--chunks',
        type=int,
        help='Rechunk zarr files to size')
    #-- run conversion in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file conversion')
    #-- clobber will overwrite the existing data
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data')
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the converted files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- convert HDF5 files to zarr files for each data product
    convert_ICESat2_zarr(args.directory, args.products, args.release,
        args.version, args.granule, args.track, YEARS=args.year,
        SUBDIRECTORY=args.subdirectory, CHUNKS=args.chunks, PROCESSES=args.np,
        CLOBBER=args.clobber, VERBOSE=args.verbose, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
