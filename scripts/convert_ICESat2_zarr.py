#!/usr/bin/env python
u"""
convert_ICESat2_zarr.py
Written by Tyler Sutterley (06/2020)

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
    -D X, --directory: working data directory
    -Y X, --year=X: years to run separated by commas
    -S X, --subdirectory=X: subdirectories to run separated by commas
    --release=X: ICESat-2 data release to run
    --version=X: ICESat-2 data version to run
    --track=X: ICESat-2 reference ground tracks to run
    --granule=X: ICESat-2 granule regions to run
    -P X, --np=X: Number of processes to use in file conversion
    -M X, --mode=X: Local permissions mode of the converted files
    -C, --clobber: Overwrite existing zarr files

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
    Written 06/2020
"""
from __future__ import print_function

import sys
import os
import re
import zarr
import h5py
import getopt
import posixpath
import traceback
import itertools
import numpy as np
import calendar, time
import multiprocessing as mp

#-- PURPOSE: convert the ICESat-2 elevation data from HDF5 to zarr
def convert_ICESat2_zarr(DIRECTORY,PRODUCTS,RELEASE,VERSIONS,GRANULES,TRACKS,
    YEARS=None,SUBDIRECTORY=None,PROCESSES=0,MODE=0o775,CLOBBER=False):

    #-- regular expression operator for finding files of a particular granule
    #-- find ICESat-2 HDF5 files in the subdirectory for product and release
    regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACKS])
    regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULES])
    regex_version = '|'.join(['{0:02d}'.format(V) for V in VERSIONS])
    file_regex_pattern = ('{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})(\d{{2}})'
        '(\d{{2}})(\d{{2}})_({1})(\d{{2}})({2})_({3})_({4})(.*?).(h5)$')

    #-- regular expression operator for finding subdirectories
    if SUBDIRECTORY:
        #-- convert particular subdirectories for product
        R2 = re.compile('('+'|'.join(SUBDIRECTORY)+')', re.VERBOSE)
    elif YEARS:
        #-- convert particular years for product
        regex_pattern = '|'.join('{0:d}'.format(y) for y in YEARS)
        R2 = re.compile('({0}).(\d+).(\d+)'.format(regex_pattern), re.VERBOSE)
    else:
        #-- convert all available subdirectories for product
        R2 = re.compile('(\d+).(\d+).(\d+)', re.VERBOSE)

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
        for hdf5_file in hdf5_file_list:
            #-- convert ICESat-2 file to zarr
            output = HDF5_to_zarr(hdf5_file,CLOBBER=CLOBBER,MODE=MODE)
            #-- print the output string
            print(output)
    else:
        #-- convert in parallel with multiprocessing Pool
        pool = mp.Pool(processes=PROCESSES)
        #-- convert each ICESat-2 data file
        output = []
        for hdf5_file in hdf5_file_list:
            #-- convert ICESat-2 file to zarr
            kwds = dict(CLOBBER=CLOBBER, MODE=MODE)
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
            print(out.get())

#-- PURPOSE: wrapper for running conversion program in multiprocessing mode
def multiprocess_convert(hdf5_file, CLOBBER=False, MODE=0o775):
    try:
        output = HDF5_to_zarr(hdf5_file,CLOBBER=CLOBBER,MODE=MODE)
    except:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        print('process id {0:d} failed'.format(os.getpid()))
        traceback.print_exc()
    else:
        return output

#-- PURPOSE: convert the HDF5 file to zarr and change permissions
def HDF5_to_zarr(hdf5_file,CLOBBER=False,MODE=0o775):
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
                copy_from_HDF5(source[k], dest, name=k)
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
            raise CopyError('an object {!r} already exists in destination '
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
            raise CopyError('an array {!r} already exists in destination '
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

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -Y X, --year=X\t\tYears to run separated by commas')
    print(' -S X, --subdirectory=X\tSubdirectories to run separated by commas')
    print(' --release=X\t\tICESat-2 data release to run')
    print(' --version=X\t\tICESat-2 data version to run')
    print(' --granule=X\t\tICESat-2 granule regions to run')
    print(' --track=X\t\tICESat-2 reference ground tracks to run')
    print(' -P X, --np=X\t\tNumber of processes to use in file conversion')
    print(' -M X, --mode=X\t\tPermission mode of the converted files')
    print(' -C, --clobber\t\tOverwrite existing zarr files\n')

#-- Main program that calls convert_ICESat2_zarr()
def main():
    #-- Read the system arguments listed after the program
    short_options = 'hD:Y:S:P:M:C'
    long_options=['help','directory=','year=','subdirectory=','release=',
        'version=','granule=','track=','np=','mode=','clobber']
    optlist,arglist=getopt.getopt(sys.argv[1:],short_options,long_options)

    #-- command line parameters
    #-- Working data directory
    DIRECTORY = os.getcwd()
    YEARS = None
    SUBDIRECTORY = None
    VERSIONS = np.arange(1,10)
    RELEASE = '003'
    GRANULES = np.arange(1,15)
    TRACKS = np.arange(1,1388)
    #-- run conversion in series if processes is 0
    PROCESSES = 0
    #-- permissions mode of the converted files (number in octal)
    MODE = 0o775
    CLOBBER = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-Y","--year"):
            YEARS = np.array(arg.split(','), dtype=np.int)
        elif opt in ("-S","--subdirectory"):
            SUBDIRECTORY = arg.split(',')
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--release",):
            RELEASE = '{0:03d}'.format(int(arg))
        elif opt in ("--version",):
            VERSIONS = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--granule",):
            GRANULES = np.array(arg.split(','), dtype=np.int)
        elif opt in ("--track",):
            TRACKS = np.sort(arg.split(',')).astype(np.int)
        elif opt in ("-P","--np"):
            PROCESSES = int(arg)
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-C","--clobber"):
            CLOBBER = True

    #-- Pre-ICESat-2 and IceBridge Products
    PROD = {}
    PROD['ATL03'] = 'Global Geolocated Photon Data'
    PROD['ATL04'] = 'Normalized Relative Backscatter'
    PROD['ATL06'] = 'Land Ice Height'
    PROD['ATL07'] = 'Sea Ice Height'
    PROD['ATL08'] = 'Land and Vegetation Height'
    PROD['ATL09'] = 'Atmospheric Layer Characteristics'
    PROD['ATL10'] = 'Sea Ice Freeboard'
    PROD['ATL12'] = 'Ocean Surface Height'
    PROD['ATL13'] = 'Inland Water Surface Height'

    #-- enter dataset to transfer as system argument
    if not arglist:
        for key,val in PROD.items():
            print('{0}: {1}'.format(key, val))
        raise Exception('No System Arguments Listed')

    #-- check that each data product entered was correctly typed
    keys = ','.join(sorted([key for key in PROD.keys()]))
    for p in arglist:
        if p not in PROD.keys():
            raise IOError('Incorrect Data Product Entered ({0})'.format(keys))

    #-- convert HDF5 files to zarr files for each data product
    convert_ICESat2_zarr(DIRECTORY, arglist, RELEASE, VERSIONS,
        GRANULES, TRACKS, YEARS=YEARS, SUBDIRECTORY=SUBDIRECTORY,
        PROCESSES=PROCESSES, MODE=MODE, CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
    main()
