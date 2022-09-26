#!/usr/bin/env python
u"""
test_download_and_read.py (08/2022)

UPDATE HISTORY:
    Updated 08/2022: added backup AWS access to files
    Updated 02/2022: add CMR query tests for cycles, tracks and granules
    Updated 01/2022: update data releases for all products
    Updated 03/2021: modify ATL11 to read from explicitly named groups
    Updated 11/2020: use output error string returned by from_nsidc
    Written 08/2020
"""
import warnings
import pytest
import posixpath
import icesat2_toolkit.utilities
from icesat2_toolkit.read_ICESat2_ATL03 import read_HDF5_ATL03
from icesat2_toolkit.read_ICESat2_ATL06 import read_HDF5_ATL06
from icesat2_toolkit.read_ICESat2_ATL07 import read_HDF5_ATL07
from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11
from icesat2_toolkit.read_ICESat2_ATL12 import read_HDF5_ATL12

#-- PURPOSE: Download an ATL03 file from NSIDC and check that read program runs
def test_ATL03_download_and_read(username,password):
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL03.005','2018.10.14',
        'ATL03_20181014000347_02350101_005_01.h5']
    buffer, error = icesat2_toolkit.utilities.from_nsidc(HOST,
        username=username, password=password,
        local=HOST[-1], verbose=True)
    #-- attempt download from AWS
    if not buffer:
        url = posixpath.dirname(icesat2_toolkit.utilities._s3_endpoints['nsidc'])
        bucket = icesat2_toolkit.utilities._s3_buckets['nsidc']
        key = icesat2_toolkit.utilities.s3_key(posixpath.join(*HOST))
        buffer, error = icesat2_toolkit.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=HOST[-1],
            verbose=True)
    #-- raise exception if download error
    if not buffer:
        raise Exception(error)
    #-- read ATL03 data from downloaded HDF5 file
    IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams = read_HDF5_ATL03(HOST[-1],
        ATTRIBUTES=False, VERBOSE=True)
    assert all(gtx in IS2_atl03_mds.keys() for gtx in IS2_atl03_beams)

#-- PURPOSE: Download an ATL06 file from NSIDC and check that read program runs
def test_ATL06_download_and_read(username,password):
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL06.005','2018.10.14',
        'ATL06_20181014001049_02350102_005_01.h5']
    buffer, error = icesat2_toolkit.utilities.from_nsidc(HOST,
        username=username, password=password,
        local=HOST[-1], verbose=True)
    #-- attempt download from AWS
    if not buffer:
        url = posixpath.dirname(icesat2_toolkit.utilities._s3_endpoints['nsidc'])
        bucket = icesat2_toolkit.utilities._s3_buckets['nsidc']
        key = icesat2_toolkit.utilities.s3_key(posixpath.join(*HOST))
        buffer, error = icesat2_toolkit.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=HOST[-1],
            verbose=True)
    #-- raise exception if download error
    if not buffer:
        raise Exception(error)
    #-- read ATL06 data from downloaded HDF5 file
    IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_beams = read_HDF5_ATL06(HOST[-1],
        ATTRIBUTES=False, HISTOGRAM=False, QUALITY=False, VERBOSE=True)
    assert all(gtx in IS2_atl06_mds.keys() for gtx in IS2_atl06_beams)

#-- PURPOSE: Download an ATL07 file from NSIDC and check that read program runs
def test_ATL07_download_and_read(username,password):
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL07.005','2018.10.14',
        'ATL07-01_20181014000347_02350101_005_03.h5']
    buffer, error = icesat2_toolkit.utilities.from_nsidc(HOST,
        username=username, password=password,
        local=HOST[-1], verbose=True)
    #-- attempt download from AWS
    if not buffer:
        url = posixpath.dirname(icesat2_toolkit.utilities._s3_endpoints['nsidc'])
        bucket = icesat2_toolkit.utilities._s3_buckets['nsidc']
        key = icesat2_toolkit.utilities.s3_key(posixpath.join(*HOST))
        buffer, error = icesat2_toolkit.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=HOST[-1],
            verbose=True)
    #-- raise exception if download error
    if not buffer:
        raise Exception(error)
    #-- read ATL07 data from downloaded HDF5 file
    IS2_ATL07_mds,IS2_ATL07_attrs,IS2_ATL07_beams = read_HDF5_ATL07(HOST[-1],
        ATTRIBUTES=False, VERBOSE=True)
    assert all(gtx in IS2_ATL07_mds.keys() for gtx in IS2_ATL07_beams)

#-- PURPOSE: Download an ATL11 file from NSIDC and check that read program runs
def test_ATL11_download_and_read(username,password):
    #-- attempt download from on-prem NSIDC
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL11.004','2019.03.29',
        'ATL11_000103_0311_004_01.h5']
    buffer, error = icesat2_toolkit.utilities.from_nsidc(HOST,
        username=username, password=password,
        local=HOST[-1], verbose=True)
    #-- attempt download from AWS
    if not buffer:
        url = posixpath.dirname(icesat2_toolkit.utilities._s3_endpoints['nsidc'])
        bucket = icesat2_toolkit.utilities._s3_buckets['nsidc']
        key = icesat2_toolkit.utilities.s3_key(posixpath.join(*HOST))
        buffer, error = icesat2_toolkit.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=HOST[-1],
            verbose=True)
    #-- raise exception if download error
    if not buffer:
        raise Exception(error)
    #-- read ATL12 data from downloaded HDF5 file
    GROUPS = ['cycle_stats','ref_surf','crossing_track_data']
    IS2_ATL11_mds,IS2_ATL11_attrs,IS2_ATL11_pairs = read_HDF5_ATL11(HOST[-1],
        ATTRIBUTES=False, GROUPS=GROUPS, VERBOSE=True)
    assert all(ptx in IS2_ATL11_mds.keys() for ptx in IS2_ATL11_pairs)
    ptx = IS2_ATL11_pairs[0]
    assert all(group in IS2_ATL11_mds[ptx].keys() for group in GROUPS)

#-- PURPOSE: Download an ATL12 file from NSIDC and check that read program runs
def test_ATL12_download_and_read(username,password):
    HOST = ['https://n5eil01u.ecs.nsidc.org','ATLAS','ATL12.004','2018.10.14',
        'ATL12_20181014031222_02370101_004_01.h5']
    buffer, error = icesat2_toolkit.utilities.from_nsidc(HOST,
        username=username, password=password,
        local=HOST[-1], verbose=True)
    #-- attempt download from AWS
    if not buffer:
        url = posixpath.dirname(icesat2_toolkit.utilities._s3_endpoints['nsidc'])
        bucket = icesat2_toolkit.utilities._s3_buckets['nsidc']
        key = icesat2_toolkit.utilities.s3_key(posixpath.join(*HOST))
        buffer, error = icesat2_toolkit.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=HOST[-1],
            verbose=True)
    #-- raise exception if download error
    if not buffer:
        raise Exception(error)
    #-- read ATL12 data from downloaded HDF5 file
    IS2_ATL12_mds,IS2_ATL12_attrs,IS2_ATL12_beams = read_HDF5_ATL12(HOST[-1],
        ATTRIBUTES=False, VERBOSE=True)
    assert all(gtx in IS2_ATL12_mds.keys() for gtx in IS2_ATL12_beams)

#-- PURPOSE: test CMR queries for specific cycles
def test_cmr_query_cycles():
    ids,urls = icesat2_toolkit.utilities.cmr(product='ATL06',
        release='005',cycles=[2,3],tracks=752,granules=10,
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_005_01.h5',
        'ATL06_20190517125119_07520310_005_01.h5']
    assert all(id in valid for id in ids)

#-- PURPOSE: test CMR queries for specific tracks
def test_cmr_query_tracks():
    ids,urls = icesat2_toolkit.utilities.cmr(product='ATL06',
        release='005',cycles=2,tracks=[752,753],granules=10,
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_005_01.h5',
        'ATL06_20190215184558_07530210_005_01.h5']
    assert all(id in valid for id in ids)

#-- PURPOSE: test CMR queries for specific granules
def test_cmr_query_granules():
    ids,urls = icesat2_toolkit.utilities.cmr(product='ATL06',
        release='005',cycles=2,tracks=752,granules=[10,11,12],
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_005_01.h5',
        'ATL06_20190215171921_07520211_005_01.h5',
        'ATL06_20190215172504_07520212_005_01.h5']
    assert all(id in valid for id in ids)
