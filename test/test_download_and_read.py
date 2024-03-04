#!/usr/bin/env python
u"""
test_download_and_read.py (03/2024)

UPDATE HISTORY:
    Updated 03/2024: use CMR queries to find granules in each download test
    Updated 03/2023: replace deprecated functions with new io functions
    Updated 08/2022: added backup AWS access to files
    Updated 02/2022: add CMR query tests for cycles, tracks and granules
    Updated 01/2022: update data releases for all products
    Updated 03/2021: modify ATL11 to read from explicitly named groups
    Updated 11/2020: use output error string returned by from_nsidc
    Written 08/2020
"""
import pytest
import posixpath
import icesat2_toolkit as is2tk

# PURPOSE: Download an ATL03 file from NSIDC and check that read program runs
def test_ATL03_download_and_read(username, password):
    # query CMR for ATL03 data
    ids, urls = is2tk.utilities.cmr(product='ATL03',
        release='006', cycles=1, tracks=235, granules=1,
        verbose=False)
    # attempt download from on-prem NSIDC
    buffer, error = is2tk.utilities.from_nsidc(urls[0],
        username=username, password=password,
        local=ids[0], verbose=True)
    # attempt download from AWS
    if not buffer:
        url = posixpath.dirname(is2tk.utilities._s3_endpoints['nsidc'])
        bucket = is2tk.utilities._s3_buckets['nsidc']
        key = is2tk.utilities.s3_key(urls[0])
        buffer, error = is2tk.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=ids[0],
            verbose=True)
    # raise exception if download error
    if not buffer:
        raise Exception(error)
    # read ATL03 data from downloaded HDF5 file
    mds,attrs,beams = is2tk.io.ATL03.read_granule(ids[0],
        ATTRIBUTES=False,
        VERBOSE=True
    )
    assert all(gtx in mds.keys() for gtx in beams)

# PURPOSE: Download an ATL06 file from NSIDC and check that read program runs
def test_ATL06_download_and_read(username, password):
    # query CMR for ATL06 data
    ids, urls = is2tk.utilities.cmr(product='ATL06',
        release='006', cycles=1, tracks=235, granules=2,
        verbose=False)
    # attempt download from on-prem NSIDC
    buffer, error = is2tk.utilities.from_nsidc(urls[0],
        username=username, password=password,
        local=ids[0], verbose=True)
    # attempt download from AWS
    if not buffer:
        url = posixpath.dirname(is2tk.utilities._s3_endpoints['nsidc'])
        bucket = is2tk.utilities._s3_buckets['nsidc']
        key = is2tk.utilities.s3_key(urls[0])
        buffer, error = is2tk.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=ids[0],
            verbose=True)
    # raise exception if download error
    if not buffer:
        raise Exception(error)
    # read ATL06 data from downloaded HDF5 file
    mds,attrs,beams = is2tk.io.ATL06.read_granule(ids[0],
        ATTRIBUTES=False,
        VERBOSE=True
    )
    assert all(gtx in mds.keys() for gtx in beams)

# PURPOSE: Download an ATL07 file from NSIDC and check that read program runs
def test_ATL07_download_and_read(username,password):
    # query CMR for ATL07 data
    ids, urls = is2tk.utilities.cmr(product='ATL07',
        release='006', cycles=1, tracks=235, regions=1,
        verbose=False)
    # attempt download from on-prem NSIDC
    buffer, error = is2tk.utilities.from_nsidc(urls[0],
        username=username, password=password,
        local=ids[0], verbose=True)
    # attempt download from AWS
    if not buffer:
        url = posixpath.dirname(is2tk.utilities._s3_endpoints['nsidc'])
        bucket = is2tk.utilities._s3_buckets['nsidc']
        key = is2tk.utilities.s3_key(urls[0])
        buffer, error = is2tk.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=ids[0],
            verbose=True)
    # raise exception if download error
    if not buffer:
        raise Exception(error)
    # read ATL07 data from downloaded HDF5 file
    mds,attrs,beams = is2tk.io.ATL07.read_granule(ids[0],
        ATTRIBUTES=False,
        VERBOSE=True
    )
    assert all(gtx in mds.keys() for gtx in beams)

# PURPOSE: Download an ATL11 file from NSIDC and check that read program runs
def test_ATL11_download_and_read(username,password):
    # query CMR for ATL11 data
    ids, urls = is2tk.utilities.cmr(product='ATL11',
        release='006', tracks=1, granules=3, verbose=False)
    # attempt download from on-prem NSIDC
    buffer, error = is2tk.utilities.from_nsidc(urls[0],
        username=username, password=password,
        local=ids[0], verbose=True)
    # attempt download from AWS
    if not buffer:
        url = posixpath.dirname(is2tk.utilities._s3_endpoints['nsidc'])
        bucket = is2tk.utilities._s3_buckets['nsidc']
        key = is2tk.utilities.s3_key(urls[0])
        buffer, error = is2tk.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=ids[0],
            verbose=True)
    # raise exception if download error
    if not buffer:
        raise Exception(error)
    # read ATL11 data from downloaded HDF5 file
    GROUPS = ['cycle_stats','ref_surf','crossing_track_data']
    mds,attrs,pairs = is2tk.io.ATL11.read_granule(ids[0],
        ATTRIBUTES=False,
        GROUPS=GROUPS,
        VERBOSE=True
    )
    assert all(ptx in mds.keys() for ptx in pairs)
    ptx = pairs[0]
    assert all(group in mds[ptx].keys() for group in GROUPS)

# PURPOSE: Download an ATL12 file from NSIDC and check that read program runs
def test_ATL12_download_and_read(username,password):
    # query CMR for ATL12 data
    ids, urls = is2tk.utilities.cmr(product='ATL12',
        release='006', cycles=1, tracks=237, granules=1,
        verbose=False)
    # attempt download from on-prem NSIDC
    buffer, error = is2tk.utilities.from_nsidc(urls[0],
        username=username, password=password,
        local=ids[0], verbose=True)
    # attempt download from AWS
    if not buffer:
        url = posixpath.dirname(is2tk.utilities._s3_endpoints['nsidc'])
        bucket = is2tk.utilities._s3_buckets['nsidc']
        key = is2tk.utilities.s3_key(urls[0])
        buffer, error = is2tk.utilities.from_nsidc([url,bucket,key],
            username=username, password=password, local=ids[0],
            verbose=True)
    # raise exception if download error
    if not buffer:
        raise Exception(error)
    # read ATL12 data from downloaded HDF5 file
    mds,attrs,beams = is2tk.io.ATL12.read_granule(ids[0],
        ATTRIBUTES=False,
        VERBOSE=True
    )
    assert all(gtx in mds.keys() for gtx in beams)

# PURPOSE: test CMR queries for specific cycles
def test_cmr_query_cycles():
    ids,urls = is2tk.utilities.cmr(product='ATL06',
        release='006',cycles=[2,3],tracks=752,granules=10,
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_006_02.h5',
        'ATL06_20190517125119_07520310_006_02.h5']
    assert all(id in valid for id in ids)

# PURPOSE: test CMR queries for specific tracks
def test_cmr_query_tracks():
    ids,urls = is2tk.utilities.cmr(product='ATL06',
        release='006',cycles=2,tracks=[752,753],granules=10,
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_006_02.h5',
        'ATL06_20190215184558_07530210_006_02.h5']
    assert all(id in valid for id in ids)

# PURPOSE: test CMR queries for specific granules
def test_cmr_query_granules():
    ids,urls = is2tk.utilities.cmr(product='ATL06',
        release='006',cycles=2,tracks=752,granules=[10,11,12],
        verbose=False)
    valid = ['ATL06_20190215171140_07520210_006_02.h5',
        'ATL06_20190215171921_07520211_006_02.h5',
        'ATL06_20190215172504_07520212_006_02.h5']
    assert all(id in valid for id in ids)
