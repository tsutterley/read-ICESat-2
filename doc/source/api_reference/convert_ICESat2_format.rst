=========================
convert_ICESat2_format.py
=========================

- Converts ICESat-2 HDF5 datafiles to zarr or rechunked HDF5 datafiles for a specified data product, release, granule and track.
- zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores
- rechunked HDF5 files can be more optimized for cloud-based object stores

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/convert_ICESat2_format.py

Calling Sequence
################

.. argparse::
    :filename: convert_ICESat2_format.py
    :func: arguments
    :prog: convert_ICESat2_format.py
    :nodescription:
    :nodefault:

    --granule -g : @replace
        Possible choices: 1--14

        ICESat-2 Granule Region

    --track -t : @replace
        Possible choices: 1--1387

        ICESat-2 Reference Ground Tracks (RGTs)
