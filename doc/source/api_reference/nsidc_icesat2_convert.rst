========================
nsidc_icesat2_convert.py
========================

- Syncs all available ICESat-2 data converting from HDF5 to zarr or rechunked HDF5 files for a specified data product, release, granule and track.
- zarr files make large datasets easily accessible to distributed computing on both local filesystems and cloud-based object stores
- rechunked HDF5 files can be more optimized for cloud-based object stores
- The first time we run the script, it will copy the necessary dataset in the selected local directory.
- If we already have all the data, and we run the script again: only files added or modified on the remote server will downloaded.

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_convert.py

Calling Sequence
################

.. argparse::
    :filename: nsidc_icesat2_convert.py
    :func: arguments
    :prog: nsidc_icesat2_convert.py
    :nodescription:
    :nodefault:

    --granule -g : @replace
        Possible choices: 1--14

        ICESat-2 Granule Region

    --track -t : @replace
        Possible choices: 1--1387

        ICESat-2 Reference Ground Tracks (RGTs)
