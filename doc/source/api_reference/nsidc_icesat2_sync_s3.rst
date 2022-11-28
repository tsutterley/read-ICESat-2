========================
nsidc_icesat2_sync_s3.py
========================

- Syncs all available ICESat-2 data for a specified data product, release, granule and track and transfers to an AWS S3 bucket using a local machine as pass through

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_sync_s3.py

Calling Sequence
################

.. argparse::
    :filename: nsidc_icesat2_sync_s3.py
    :func: arguments
    :prog: nsidc_icesat2_sync_s3.py
    :nodescription:
    :nodefault:

    --granule -g : @replace
        Possible choices: 1--14

        ICESat-2 Granule Region

    --track -t : @replace
        Possible choices: 1--1387

        ICESat-2 Reference Ground Tracks (RGTs)
