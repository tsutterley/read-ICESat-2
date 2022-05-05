========================
nsidc_icesat2_dragann.py
========================

- Acquires the `ATL03 geolocated photon height product <https://nsidc.org/data/ATL03>`_ and appends the `ATL08 DRAGANN classifications <https://nsidc.org/sites/nsidc.org/files/technical-references/ICESat2_ATL08_ATBD_r003.pdf>`_ from NSIDC.
- The first time we run the script, it will copy the necessary dataset in the selected local directory.
- If we already have all the data, and we run the script again: only files added or modified on the remote server will downloaded.

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_dragann.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/nsidc_icesat2_dragann.py
    :func: arguments
    :prog: nsidc_icesat2_dragann.py
    :nodescription:
    :nodefault:

    --granule -g : @replace
        Possible choices: 1--14

        ICESat-2 Granule Region

    --track -t : @replace
        Possible choices: 1--1387

        ICESat-2 Reference Ground Tracks (RGTs)
