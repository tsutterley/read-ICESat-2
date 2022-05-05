=========================
copy_scf_ICESat2_files.py
=========================

- Copies ICESat-2 HDF5 files from the SCF server

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/copy_scf_ICESat2_files.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/copy_scf_ICESat2_files.py
    :func: arguments
    :prog: copy_scf_ICESat2_files.py
    :nodescription:
    :nodefault:

    --granule -g : @replace
        Possible choices: 1--14

        ICESat-2 Granule Region

    --track -t : @replace
        Possible choices: 1--1387

        ICESat-2 Reference Ground Tracks (RGTs)
