===============
Getting Started
===============

- `Register at NASA Earthdata and add NSIDC applications <./NASA-Earthdata.html>`_
- Retrieve ``ATL03_20191115042423_07520512_004_01.h5`` from NSIDC

.. code-block:: bash

    python nsidc_icesat2_sync.py --user <username> --directory <path_to_directory> \
        --subdirectory 2019.11.15 --release 004 --version 01 --track 752 \
        --granule 12 --flatten --mode 0o775 ATL03

- Run Jupyter notebook ``Read ICESat-2 ATL03.ipynb`` to learn about ICESat-2 and ATL03
