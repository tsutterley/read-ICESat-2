=====================
read_ICESat2_ATL03.py
=====================

- Reads ICESat-2 ATL03 Global Geolocated Photons and ATL09 Atmospheric Characteristics data files

Calling Sequence
################

.. code-block:: python

    from icesat2_toolkit.read_ICESat2_ATL03 import read_HDF5_ATL03
    IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams = read_HDF5_ATL03(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL03.py

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL03.read_HDF5_ATL03

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL03.read_HDF5_ATL09

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL03.find_HDF5_ATL03_beams

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL03.read_HDF5_ATL03_main

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL03.read_HDF5_ATL03_beam
