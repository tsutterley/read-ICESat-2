=====================
read_ICESat2_ATL11.py
=====================

- Reads ICESat-2 ATL11 Annual Land Ice Height data files

.. code-block:: python

    from icesat2_toolkit.read_ICESat2_ATL11 import read_HDF5_ATL11
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = read_HDF5_ATL11(FILENAME,
        GROUPS=['cycle_stats','crossing_track_data'], VERBOSE=True)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/read_ICESat2_ATL11.py

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL11.read_HDF5_ATL11

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL11.find_HDF5_ATL11_pairs

.. autofunction:: icesat2_toolkit.read_ICESat2_ATL11.read_HDF5_ATL11_pair
