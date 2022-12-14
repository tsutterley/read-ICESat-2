========
io.ATL11
========

- Reads ICESat-2 ATL11 Annual Land Ice Height data files

.. code-block:: python

    from icesat2_toolkit.io import ATL11
    IS2_atl11_mds,IS2_atl11_attrs,IS2_atl11_pairs = ATL11.read_granule(FILENAME,
        GROUPS=['cycle_stats','crossing_track_data'], VERBOSE=True)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL11.py

.. autofunction:: icesat2_toolkit.io.ATL11.read_granule

.. autofunction:: icesat2_toolkit.io.ATL11.find_pairs

.. autofunction:: icesat2_toolkit.io.ATL11.read_pair
