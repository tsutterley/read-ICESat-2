========
io.ATL12
========

- Reads ICESat-2 ATL12 Ocean Surface Height data files

.. code-block:: python

    from icesat2_toolkit.io import ATL12
    IS2_atl12_mds,IS2_atl12_attrs,IS2_atl12_beams = ATL12.read_granule(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL12.py

.. autofunction:: icesat2_toolkit.io.ATL12.read_granule

.. autofunction:: icesat2_toolkit.io.ATL12.find_beams

.. autofunction:: icesat2_toolkit.io.ATL12.read_beam
