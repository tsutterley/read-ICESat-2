========
io.ATL06
========

- Reads ICESat-2 ATL06 Land Ice Elevation data files

Calling Sequence
################

.. code-block:: python

    from icesat2_toolkit.io import ATL06
    IS2_atl06_mds,IS2_atl06_attrs,IS2_atl06_beams = ATL06.read_granule(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL06.py

.. autofunction:: icesat2_toolkit.io.ATL06.read_granule

.. autofunction:: icesat2_toolkit.io.ATL06.find_beams

.. autofunction:: icesat2_toolkit.io.ATL06.read_beam
