========
io.ATL07
========

- Reads ICESat-2 ATL07 Sea Ice Height data files

Calling Sequence
################

.. code-block:: python

    from icesat2_toolkit.io import ATL07
    IS2_atl07_mds,IS2_atl07_attrs,IS2_atl07_beams = ATL07.read_granule(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL07.py

.. autofunction:: icesat2_toolkit.io.ATL07.read_granule

.. autofunction:: icesat2_toolkit.io.ATL07.find_beams

.. autofunction:: icesat2_toolkit.io.ATL07.read_beam
