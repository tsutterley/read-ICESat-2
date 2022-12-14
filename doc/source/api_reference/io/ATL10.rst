========
io.ATL10
========

- Reads ICESat-2 ATL10 Sea Ice Freeboard data files

Calling Sequence
################

.. code-block:: python

    from icesat2_toolkit.io import ATL10
    IS2_atl10_mds,IS2_atl10_attrs,IS2_atl10_beams = ATL10.read_granule(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL10.py

.. autofunction:: icesat2_toolkit.io.ATL10.read_granule

.. autofunction:: icesat2_toolkit.io.ATL10.find_beams

.. autofunction:: icesat2_toolkit.io.ATL10.read_beam
