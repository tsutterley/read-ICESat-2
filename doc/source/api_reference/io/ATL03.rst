========
io.ATL03
========

- Reads ICESat-2 ATL03 Global Geolocated Photons
- Interpolates ATL09 Atmospheric Characteristics to ATL03 segments

Calling Sequence
################

.. code-block:: python

    from icesat2_toolkit.io import ATL03
    IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams = ATL03.read_granule(FILENAME)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/io/ATL03.py

.. autofunction:: icesat2_toolkit.io.ATL03.read_granule

.. autofunction:: icesat2_toolkit.io.ATL03.interpolate_ATL09

.. autofunction:: icesat2_toolkit.io.ATL03.find_beams

.. autofunction:: icesat2_toolkit.io.ATL03.read_main

.. autofunction:: icesat2_toolkit.io.ATL03.read_beam
