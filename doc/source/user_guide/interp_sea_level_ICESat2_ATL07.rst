=================================
interp_sea_level_ICESat2_ATL07.py
=================================

- Interpolates `AVISO sea level height estimates <https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/>`_ for correcting ICESat-2 sea ice height data

    * sea level anomalies (``sla``)
    * absolute dynamic topography (``adt``)
    * mean dynamic topography (``mdt``)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/interp_sea_level_ICESat2_ATL07.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/interp_sea_level_ICESat2_ATL07.py
    :func: arguments
    :prog: interp_sea_level_ICESat2_ATL07.py
    :nodescription:
    :nodefault:
