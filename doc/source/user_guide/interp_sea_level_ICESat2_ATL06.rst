=================================
interp_sea_level_ICESat2_ATL06.py
=================================

- Interpolates `AVISO sea level height estimates <https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/>`_ for correcting ICESat-2 land ice elevation data

    * sea level anomalies (``sla``)
    * absolute dynamic topography (``adt``)
    * mean dynamic topography (``mdt``)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/interp_sea_level_ICESat2_ATL06.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/interp_sea_level_ICESat2_ATL06.py
    :func: arguments
    :prog: interp_sea_level_ICESat2_ATL06.py
    :nodescription:
    :nodefault:
