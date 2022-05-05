=================================
interp_sea_level_ICESat2_ATL11.py
=================================

- Interpolates `AVISO sea level height estimates <https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/>`_ for correcting ICESat-2 annual land ice height data

    * sea level anomalies (``sla``)
    * absolute dynamic topography (``adt``)
    * mean dynamic topography (``mdt``)
- Interpolates sea level estimates for both along-track and across-track locations

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/interp_sea_level_ICESat2_ATL11.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/interp_sea_level_ICESat2_ATL11.py
    :func: arguments
    :prog: interp_sea_level_ICESat2_ATL11.py
    :nodescription:
    :nodefault:
