===================================
interp_IB_response_ICESat2_ATL07.py
===================================

- Calculates inverse-barometer (IB) responses for correcting ICESat-2 sea ice height data [Wunsch1997]_ [HofmannWellenhof2006]_
- Can use mean sea level pressure outputs from `ERA-Interim <http://apps.ecmwf.int/datasets/data/interim-full-moda>`_, `ERA5 <http://apps.ecmwf.int/data-catalogues/era5/?class=ea>`_ and `MERRA-2 <https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/>`_ reanalyses

Calling Sequence
################

.. code-block:: bash

    python interp_IB_response_ICESat2_ATL07.py --directory <path_to_directory> --reanalysis <model> input_file

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/interp_IB_response_ICESat2_ATL07.py

Inputs
######

1. ``input_file``: input ICESat-2 ATL07 file

Command Line Options
####################

- ``-D X``, ``--directory X``: Working data directory
- ``-R X``, ``--reanalysis X``: Reanalysis model to run
    * ``'ERA-Interim'``
    * ``'ERA5'``
    * ``'MERRA-2'``
- ``-m X``, ``--mean X``: Start and end year range for mean
- ``-d X``, ``--density X``: Density of seawater in kg/m\ :sup:`3`
- ``-V``, ``--verbose``: Output information about each created file
- ``-M X``, ``--mode X``: Permission mode of output file

References
##########

.. [Wunsch1997] Wunsch and Stammer. "Atmospheric loading and the oceanic "inverted barometer" effect", *Reviews of Geophysics*, 35(1), 79-107, (1997). `doi:10.1029/96RG03037 <https://doi.org/10.1029/96RG03037>`_
.. [HofmannWellenhof2006] B. Hofmann-Wellenhof and H. Moritz, *Physical Geodesy*, 2nd Edition, 403 pp., (2006). `doi: 10.1007/978-3-211-33545-1 <https://doi.org/10.1007/978-3-211-33545-1>`_
