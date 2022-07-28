========================
MPI_DEM_ICESat2_ATL11.py
========================

- Determines which digital elevation model tiles to read for a given ATL11 file
- Reads 3\ |times|\ 3 array of tiles for points within bounding box of central mosaic tile
- Interpolates digital elevation model to locations of ICESat-2 ATL11 segments

- ArcticDEM 2m digital elevation model tiles

    * `http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/ <http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/>`_
    * `http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ <http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/>`_

- REMA 8m digital elevation model tiles

    * `http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/ <http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/>`_
    * `http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/ <http://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/>`_

- GIMP 30m digital elevation model tiles

    * `https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/ <https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0645.001/>`_

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_DEM_ICESat2_ATL11.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/MPI_DEM_ICESat2_ATL11.py
    :func: arguments
    :prog: MPI_DEM_ICESat2_ATL11.py
    :nodescription:
    :nodefault:

.. |times|      unicode:: U+00D7 .. MULTIPLICATION SIGN
