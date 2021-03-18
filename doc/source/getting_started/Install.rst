======================
Setup and Installation
======================

Dependencies
############
This software is dependent on open source programs that can be installed using OS-specific package management systems,
`conda <https://anaconda.org/conda-forge/repo>`_ or from source:

- `MPI <https://www.open-mpi.org/>`_
- `GDAL <https://gdal.org/index.html>`_
- `GEOS <https://trac.osgeo.org/geos>`_
- `PROJ <https://proj.org/>`_
- `HDF5 <https://www.hdfgroup.org>`_
- `libxml2 <http://xmlsoft.org/>`_
- `libxslt <http://xmlsoft.org/XSLT/>`_

Installation
############

Presently ``read-ICESat-2`` is only available for use as a
`GitHub repository <https://github.com/tsutterley/read-ICESat-2>`_.
The contents of the repository can be download as a
`zipped file <https://github.com/tsutterley/read-ICESat-2/archive/main.zip>`_  or cloned.
To use this repository, please fork into your own account and then clone onto your system.

.. code-block:: bash

    git clone https://github.com/tsutterley/read-ICESat-2.git

Can then install using ``setuptools``

.. code-block:: bash

    python setup.py install

or ``pip``

.. code-block:: bash

    python3 -m pip install --user .

Alternatively can install the utilities directly from GitHub with ``pip``:

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/read-ICESat-2.git

Executable versions of this repository can also be tested using
`Binder <https://mybinder.org/v2/gh/tsutterley/read-ICESat-2/main>`_ or
`Pangeo <https://binder.pangeo.io/v2/gh/tsutterley/read-ICESat-2/main>`_.
