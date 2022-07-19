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

The version of GDAL used within ``read-ICESat-2`` will match the version of the installed C program.
The path to the C program that will be used with ``read-ICESat-2`` is given by:

.. code-block:: bash

    gdal-config --datadir

The ``read-ICESat-2`` installation uses the ``gdal-config`` routines to set the GDAL package version.

Installation
############

``read-ICESat-2`` is available for download from the `GitHub repository <https://github.com/tsutterley/read-ICESat-2>`_,
and the `Python Package Index (pypi) <https://pypi.org/project/icesat2-toolkit/>`_.
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

| This repository can be also tested using `BinderHub <https://github.com/jupyterhub/binderhub>`_ platforms:
| |Binder| |Pangeo|

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/tsutterley/read-ICESat-2/main

.. |Pangeo| image:: https://img.shields.io/static/v1.svg?logo=Jupyter&label=PangeoBinderAWS&message=us-west-2&color=orange
   :target: https://aws-uswest2-binder.pangeo.io/v2/gh/tsutterley/read-ICESat-2/main?urlpath=lab
