=============
Parallel HDF5
=============

- Parallel HDF5 is a way to use HDF5 in a parallel computing environment
- Configuration of the HDF5 library that uses MPI (Message Passing Interface) to share open files across multiple parallel processes
- Some software provided here has been build using Parallel HDF5 to distribute the computational load across nodes

Binary Installation (Linux)
###########################

.. code-block:: bash

    sudo apt-get install openmpi-bin mpi-default-bin openmpi-doc libopenmpi-dev
    sudo apt-get install libhdf5-openmpi-dev
    sudo apt-get install python3-setuptools python3-pip
    env MPICC=/usr/bin/mpicc python3 -m pip install --user mpi4py
    export CC=mpicc
    export HDF5_MPI="ON"
    python3 -m pip install --user --no-binary=h5py h5py


Binary Installation (MacPorts)
##############################

.. code-block:: bash

    sudo port install hdf5 +openmpi_devel
    sudo port install py37-h5py +openmpi_devel

Source Installation
###################

The versions of each program may be different based on software upgrades and your computer version.
This is only meant to guide you along your way.

Dependencies
------------

- `zlib <https://www.zlib.net/>`_

.. code-block:: bash

    curl -O http://zlib.net/zlib-1.2.11.tar.gz
    export CFLAGS=-fPIC
    mkdir -p $HOME/packages/zlib/1.2.11
    ./configure --prefix=$HOME/packages/zlib/1.2.11
    make
    make install

modulefile for local installation of zlib (``~/privatemodules/zlib/1.2.11``):

.. code-block:: tcl

    #%Module 1.0
    #
    #  zlib module for use with 'environment-modules' package:
    #
    module-whatis		"Provides zlib 1.2.11 (local)"
    global              env
    prepend-path        PATH            $env(HOME)/packages/zlib/1.2.11/bin
    prepend-path        LD_LIBRARY_PATH $env(HOME)/packages/zlib/1.2.11/lib
    prepend-path        MANPATH         $env(HOME)/packages/zlib/1.2.11/share/man/
    append-path         ZLIB_DIR        $env(HOME)/packages/zlib/1.2.11/


- `szip <https://support.hdfgroup.org/doc_resource/SZIP/>`_

.. code-block:: bash

    curl -O https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
    export CFLAGS=-fPIC
    mkdir -p $HOME/packages/szip/2.1.1
    ./configure --prefix=$HOME/packages/szip/2.1.1
    make
    make install

modulefile for local installation of szip (``~/privatemodules/szip/2.1.1``):

.. code-block:: tcl

    #%Module 1.0
    #
    #  szip module for use with 'environment-modules' package:
    #
    module-whatis		"Provides szip 2.1.1 (local)"
    global              env
    prepend-path        PATH            $env(HOME)/packages/szip/2.1.1/bin
    prepend-path        LD_LIBRARY_PATH $env(HOME)/packages/szip/2.1.1/lib
    prepend-path        MANPATH         $env(HOME)/packages/szip/2.1.1/share/man/


- `OpenMPI <https://www.open-mpi.org/>`_

    * https://www.open-mpi.org/software/ompi/v4.0/
    * https://www.open-mpi.org/faq/?category=building

.. code-block:: bash

    curl -O https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz
    mkdir -p $HOME/packages/mpi/openmpi/4.0.3
    ./configure --prefix=$HOME/packages/mpi/openmpi/4.0.3
    make all install

modulefile for local installation of OpenMPI (``~/privatemodules/mpi/openmpi/4.0.3``):

.. code-block:: tcl

    #%Module 1.0
    #
    #  OpenMPI module for use with 'environment-modules' package:
    #
    module-whatis		"Provides openmpi 4.0.3 (local)"
    global              env
    prepend-path        OPAL_PREFIX      $env(HOME)/packages/mpi/openmpi/4.0.3
    prepend-path        PATH             $env(HOME)/packages/mpi/openmpi/4.0.3/bin
    prepend-path        LD_LIBRARY_PATH  $env(HOME)/packages/mpi/openmpi/4.0.3/lib
    prepend-path        MANPATH          $env(HOME)/packages/mpi/openmpi/4.0.3/share/man/


- `HDF5 <https://www.hdfgroup.org>`_

.. code-block:: bash

    curl -O https://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
    mkdir -p $HOME/packages/hdf5/1.10.5
    CC=~/packages/mpi/openmpi/4.0.3/bin/mpicc ./configure \
        --enable-parallel --enable-hl --enable-shared \
        --prefix=$HOME/packages/hdf5/1.10.5 \
        --with-zlib=$HOME/packages/zlib/1.2.11 \
        --with-szip=$HOME/packages/szip/2.1.1
    make
    make check
    make install

modulefile for local installation of HDF5 (``~/privatemodules/hdf5/1.10.5``):

.. code-block:: tcl

    #%Module 1.0
    #
    #  HDF5 module for use with 'environment-modules' package:
    #
    module-whatis		"Provides hdf5 1.10.5 (local)"
    global              env
    prereq	$env(HOME)/privatemodules/zlib/1.2.11	$env(HOME)/privatemodules/szip/2.1.1	$env(HOME)/privatemodules/mpi/openmpi/4.0.3
    prepend-path        PATH            $env(HOME)/packages/hdf5/1.10.5/bin
    prepend-path        LD_LIBRARY_PATH $env(HOME)/packages/hdf5/1.10.5/lib
    prepend-path        MANPATH         $env(HOME)/packages/hdf5/1.10.5/share/man/
    append-path         HDF5_DIR        $env(HOME)/packages/hdf5/1.10.5/


- `Python <https://www.python.org/>`_

modulefile for root installation of Python3 (``~/privatemodules/python/3.6.9``):

.. code-block:: tcl

    #%Module 1.0
    #
    #  Python module for use with 'environment-modules' package:
    #
    module-whatis		"Provides Python 3.6.9"
    global              env
    set-alias           python      /usr/bin/python3
    set-alias           pip         /usr/bin/pip3
    prepend-path    PYTHONPATH      $env(HOME)/local/lib/python3.6/site-packages
    setenv          PYTHONUSERBASE  $env(HOME)/local/


Python Packages
---------------

- `mpi4py: Python bindings of the Message Passing Interface (MPI) <https://mpi4py.readthedocs.io/en/stable/>`_
- `h5py: Read and write HDF5 files from Python <http://www.h5py.org>`_

	* http://docs.h5py.org/en/stable/mpi.html
	* http://docs.h5py.org/en/stable/build.html

.. code-block:: bash

    python3 -m pip install --upgrade pip
    env MPICC=~/packages/mpi/openmpi/4.0.3/bin/mpicc python3 -m pip install --user mpi4py
    export CC=$HOME/packages/mpi/openmpi/4.0.3/bin/mpicc
    export HDF5_MPI="ON"
    python3 -m pip install --user --no-binary=h5py h5py

- If getting error showing Open MPI removed legacy modules from the MPI-1 standard: `install maintained mpi4py from bitbucket <https://bitbucket.org/mpi4py/mpi4py/issues/115/cannot-build-against-openmpi-400>`_

.. code-block:: bash

    env MPICC=~/packages/mpi/openmpi/4.0.3/bin/mpicc \
    python3 -m pip install --user \
    https://bitbucket.org/mpi4py/mpi4py/get/maint.zip
