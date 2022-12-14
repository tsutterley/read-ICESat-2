=======
spatial
=======

Utilities for reading and operating on spatial data

 - Can read netCDF4, HDF5 or geotiff files

Calling Sequence
================

Reading a netCDF4 file

.. code-block:: python

    import icesat2_toolkit.spatial
    dinput = icesat2_toolkit.spatial.from_netCDF4(path_to_netCDF4_file)

Reading a HDF5 file

.. code-block:: python

    import icesat2_toolkit.spatial
    dinput = icesat2_toolkit.spatial.from_HDF5(path_to_HDF5_file)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/spatial.py

General Methods
===============

.. autofunction:: icesat2_toolkit.spatial.case_insensitive_filename

.. autofunction:: icesat2_toolkit.spatial.from_file

.. autofunction:: icesat2_toolkit.spatial.from_netCDF4

.. autofunction:: icesat2_toolkit.spatial.from_HDF5

.. autofunction:: icesat2_toolkit.spatial.from_geotiff

.. autofunction:: icesat2_toolkit.spatial.convert_ellipsoid

.. autofunction:: icesat2_toolkit.spatial.compute_delta_h

.. autofunction:: icesat2_toolkit.spatial.wrap_longitudes

.. autofunction:: icesat2_toolkit.spatial.to_cartesian

.. autofunction:: icesat2_toolkit.spatial.to_sphere

.. autofunction:: icesat2_toolkit.spatial.to_geodetic

.. autofunction:: icesat2_toolkit.spatial.scale_areas

.. autofunction:: icesat2_toolkit.spatial.inside_polygon
