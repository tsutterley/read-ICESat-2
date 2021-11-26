==========
spatial.py
==========

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


.. method:: icesat2_toolkit.spatial.case_insensitive_filename(filename)

    Searches a directory for a filename without case dependence


.. method:: icesat2_toolkit.spatial.from_file(filename, format, **kwargs):

    Wrapper function for reading data from an input format

    Arguments:

        full path of input file

        format of input file


.. method:: icesat2_toolkit.spatial.from_netCDF4(filename, compression=None, verbose=False, xname='x', yname='y', varname='data')

    Read data from a netCDF4 file

    Arguments: full path of input netCDF4 file

    Keyword arguments:

        ``compression`` netCDF4 file is compressed or streamed from memory

        ``verbose`` print netCDF4 file information

        ``xname`` input x variable name in netCDF4 file

        ``yname`` input y variable units in netCDF4 file

        ``varname`` input data variable units in netCDF4 file


.. method:: icesat2_toolkit.spatial.from_HDF5(filename, compression=None, verbose=False, xname='x', yname='y', varname='data')

    Read data from a HDF5 file

    Arguments: full path of input HDF5 file

    Keyword arguments:

        ``compression`` HDF5 file is compressed or streamed from memory

        ``verbose`` print HDF5 file information

        ``xname`` input x variable name in HDF5 file

        ``yname`` input y variable units in HDF5 file

        ``varname`` input data variable units in HDF5 file


.. method:: icesat2_toolkit.spatial.from_geotiff(filename, compression=None, verbose=False)

    Read data from a geotiff file

    Arguments: full path of input geotiff file

    Keyword arguments:

        ``compression`` geotiff file is compressed using gzip

        ``verbose`` print geotiff filename


.. method:: icesat2_toolkit.spatial.convert_ellipsoid(phi1, h1, a1, f1, a2, f2, eps=1e-12, itmax=10)

    Convert latitudes and heights to a different ellipsoid using Newton-Raphson

    Arguments:

        ``phi1``: latitude of input ellipsoid in degrees

        ``h1``: height above input ellipsoid in meters

        ``a1``: semi-major axis of input ellipsoid

        ``f1``: flattening of input ellipsoid

        ``a2``: semi-major axis of output ellipsoid

        ``f2``: flattening of output ellipsoid

    Keyword arguments:

        ``eps``: tolerance to prevent division by small numbers and to determine convergence

        ``itmax``: maximum number of iterations to use in Newton-Raphson

    Returns:

        ``phi2``: latitude of output ellipsoid in degrees

        ``h2``: height above output ellipsoid in meters


.. method:: icesat2_toolkit.spatial.compute_delta_h(a1, f1, a2, f2, lat)

    Compute difference in elevation for two ellipsoids at a given latitude using a simplified empirical equation

    Arguments:

        ``a1``: semi-major axis of input ellipsoid

        ``f1``: flattening of input ellipsoid

        ``a2``: semi-major axis of output ellipsoid

        ``f2``: flattening of output ellipsoid

        ``lat``: array of latitudes in degrees

    Returns:

        ``delta_h``: difference in elevation for two ellipsoids


.. method:: icesat2_toolkit.spatial.wrap_longitudes(lon):

    Wraps longitudes to range from -180 to +180

    Arguments:

        ``lon``: longitude


.. method:: icesat2_toolkit.spatial.to_cartesian(lon,lat,a_axis=6378137.0,flat=1.0/298.257223563)

    Converts geodetic coordinates to Cartesian coordinates

    Arguments:

        ``lon``: longitude

        ``lat``: latitude

    Keyword arguments:

        ``h``: height

        ``a_axis``: semimajor axis of the ellipsoid

        ``flat``: ellipsoidal flattening

    Returns:

        ``x``, ``y``, ``z`` in Cartesian coordinates


.. method:: icesat2_toolkit.spatial.to_sphere(x,y,z)

    Convert from Cartesian coordinates to spherical coordinates

    Arguments:

        ``x``, ``y``, ``z`` in Cartesian coordinates

    Returns:

        ``lon``: longitude

        ``lat``: latitude

        ``rad``: radius


.. method:: icesat2_toolkit.spatial.to_geodetic(x,y,z,a_axis=6378137.0,flat=1.0/298.257223563)

    Convert from Cartesian coordinates to geodetic coordinates using `a closed form solution <https://arc.aiaa.org/doi/abs/10.2514/3.21016>`_

    Arguments:

        ``x``, ``y``, ``z`` in Cartesian coordinates

    Keyword arguments:

        ``a_axis``: semimajor axis of the ellipsoid

        ``flat``: ellipsoidal flattening

    Returns:

        ``lon``: longitude

        ``lat``: latitude

        ``h``: height


.. method:: icesat2_toolkit.spatial.scale_areas(lat, flat=1.0/298.257223563, ref=70.0)

    Calculates area scaling factors for a polar stereographic projection

    Arguments:

        ``lat``: latitude

    Keyword arguments:

        ``flat``: ellipsoidal flattening

        ``ref``: reference latitude (true scale latitude)

    Returns:

        ``scale``: area scaling factors at input latitudes
