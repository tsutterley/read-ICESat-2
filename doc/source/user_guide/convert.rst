==========
convert.py
==========

Utilities for converting ICESat-2 HDF5 files into different formats

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/convert.py


General Methods
===============

.. class:: convert(object)

    .. attribute:: object.filename

        input HDF5 filename or BytesIO object


    .. attribute:: object.reformat

        output format (``'zarr'``, ``'HDF'``, ``'JPL'``, ``'csv'``, ``'txt'``, ``'dataframe'``)


    .. method:: object.file_converter(**kwds)

        Convert a HDF5 file to another format


    .. method:: object.HDF5_to_zarr(**kwds)

        Convert a HDF5 file to zarr copying all file data


    .. method:: object.HDF5_to_HDF5(**kwds)

        Rechunk a HDF5 file copying all file data


    .. method:: object.copy_from_HDF5(source, dest, name=None, **kwds)

        Copy a named variable from the `source` HDF5 into the `dest` file


    .. method:: object.HDF5_to_ascii(**kwds)

        Convert a HDF5 file to beam-level ascii files copying reduced sets of data


    .. method:: object.HDF5_to_dataframe(**kwds)

        Convert a HDF5 file to a pandas dataframe copying reduced sets of data
