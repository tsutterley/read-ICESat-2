=======================================
MPI_reduce_ICESat2_ATL06_ice_shelves.py
=======================================

- Create masks for reducing ICESat-2 ATL06 data into floating ice shelves
- `IceBridge BedMachine Greenland, Version 4 <https://doi.org/10.5067/VLJ5YXKCNGXO>`_
- `MEaSUREs Antarctic Boundaries for IPY 2007-2009 <https://doi.org/10.5067/AXE4121732AD>`_

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL06_ice_shelves.py

Calling Sequence
################

.. argparse::
    :filename: MPI_reduce_ICESat2_ATL06_ice_shelves.py
    :func: arguments
    :prog: MPI_reduce_ICESat2_ATL06_ice_shelves.py
    :nodescription:
    :nodefault:

Creating BedMachine floating ice mask
#####################################

.. code-block:: bash

    # convert mask variable to separate geotiff
    gdal_translate -co "COMPRESS=LZW" -a_srs EPSG:3413 \
        NETCDF:BedMachineGreenland-2021-04-20.nc:mask \
        -of 'Gtiff' BedMachineGreenlandMaskv4.tif

    # create floating ice mask (1=Greenland floating, 0=other)
    gdalwarp -co "COMPRESS=LZW" -srcnodata 1 -dstnodata 0 -ot Byte BedMachineGreenlandMaskv4.tif tmp1.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 2 -dstnodata 0 -ot Byte tmp1.tif tmp2.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 4 -dstnodata 0 -ot Byte tmp2.tif tmp3.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 3 -dstnodata 1 -ot Byte tmp3.tif tmp4.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 0 -dstnodata 0 -ot Byte tmp4.tif BedMachineGreenlandFloatingv4.tif
    rm tmp*.tif

    # create vectorized floating ice mask
    gdal_polygonize.py BedMachineGreenlandFloatingv4.tif tmp.shp
    # add area attribute and convert units from m^2 to km^2
    ogr2ogr BedMachineGreenlandFloatingv4.shp tmp.shp -sql "SELECT *, OGR_GEOM_AREA/1000000 AS area FROM tmp"
    rm tmp.dbf tmp.prj tmp.shp tmp.shx
