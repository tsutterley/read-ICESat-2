=======================================
MPI_reduce_ICESat2_ATL11_ice_shelves.py
=======================================

- Create masks for reducing ICESat-2 ATL11 data into grounded ice regions
- `IceBridge BedMachine Greenland, Version 4 <https://doi.org/10.5067/VLJ5YXKCNGXO>`_
- `MEaSUREs Antarctic Boundaries for IPY 2007-2009 <https://doi.org/10.5067/AXE4121732AD>`_

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL11_grounded.py

Calling Sequence
################

.. argparse::
    :filename: ../../scripts/MPI_reduce_ICESat2_ATL11_grounded.py
    :func: arguments
    :prog: MPI_reduce_ICESat2_ATL11_grounded.py
    :nodescription:
    :nodefault:

Creating BedMachine grounded ice mask
#####################################

.. code-block:: bash

    # convert mask variable to separate geotiff
    gdal_translate -co "COMPRESS=LZW" -a_srs EPSG:3413 \
        NETCDF:BedMachineGreenland-2021-04-20.nc:mask \
        -of 'Gtiff' BedMachineGreenlandMaskv4.tif

    # create grounded ice mask (1=Greenland grounded, 0=other)
    gdalwarp -co "COMPRESS=LZW" -srcnodata 1 -dstnodata 0 -ot Byte BedMachineGreenlandMaskv4.tif tmp1.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 4 -dstnodata 0 -ot Byte tmp1.tif tmp2.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 3 -dstnodata 0 -ot Byte tmp2.tif tmp3.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 2 -dstnodata 1 -ot Byte tmp3.tif tmp4.tif
    gdalwarp -co "COMPRESS=LZW" -srcnodata 0 -dstnodata 0 -ot Byte tmp4.tif BedMachineGreenlandGroundedv4.tif
    rm tmp*.tif

    # create vectorized grounded ice mask
    gdal_polygonize.py BedMachineGreenlandGroundedv4.tif tmp.shp
    # add area attribute and convert units from m^2 to km^2
    ogr2ogr BedMachineGreenlandGroundedv4.shp tmp.shp -sql "SELECT *, OGR_GEOM_AREA/1000000 AS area FROM tmp"
    rm tmp.dbf tmp.prj tmp.shp tmp.shx
