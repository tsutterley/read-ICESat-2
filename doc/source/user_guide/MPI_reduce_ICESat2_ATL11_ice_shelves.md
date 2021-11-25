MPI_reduce_ICESat2_ATL11_ice_shelves.py
=======================================

- Create masks for reducing ICESat-2 data into floating ice shelves
- [IceBridge BedMachine Greenland, Version 4](https://doi.org/10.5067/VLJ5YXKCNGXO)
- [MEaSUREs Antarctic Boundaries for IPY 2007-2009](https://doi.org/10.5067/AXE4121732AD)

#### Calling Sequence
```bash
mpiexec -np <processes> python3 MPI_reduce_ICESat2_ATL11_ice_shelves.py <path_to_ATL11_file>
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL11_ice_shelves.py)

#### Inputs
1. `ATL11_file`: full path to ATL11 file

#### Command Line Options
- `-D X`, `--directory X`: Working data directory for ice shelf shapefiles
- `-B X`, `--buffer X`: Distance in kilometers to buffer ice shelves mask
- `-V`, `--verbose`: output module information for process
- `-M X`, `--mode X`: permissions mode of output HDF5 datasets

#### Creating BedMachine floating ice mask
```bash
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
```