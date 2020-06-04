ICESat-2 Data Products
======================

| Product | Name | Description |
| ------- | ---- | ----------- |
| **ATL00** | Telemetry Data | Raw ATLAS telemetry in packet format|
| **ATL01** | Reformatted Telemetry | Parsed, partially reformatted into HDF5, generated daily, segmented into several minute granules|
| **ATL02** | Science Unit Converted Telemetry | Photon time of flight, corrected for instrument effects. Includes all photons, pointing data, spacecraft position, housekeeping data, engineering data, and raw atmospheric profiles, segmented into several minute granules.|
| **ATL03** | Global Geolocated Photon Data | Precise latitude, longitude and elevation for every received photon, arranged by beam in the along-track direction. Photons classified by signal vs. background, as well as by surface type (land ice, sea ice, land, ocean), including all geophysical corrections. Segmented into several minute granules.|
| **ATL04** | Uncalibrated Backscatter Profiles | Along-track atmospheric backscatter data, 25 times per second. Includes calibration coefficients for polar regions. Segmented into several minute granules.|
| **ATL06** | Land Ice Elevation | Surface height for each beam with along- and across-track slopes calculated for each beam pair. Posted at 40m along-track; segmented into several minute granules.|
| **ATL07** | Arctic/Antarctic Sea Ice Elevation | Height of sea ice and open water leads at varying length scale based on returned photon rate for each beam presented along-track|
| **ATL08** | Land Water Vegetation Elevation | Height of ground including canopy surface posted at variable length scales relative to signal level, for each beam presented along-track. Where data permits include canopy height, canopy cover percentage, surface slope and roughness, and apparent reflectance.|
| **ATL09** | Calibrated Backscatter and Cloud Characteristics  | Along-track cloud and other significant atmosphere layer heights, blowing snow, integrated backscatter, and optical depth.|
| **ATL10** | Arctic/Antarctic Sea Ice Freeboard | Estimate of sea ice freeboard over specific spatial scales using all available sea surface height measurements. Contains statistics of sea surface and sea ice heights.|
| **ATL11** | Antarctica / Greenland Ice Sheet H(t) Series | Time series of height at points on the ice sheet, calculated based on repeat tracks and/or cross-overs.|
| **ATL12** | Ocean Elevation | Surface height at specific length scale.|
| **ATL13** | Inland Water Height | Along-track inland and near shore water surface height distribution within water mask.|
| **ATL14** | Antarctica / Greenland Ice Sheet H(t) Gridded | Height maps of each ice sheet for each year based on all available elevation data.|
| **ATL15** | Antarctica / Greenland Ice Sheet dh/dt Gridded | Height change maps for each ice sheet, for each mission year, and for the whole mission.|
| **ATL16** | ATLAS Atmosphere Weekly | Polar cloud fraction, blowing snow frequency, ground detection frequency.|
| **ATL17** | ATLAS Atmosphere Monthly | Polar cloud fraction, blowing snow frequency, ground detection frequency.|
| **ATL18** | Land/Canopy Gridded | Gridded ground surface height, canopy height, and canopy cover estimates.|
| **ATL19** | Mean Sea Surface (MSS) | Gridded ocean height product.|
| **ATL20** | Arctic / Antarctic Gridded Sea Ice Freeboard | Gridded sea ice freeboard.|
| **ATL21** | Arctic/Antarctic Gridded Sea Surface Height w/in Sea Ice | Gridded monthly sea surface height inside the sea ice cover.|

### ICESat-2 Granules
Each orbit of ICESat-2 data is broken up into 14 granule regions.  The granule boundaries limit the size of each ATL03 file and simplify the formation of higher level data products.  
![ICESat-2-global-granules](https://raw.githubusercontent.com/tsutterley/read-ICESat-2/master/data/ICESat-2_granules_global.png)  
![ICESat-2-polar-granules](https://raw.githubusercontent.com/tsutterley/read-ICESat-2/master/data/ICESat-2_granules_polar.png)  
