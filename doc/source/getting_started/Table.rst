+-----------+------------------------------+-------------------------------------------+
|  Product  |             Name             |                Description                |
+===========+==============================+===========================================+
| **ATL00** | Telemetry Data               | Raw ATLAS telemetry in packet format      |
+-----------+------------------------------+-------------------------------------------+
| **ATL01** | Reformatted Telemetry        | Parsed, partially reformatted into HDF5,  |
|           |                              | generated daily, segmented into several   |
|           |                              | minute granules                           |
+-----------+------------------------------+-------------------------------------------+
| **ATL02** | Science Unit Converted       | Photon time of flight, corrected for      |
|           | Telemetry                    | instrument effects. Includes all photons, |
|           |                              | pointing data, spacecraft position,       |
|           |                              | housekeeping data, engineering data, and  |
|           |                              | raw atmospheric profiles, segmented into  |
|           |                              | several minute granules.                  |
+-----------+------------------------------+-------------------------------------------+
| **ATL03** | Global Geolocated Photon     | Precise latitude, longitude and elevation |
|           | Data                         | for every received photon, arranged by    |
|           |                              | beam in the along-track direction.        |
|           |                              | Photons classified by signal vs.          |
|           |                              | background, as well as by surface type    |
|           |                              | (land ice, sea ice, land, ocean),         |
|           |                              | including all geophysical corrections.    |
|           |                              | Segmented into several minute granules.   |
+-----------+------------------------------+-------------------------------------------+
| **ATL04** | Uncalibrated Backscatter     | Along-track atmospheric backscatter data, |
|           | Profiles                     | 25 times per second. Includes calibration |
|           |                              | coefficients for polar regions. Segmented |
|           |                              | into several minute granules.             |
+-----------+------------------------------+-------------------------------------------+
| **ATL06** | Land Ice Elevation           | Surface height for each beam with along-  |
|           |                              | and across-track slopes calculated for    |
|           |                              | each beam pair. Posted at 40 meters       |
|           |                              | along-track. Segmented into several       |
|           |                              | minute granules.                          |
+-----------+------------------------------+-------------------------------------------+
| **ATL07** | Arctic/Antarctic Sea Ice     | Height of sea ice and open water leads at |
|           | Elevation                    | varying length scale based on returned    |
|           |                              | photon rate for each beam presented       |
|           |                              | along-track. Segmented into several       |
|           |                              | minute granules.                          |
+-----------+------------------------------+-------------------------------------------+
| **ATL08** | Land Water Vegetation        | Height of ground including canopy surface |
|           | Elevation                    | posted at variable length scales relative |
|           |                              | to signal level, for each beam presented  |
|           |                              | along-track. Where data permits include   |
|           |                              | canopy height, canopy cover percentage,   |
|           |                              | surface slope and roughness, and apparent |
|           |                              | reflectance.                              |
+-----------+------------------------------+-------------------------------------------+
| **ATL09** | Calibrated Backscatter and   | Along-track cloud and other significant   |
|           | Cloud Characteristics        | atmosphere layer heights, blowing snow,   |
|           |                              | integrated backscatter, and optical       |
|           |                              | depth.                                    |
+-----------+------------------------------+-------------------------------------------+
| **ATL10** | Arctic/Antarctic Sea Ice     | Estimate of sea ice freeboard over        |
|           | Freeboard                    | specific spatial scales using all         |
|           |                              | available sea surface height              |
|           |                              | measurements. Contains statistics of sea  |
|           |                              | surface and sea ice heights.              |
+-----------+------------------------------+-------------------------------------------+
| **ATL11** | Antarctic/Greenland Ice      | Time series of height at points on the    |
|           | Sheet H(t) Series            | ice sheet, calculated based on repeat     |
|           |                              | tracks and/or cross-overs.                |
+-----------+------------------------------+-------------------------------------------+
| **ATL12** | Ocean Elevation              | Surface height at specific length scale.  |
+-----------+------------------------------+-------------------------------------------+
| **ATL13** | Inland Water Height          | Along-track inland and near shore water   |
|           |                              | surface height distribution within water  |
|           |                              | mask.                                     |
+-----------+------------------------------+-------------------------------------------+
| **ATL14** | Antarctic/Greenland Ice      | Height maps of each ice sheet for each    |
|           | Sheet H(t) Gridded           | year based on all available elevation     |
|           |                              | data.                                     |
+-----------+------------------------------+-------------------------------------------+
| **ATL15** | Antarctic/Greenland Ice      | Height change maps for each ice sheet,    |
|           | Sheet dh/dt Gridded          | for each mission year, and for the whole  |
|           |                              | mission.                                  |
+-----------+------------------------------+-------------------------------------------+
| **ATL16** | ATLAS Atmosphere Weekly      | Polar cloud fraction, blowing snow        |
|           |                              | frequency, ground detection frequency.    |
+-----------+------------------------------+-------------------------------------------+
| **ATL17** | ATLAS Atmosphere Monthly     | Polar cloud fraction, blowing snow        |
|           |                              | frequency, ground detection frequency.    |
+-----------+------------------------------+-------------------------------------------+
| **ATL18** | Land/Canopy Gridded          | Gridded ground surface height, canopy     |
|           |                              | height, and canopy cover estimates.       |
+-----------+------------------------------+-------------------------------------------+
| **ATL19** | Mean Sea Surface (MSS)       | Gridded ocean height product.             |
+-----------+------------------------------+-------------------------------------------+
| **ATL20** | Arctic/Antarctic Gridded Sea | Gridded sea ice freeboard.                |
|           | Ice Freeboard                |                                           |
+-----------+------------------------------+-------------------------------------------+
| **ATL21** | Arctic/Antarctic Gridded Sea | Gridded monthly sea surface height inside |
|           | Surface Height w/in Sea Ice  | the sea ice cover.                        |
+-----------+------------------------------+-------------------------------------------+
