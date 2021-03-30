interp_sea_level_ICESat2_ATL07.py
=================================

- Interpolates [AVISO sea level height estimates](https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/) for correcting ICESat-2 sea ice height data
    * sea level anomalies (sla)
    * absolute dynamic topography (adt)
    * mean dynamic topography (mdt)

#### Calling Sequence
```bash
python interp_sea_level_ICESat2_ATL07.py --directory <path_to_directory> input_file
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/interp_sea_level_ICESat2_ATL07.py)

#### Inputs
1. `input_file`: input ICESat-2 ATL07 file

#### Command Line Options
- `-D X`, `--directory X`: Working data directory
- `-V`, `--verbose`: Output information about each created file
- `-M X`, `--mode X`: Permission mode of output file
