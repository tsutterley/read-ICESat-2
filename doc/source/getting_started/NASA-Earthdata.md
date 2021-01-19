NASA Earthdata
==============

#### NASA Data Distribution Centers
The NASA Earth Science Data Information Systems Project funds and operates [12 Distributed Active Archive Centers (DAACs)](https://earthdata.nasa.gov/about/daacs) throughout the United States.  These centers have recently transitioned from ftp to https servers.
The https updates are designed to increase performance and improve security during data retrieval. NASA Earthdata uses [OAuth2](https://wiki.earthdata.nasa.gov/pages/viewpage.action?pageId=71700485), an approach to authentication that protects your personal information.
- https://urs.earthdata.nasa.gov/documentation
- https://wiki.earthdata.nasa.gov/display/EL/Knowledge+Base
- https://nsidc.org/support/how/how-do-i-programmatically-request-data-services
- https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-https-earthdata-login-enabled

#### NSIDC
The [National Snow and Ice Data Center (NSIDC)](https://nsidc.org/daac/) DAAC provides data and information for snow and ice processes, particularly interactions among snow, ice, atmosphere, and ocean, in support of research in global change detection and model validation. If any problems contact NSIDC support at [nsidc@nsidc.org](mailto:nsidc@nsidc.org) or the NASA EOSDIS support team [support@earthdata.nasa.gov](mailto:support@earthdata.nasa.gov).

#### Steps to Sync from NSIDC
1. [Register with NASA Earthdata Login system](https://urs.earthdata.nasa.gov/users/new)
2. [After registering, login to the system](https://urs.earthdata.nasa.gov/home)
3. Add `NSIDC_DATAPOOL_OPS` and `nsidc-daacdata` [applications to Earthdata](https://wiki.earthdata.nasa.gov/display/EL/How+To+Pre-authorize+an+application)
4. Copy your NASA Earthdata credentials or [create a .netrc file](https://nsidc.org/support/how/v0-programmatic-data-access-guide) to store your credentials permanently
```bash
echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> ~/.netrc
chmod 0600 ~/.netrc
```
5. Sync data from NSIDC as [HDF5](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_sync.py) or [zarr](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_convert.py), or use the [NSIDC subsetting API](https://github.com/tsutterley/nsidc-subsetter) to gather data

#### Other Data Access Examples
- [Curl and Wget](https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget)
- [Python](https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python)
- [Jupyter](https://github.com/nsidc/NSIDC-Data-Access-Notebook)
- [icepyx](https://github.com/icesat2py/icepyx)
