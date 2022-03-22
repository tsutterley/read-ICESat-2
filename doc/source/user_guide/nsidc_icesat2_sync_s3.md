nsidc_icesat2_sync_s3.py
========================

- Syncs all available ICESat-2 data for a specified data product, release, granule and track and transfers to an AWS S3 bucket using a local machine as pass through

#### Calling Sequence
```bash
python nsidc_icesat2_sync_s3.py --user <username> --aws-access-key-id <access-key> \
	--aws-secret-access-key <secret-key> --s3-bucket-name <bucket> --release 004 \
	--granule 10 11 12 --flatten ATL06
```
[Source code](https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/nsidc_icesat2_sync_s3.py)

#### Command Line Options
- `-U X`, `--user X`: username for NASA Earthdata Login
- `-W X`, `--password X`: password for NASA Earthdata Login
- `-N X`, `--netrc X`: path to .netrc file for alternative authentication
- `--aws-access-key-id X`: AWS Access Key ID
- `--aws-secret-access-key X`: AWS Secret Key
- `--aws-region-name X`: AWS Region Name
- `--s3-bucket-name X`: AWS S3 Bucket Name
- `--s3-bucket-path X`: AWS S3 Bucket Path
- `-Y X`, `--year X`: years to sync
- `-S X`, `--subdirectory X`: specific subdirectories to sync
- `-r X`, `--release X`: ICESat-2 data release to sync
- `-v X`, `--version X:` ICESat-2 data version to sync
- `-t X`, `--track X`: ICESat-2 reference ground tracks to sync
- `-g X`, `--granule X`: ICESat-2 granule regions to sync
- `-c X`, `--cycle X`: ICESat-2 orbital cycles to sync
- `-n X`, `--region X`: ICESat-2 Named Region to sync for ATL14/ATL15
- `-a`, `--auxiliary`: Sync ICESat-2 auxiliary files for each HDF5 file
- `-I X, --index X`: Input index of ICESat-2 files to sync
- `-F`, `--flatten`: Do not create subdirectories
- `-P X`, `--np X`: Number of processes to use in file downloads
- `-T X`, `--timeout X`: Timeout in seconds for blocking operations
- `-R X`, `--retry X`: Connection retry attempts
- `-C`, `--clobber`: Overwrite existing data in transfer
