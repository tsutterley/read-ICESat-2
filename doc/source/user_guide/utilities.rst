============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from NSIDC via https when Earthdata credentials are supplied
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/utilities.py


General Methods
===============

.. autofunction:: icesat2_toolkit.utilities.get_data_path

.. autofunction:: icesat2_toolkit.utilities.get_hash

.. autofunction:: icesat2_toolkit.utilities.url_split

.. autofunction:: icesat2_toolkit.utilities.get_unix_time

.. autofunction:: icesat2_toolkit.utilities.even

.. autofunction:: icesat2_toolkit.utilities.ceil

.. autofunction:: icesat2_toolkit.utilities.copy

.. autofunction:: icesat2_toolkit.utilities.check_ftp_connection

.. autofunction:: icesat2_toolkit.utilities.ftp_list

.. autofunction:: icesat2_toolkit.utilities.from_ftp

.. autofunction:: icesat2_toolkit.utilities.check_connection

.. autofunction:: icesat2_toolkit.utilities.http_list

.. autofunction:: icesat2_toolkit.utilities.from_http

.. autofunction:: icesat2_toolkit.utilities.attempt_login

.. autofunction:: icesat2_toolkit.utilities.build_opener

.. autofunction:: icesat2_toolkit.utilities.check_credentials

.. autofunction:: icesat2_toolkit.utilities.nsidc_list

.. autofunction:: icesat2_toolkit.utilities.from_nsidc

.. autofunction:: icesat2_toolkit.utilities.query_release

.. autofunction:: icesat2_toolkit.utilities.cycles

.. autofunction:: icesat2_toolkit.utilities.tracks

.. autofunction:: icesat2_toolkit.utilities.granules

.. autofunction:: icesat2_toolkit.utilities.regions

.. autofunction:: icesat2_toolkit.utilities.resolutions

.. autofunction:: icesat2_toolkit.utilities.readable_granules

.. autofunction:: icesat2_toolkit.utilities.cmr_filter_json

.. autofunction:: icesat2_toolkit.utilities.cmr
