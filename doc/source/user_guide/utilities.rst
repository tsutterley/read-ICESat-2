============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from NSIDC via https when Earthdata credentials are supplied
 - Checks MD5 hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/utilities.py


General Methods
===============


.. method:: icesat2_toolkit.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        `relpath`: local relative path as list or string


.. method:: icesat2_toolkit.utilities.get_hash(local)

    Get the MD5 hash value from a local file

    Arguments:

        `local`: path to file


.. method:: icesat2_toolkit.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        `time_string`: formatted time string to parse

    Keyword arguments:

        `format`: format for input time string


.. method:: icesat2_toolkit.utilities.copy(source, destination, verbose=False, move=False)

    Copy or move a file with all system information

    Arguments:

        `source`: source file

        `destination`: copied destination file

    Keyword arguments:

        `verbose`: print file transfer information

        `move`: remove the source file


.. method:: icesat2_toolkit.utilities.ftp_list(HOST,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `output`: list of items in a directory

        `mtimes`: list of last modification times for items in the directory


.. method:: icesat2_toolkit.utilities.from_ftp(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: icesat2_toolkit.utilities.from_http(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: icesat2_toolkit.utilities.build_opener(username,password,context=ssl.SSLContext(),password_manager=True,get_ca_certs=False,redirect=False,authorization_header=True,urs=None)

    build urllib opener for NASA Earthdata with supplied credentials

    Arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

    Keyword arguments:

        `context`: SSL context for opener object

        `password_manager`: create password manager context using default realm

        `get_ca_certs`: get list of loaded “certification authority” certificates

        `redirect`: create redirect handler object

        `authorization_header`: add base64 encoded authorization header to opener

        `urs`: Earthdata login URS 3 host


.. method:: icesat2_toolkit.utilities.check_credentials()

    Check that entered NASA Earthdata credentials are valid


.. method:: icesat2_toolkit.utilities.nsidc_list(HOST,username=None,password=None,build=True,timeout=None,parser=None,pattern='',sort=False)

    Download a file from a NSIDC https server

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `build`: Build opener and check Earthdata credentials

        `timeout`: timeout in seconds for blocking operations

        `parser`: HTML parser for lxml

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `colnames`: list of column names in a directory

        `collastmod`: list of last modification times for items in the directory


.. method:: icesat2_toolkit.utilities.from_nsidc(HOST,username=None,password=None,build=True,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a NSIDC https server

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `build`: Build opener and check Earthdata credentials

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file
