#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (02/2022)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

UPDATE HISTORY:
    Updated 02/2022: add NASA Common Metadata Repository (CMR) queries
        added generic list from Apache http server. verify host inputs
    Updated 10/2021: using python logging for handling verbose output
        add parser for converting file lines to arguments
    Updated 08/2021: NSIDC no longer requires authentication headers
    Updated 07/2021: return Earthdata opener from build function
    Updated 03/2021: added sha1 option for retrieving file hashes
    Updated 01/2021: added username and password to ftp functions
        added ftp connection check
    Updated 12/2020: added file object keyword for downloads if verbose
        add url split function for creating url location lists
    Updated 11/2020: nsidc_list and from_nsidc will output error strings
        normalize source and destination paths in copy
    Updated 09/2020: copy from http and https to bytesIO object in chunks
        use netrc credentials if not entered from NSIDC functions
        generalize build opener function for different Earthdata instances
    Updated 08/2020: add Earthdata opener, login and download functions
    Written 08/2020
"""
from __future__ import print_function

import sys
import os
import re
import io
import ssl
import json
import netrc
import ftplib
import shutil
import base64
import socket
import inspect
import hashlib
import logging
import datetime
import warnings
import posixpath
import lxml.etree
import calendar,time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
    import urllib.request as urllib2

def get_data_path(relpath):
    """
    Get the absolute path within a package from a relative path

    Arguments
    ---------
    relpath: relative path
    """
    #-- current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    if isinstance(relpath,list):
        #-- use *splat operator to extract from list
        return os.path.join(filepath,*relpath)
    elif isinstance(relpath,str):
        return os.path.join(filepath,relpath)

#-- PURPOSE: get the hash value of a file
def get_hash(local, algorithm='MD5'):
    """
    Get the hash value from a local file or BytesIO object

    Arguments
    ---------
    local: BytesIO object or path to file

    Keyword Arguments
    -----------------
    algorithm: hashing algorithm for checksum validation
        MD5: Message Digest
        sha1: Secure Hash Algorithm
    """
    #-- check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        if (algorithm == 'MD5'):
            return hashlib.md5(local.getvalue()).hexdigest()
        elif (algorithm == 'sha1'):
            return hashlib.sha1(local.getvalue()).hexdigest()
    elif os.access(os.path.expanduser(local),os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(os.path.expanduser(local), 'rb') as local_buffer:
            #-- generate checksum hash for a given type
            if (algorithm == 'MD5'):
                return hashlib.md5(local_buffer.read()).hexdigest()
            elif (algorithm == 'sha1'):
                return hashlib.sha1(local_buffer.read()).hexdigest()
    else:
        return ''

#-- PURPOSE: recursively split a url path
def url_split(s):
    """
    Recursively split a url path into a list

    Arguments
    ---------
    s: url string
    """
    head, tail = posixpath.split(s)
    if head in ('http:','https:','ftp:'):
        return s,
    elif head in ('', posixpath.sep):
        return tail,
    return url_split(head) + (tail,)

#-- PURPOSE: convert file lines to arguments
def convert_arg_line_to_args(arg_line):
    """
    Convert file lines to arguments

    Arguments
    ---------
    arg_line: line string containing a single argument and/or comments
    """
    #-- remove commented lines and after argument comments
    for arg in re.sub(r'\#(.*?)$',r'',arg_line).split():
        if not arg.strip():
            continue
        yield arg

#-- PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S'):
    """
    Get the Unix timestamp value for a formatted date string

    Arguments
    ---------
    time_string: formatted time string to parse

    Keyword arguments
    -----------------
    format: format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        return None
    else:
        return calendar.timegm(parsed_time)

#-- PURPOSE: rounds a number to an even number less than or equal to original
def even(value):
    """
    Rounds a number to an even number less than or equal to original

    Arguments
    ---------
    value: number to be rounded
    """
    return 2*int(value//2)

#-- PURPOSE: rounds a number upward to its nearest integer
def ceil(value):
    """
    Rounds a number upward to its nearest integer

    Arguments
    ---------
    value: number to be rounded upward
    """
    return -int(-value//1)

#-- PURPOSE: make a copy of a file with all system information
def copy(source, destination, verbose=False, move=False):
    """
    Copy or move a file with all system information

    Arguments
    ---------
    source: source file
    destination: copied destination file

    Keyword arguments
    -----------------
    verbose: print file transfer information
    move: remove the source file
    """
    source = os.path.abspath(os.path.expanduser(source))
    destination = os.path.abspath(os.path.expanduser(destination))
    logging.info('{0} -->\n\t{1}'.format(source,destination))
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    if move:
        os.remove(source)

#-- PURPOSE: check ftp connection
def check_ftp_connection(HOST,username=None,password=None):
    """
    Check internet connection with ftp host

    Arguments
    ---------
    HOST: remote ftp host

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    """
    #-- attempt to connect to ftp host
    try:
        f = ftplib.FTP(HOST)
        f.login(username, password)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        return True

#-- PURPOSE: list a directory on a ftp host
def ftp_list(HOST,username=None,password=None,timeout=None,
    basename=False,pattern=None,sort=False):
    """
    List a directory on a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    timeout: timeout in seconds for blocking operations
    basename: return the file or directory basename instead of the full path
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    output: list of items in a directory
    mtimes: list of last modification times for items in the directory
    """
    #-- verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    else:
        ftp.login(username,password)
        #-- list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        #-- get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        #-- iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                #-- try sending modification time command
                mdtm = ftp.sendcmd('MDTM {0}'.format(f))
            except ftplib.error_perm:
                #-- directories will return with an error
                pass
            else:
                #-- convert the modification time into unix time
                mtimes[i] = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        #-- reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            #-- reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            #-- sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- close the ftp connection
        ftp.close()
        #-- return the list of items and last modified times
        return (output,mtimes)

#-- PURPOSE: download a file from a ftp host
def from_ftp(HOST,username=None,password=None,timeout=None,local=None,
    hash='',chunk=8192,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp username
    password: ftp password
    timeout: timeout in seconds for blocking operations
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try downloading from ftp
    try:
        #-- try to connect to ftp host
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    else:
        ftp.login(username,password)
        #-- remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary('RETR {0}'.format(ftp_remote_path),
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd('MDTM {0}'.format(ftp_remote_path))
        remote_mtime = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
            #-- keep remote modification time of file and local access time
            os.utime(local, (os.stat(local).st_atime, remote_mtime))
        #-- close the ftp connection
        ftp.close()
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: check internet connection
def check_connection(HOST):
    """
    Check internet connection with http host

    Arguments
    ---------
    HOST: remote http host
    """
    #-- attempt to connect to http host
    try:
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: list a directory on an Apache http Server
def http_list(HOST,timeout=None,context=ssl.SSLContext(),
    parser=lxml.etree.HTMLParser(),format='%Y-%m-%d %H:%M',
    pattern='',sort=False):
    """
    List a directory on an Apache http Server

    Arguments
    ---------
    HOST: remote http host path split as list

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    context: SSL context for url opener object
    parser: HTML parser for lxml
    format: format for input time string
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    colnames: list of column names in a directory
    collastmod: list of last modification times for items in the directory
    colerror: notification for list error
    """
    #-- verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try listing from http
    try:
        #-- Create and submit request.
        request=urllib2.Request(posixpath.join(*HOST))
        response=urllib2.urlopen(request,timeout=timeout,context=context)
    except (urllib2.HTTPError, urllib2.URLError) as e:
        colerror = 'List error from {0}'.format(posixpath.join(*HOST))
        return (False,False,colerror)
    else:
        #-- read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response,parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        #-- get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            #-- reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            #-- sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- return the list of column names and last modified times
        return (colnames,collastmod,None)

#-- PURPOSE: download a file from a http host
def from_http(HOST,timeout=None,context=ssl.SSLContext(),local=None,hash='',
    chunk=16384,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a http host

    Arguments
    ---------
    HOST: remote http host path split as list

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    context: SSL context for url opener object
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try downloading from http
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout,context=context)
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: "login" to NASA Earthdata with supplied credentials
def build_opener(username, password, context=ssl.SSLContext(),
    password_manager=True, get_ca_certs=False, redirect=False,
    authorization_header=False, urs='https://urs.earthdata.nasa.gov'):
    """
    build urllib opener for NASA Earthdata with supplied credentials

    Arguments
    ---------
    username: NASA Earthdata username
    password: NASA Earthdata password

    Keyword arguments
    -----------------
    context: SSL context for opener object
    password_manager: create password manager context using default realm
    get_ca_certs: get list of loaded “certification authority” certificates
    redirect: create redirect handler object
    authorization_header: add base64 encoded authorization header to opener
    urs: Earthdata login URS 3 host
    """
    #-- https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    #-- create a password manager
    if password_manager:
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        #-- Add the username and password for NASA Earthdata Login system
        password_mgr.add_password(None,urs,username,password)
        handler.append(urllib2.HTTPBasicAuthHandler(password_mgr))
    #-- Create cookie jar for storing cookies. This is used to store and return
    #-- the session cookie given to use by the data server (otherwise will just
    #-- keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    #-- SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    handler.append(urllib2.HTTPSHandler(context=context))
    #-- redirect handler
    if redirect:
        handler.append(urllib2.HTTPRedirectHandler())
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    #-- Encode username/password for request authorization headers
    #-- add Authorization header to opener
    if authorization_header:
        b64 = base64.b64encode('{0}:{1}'.format(username,password).encode())
        opener.addheaders = [("Authorization","Basic {0}".format(b64.decode()))]
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.
    return opener

#-- PURPOSE: check that entered NASA Earthdata credentials are valid
def check_credentials():
    """
    Check that entered NASA Earthdata credentials are valid
    """
    try:
        remote_path = posixpath.join('https://n5eil01u.ecs.nsidc.org','ATLAS')
        request = urllib2.Request(url=remote_path)
        response = urllib2.urlopen(request, timeout=20)
    except urllib2.HTTPError:
        raise RuntimeError('Check your NASA Earthdata credentials')
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: list a directory on NSIDC https server
def nsidc_list(HOST,username=None,password=None,build=True,timeout=None,
    parser=lxml.etree.HTMLParser(),pattern='',sort=False):
    """
    List a directory on NSIDC

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    username: NASA Earthdata username
    password: NASA Earthdata password
    build: Build opener and check NASA Earthdata credentials
    timeout: timeout in seconds for blocking operations
    parser: HTML parser for lxml
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    colnames: list of column names in a directory
    collastmod: list of last modification times for items in the directory
    colerror: notification for list error
    """
    #-- use netrc credentials
    if build and not (username or password):
        urs = 'urs.earthdata.nasa.gov'
        username,login,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try listing from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        tree = lxml.etree.parse(urllib2.urlopen(request,timeout=timeout),parser)
    except (urllib2.HTTPError, urllib2.URLError) as e:
        colerror = 'List error from {0}'.format(posixpath.join(*HOST))
        return (False,False,colerror)
    else:
        #-- read and parse request for files (column names and modified times)
        colnames = tree.xpath('//td[@class="indexcolname"]//a/@href')
        #-- get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format='%Y-%m-%d %H:%M')
            for i in tree.xpath('//td[@class="indexcollastmod"]/text()')]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern,f)]
            #-- reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            #-- sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        #-- return the list of column names and last modified times
        return (colnames,collastmod,None)

#-- PURPOSE: download a file from a NSIDC https server
def from_nsidc(HOST,username=None,password=None,build=True,timeout=None,
    local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775):
    """
    Download a file from a NSIDC https server

    Arguments
    ---------
    HOST: remote https host path split as list

    Keyword arguments
    -----------------
    username: NASA Earthdata username
    password: NASA Earthdata password
    build: Build opener and check NASA Earthdata credentials
    timeout: timeout in seconds for blocking operations
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    fid: open file object to print if verbose
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    response_error: notification for response error
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- use netrc credentials
    if build and not (username or password):
        urs = 'urs.earthdata.nasa.gov'
        username,login,password = netrc.netrc().authenticators(urs)
    #-- build urllib2 opener and check credentials
    if build:
        #-- build urllib2 opener with credentials
        build_opener(username, password)
        #-- check credentials
        check_credentials()
    #-- verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    #-- try downloading from https
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout)
    except:
        response_error = 'Download error from {0}'.format(posixpath.join(*HOST))
        return (False,response_error)
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- convert to absolute path
            local = os.path.abspath(local)
            #-- create directory if non-existent
            if not os.access(os.path.dirname(local), os.F_OK):
                os.makedirs(os.path.dirname(local), mode)
            #-- print file information
            args = (posixpath.join(*HOST),local)
            logging.info('{0} -->\n\t{1}'.format(*args), file=fid)
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return (remote_buffer,None)

#-- PURPOSE: build formatted query string for ICESat-2 release
def query_release(release):
    """
    Build formatted query string for ICESat-2 release

    Arguments
    ---------
    release: ICESat-2 data release to query

    Returns
    -------
    query_params: formatted string for CMR queries
    """
    if release is None:
        return ''
    #-- maximum length of version in CMR queries
    desired_pad_length = 3
    if len(str(release)) > desired_pad_length:
        raise RuntimeError('Release string too long: "{0}"'.format(release))
    #-- Strip off any leading zeros
    release = str(release).lstrip('0')
    query_params = ''
    while len(release) <= desired_pad_length:
        padded_release = release.zfill(desired_pad_length)
        query_params += '&version={0}'.format(padded_release)
        desired_pad_length -= 1
    return query_params

#-- PURPOSE: check if the submitted cycles are valid
def cycles(cycle):
    """
    Check if the submitted cycles are valid

    Arguments
    ---------
    cycle: ICESat-2 orbital cycle
    """
    #-- string length of cycles in granules
    cycle_length = 2
    #-- number of GPS seconds between the GPS epoch and ATLAS SDP epoch
    atlas_sdp_gps_epoch = 1198800018.0
    #-- number of GPS seconds since the GPS epoch for first ATLAS data point
    atlas_gps_start_time = atlas_sdp_gps_epoch + 24710205.39202261
    epoch1 = datetime.datetime(1980, 1, 6, 0, 0, 0)
    epoch2 = datetime.datetime(1970, 1, 1, 0, 0, 0)
    #-- get the total number of seconds since the start of ATLAS and now
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    atlas_UNIX_start_time = atlas_gps_start_time - delta_time_epochs
    present_time = datetime.datetime.now().timestamp()
    #-- divide total time by cycle length to get the maximum number of orbital cycles
    ncycles = ceil((present_time - atlas_UNIX_start_time) / (86400 * 91))
    all_cycles = [str(c + 1).zfill(cycle_length) for c in range(ncycles)]
    if cycle is None:
        return ["??"]
    else:
        if isinstance(cycle, (str,int)):
            assert int(cycle) > 0, "Cycle number must be positive"
            cycle_list = [str(cycle).zfill(cycle_length)]
        elif isinstance(cycle, list):
            cycle_list = []
            for c in cycle:
                assert int(c) > 0, "Cycle number must be positive"
                cycle_list.append(str(c).zfill(cycle_length))
        else:
            raise TypeError("Please enter the cycle number as a list or string")
        #-- check if user-entered cycle is outside of currently available range
        if not set(all_cycles) & set(cycle_list):
            warnings.filterwarnings("always")
            warnings.warn("Listed cycle is not presently available")
        return cycle_list

#-- PURPOSE: check if the submitted RGTs are valid
def tracks(track):
    """
    Check if the submitted RGTs are valid

    Arguments
    ---------
    track: ICESat-2 reference ground track (RGT)
    """
    #-- string length of RGTs in granules
    track_length = 4
    #-- total number of ICESat-2 satellite RGTs is 1387
    all_tracks = [str(tr + 1).zfill(track_length) for tr in range(1387)]
    if track is None:
        return ["????"]
    else:
        if isinstance(track, (str,int)):
            assert int(track) > 0, "Reference Ground Track must be positive"
            track_list = [str(track).zfill(track_length)]
        elif isinstance(track, list):
            track_list = []
            for t in track:
                assert int(t) > 0, "Reference Ground Track must be positive"
                track_list.append(str(t).zfill(track_length))
        else:
            raise TypeError(
                "Reference Ground Track as a list or string"
            )
        #-- check if user-entered RGT is outside of the valid range
        if not set(all_tracks) & set(track_list):
            warnings.filterwarnings("always")
            warnings.warn("Listed Reference Ground Track is not available")
        return track_list

#-- PURPOSE: check if the submitted granule regions are valid
def granules(granule):
    """
    Check if the submitted granule regions are valid

    Arguments
    ---------
    granule: ICESat-2 granule region
    """
    #-- string length of granule regions in granule files
    granule_length = 2
    #-- total number of ICESat-2 granule regions is 14
    all_granules = [str(g).zfill(granule_length) for g in range(1,15)]
    if granule is None:
        return ["??"]
    else:
        if isinstance(granule, (str,int)):
            assert int(granule) > 0, "Cycle region must be positive"
            granule_list = [str(granule).zfill(granule_length)]
        elif isinstance(granule, list):
            granule_list = []
            for g in granule:
                assert int(g) > 0, "Granule region must be positive"
                granule_list.append(str(g).zfill(granule_length))
        else:
            raise TypeError("Please enter the cycle number as a list or string")
        #-- check if user-entered granule is outside of currently available range
        if not set(all_granules) & set(granule_list):
            warnings.filterwarnings("always")
            warnings.warn("Listed cycle is not presently available")
        return granule_list

def readable_granules(product, **kwargs):
    """
    Create list of readable granule names for CMR queries

    Arguments
    ---------
    product: ICESat-2 data product

    Keyword arguments
    -----------------
    cycles: List of 91-day orbital cycle strings to query
    tracks: List of Reference Ground Track (RGT) strings to query
    granules: List of ICESat-2 granule region strings to query

    Returns
    -------
    list of readable granule names for CMR queries
    """
    #-- default keyword arguments
    kwargs.setdefault("cycles", None)
    kwargs.setdefault("tracks", None)
    kwargs.setdefault("granules", None)
    #-- list of readable granule names
    readable_granule_list = []
    #-- for each available cycle of interest
    for c in cycles(kwargs["cycles"]):
        #-- for each available track of interest
        for t in tracks(kwargs["tracks"]):
            #-- for each available granule region of interest
            for g in granules(kwargs["granules"]):
                #-- use single character wildcards "?" for date strings,
                #-- sea ice product hemispheres, and any unset parameters
                if product in ("ATL07", "ATL10", "ATL20", "ATL21"):
                    args = (product, 14 * "?", t, c, g)
                    pattern = "{0}-??_{1}_{2}{3}{4}_*"
                elif product in ("ATL11",):
                    args = (product, t, g)
                    pattern = "{0}_{1}{2}_*"
                else:
                    args = (product, 14 * "?", t, c, g)
                    pattern = "{0}_{1}_{2}{3}{4}_*"
                #-- append the granule pattern
                readable_granule_list.append(pattern.format(*args))
    #-- return readable granules list
    return readable_granule_list

#-- PURPOSE: filter the CMR json response for desired data files
def cmr_filter_json(search_results, request_type="application/x-hdfeos"):
    """
    Filter the CMR json response for desired data files

    Arguments
    ---------
    search_results: json response from CMR query

    Keyword arguments
    -----------------
    request_type: data type for reducing CMR query

    Returns
    -------
    producer_granule_ids: list of ICESat-2 granules
    granule_urls: list of ICESat-2 granule urls from NSIDC
    """
    #-- output list of granule ids and urls
    producer_granule_ids = []
    granule_urls = []
    #-- check that there are urls for request
    if ('feed' not in search_results) or ('entry' not in search_results['feed']):
        return (producer_granule_ids,granule_urls)
    #-- iterate over references and get cmr location
    for entry in search_results['feed']['entry']:
        producer_granule_ids.append(entry['producer_granule_id'])
        for link in entry['links']:
            if (link['type'] == request_type):
                granule_urls.append(link['href'])
                break
    #-- return the list of urls and granule ids
    return (producer_granule_ids,granule_urls)

#-- PURPOSE: cmr queries for orbital parameters
def cmr(product=None, release=None, cycles=None, tracks=None,
    granules=None, verbose=False, fid=sys.stdout):
    """
    Query the NASA Common Metadata Repository (CMR) for ICESat-2 data

    Keyword arguments
    -----------------
    product: ICESat-2 data product to query
    release: ICESat-2 data release to query
    cycles: List of 91-day orbital cycle strings to query
    tracks: List of Reference Ground Track (RGT) strings to query
    granules: List of ICESat-2 granule region strings to query
    verbose: print file transfer information
    fid: open file object to print if verbose

    Returns
    -------
    producer_granule_ids: list of ICESat-2 granules
    granule_urls: list of ICESat-2 granule urls from NSIDC
    """
    #-- create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    #-- build urllib2 opener with SSL context
    build_opener(None, None, context=ssl.SSLContext(),
        password_manager=False)
    #-- build CMR query
    cmr_format = 'json'
    cmr_provider = 'NSIDC_ECS'
    cmr_page_size = 2000
    CMR_HOST = ['https://cmr.earthdata.nasa.gov','search',
        'granules.{0}'.format(cmr_format)]
    #-- build list of CMR query parameters
    CMR_KEYS = []
    CMR_KEYS.append('?provider={0}'.format(cmr_provider))
    CMR_KEYS.append('&sort_key[]=start_date')
    CMR_KEYS.append('&sort_key[]=producer_granule_id')
    CMR_KEYS.append('&scroll=true')
    CMR_KEYS.append('&page_size={0}'.format(cmr_page_size))
    #-- append product string
    CMR_KEYS.append('&short_name={0}'.format(product))
    #-- append release strings
    CMR_KEYS.append(query_release(release))
    #-- append keys for querying specific granules
    CMR_KEYS.append("&options[readable_granule_name][pattern]=true")
    CMR_KEYS.append("&options[spatial][or]=true")
    readable_granule_list = readable_granules(product,
        cycles=cycles, tracks=tracks, granules=granules)
    for gran in readable_granule_list:
        CMR_KEYS.append("&readable_granule_name[]={0}".format(gran))
    #-- full CMR query url
    cmr_query_url = "".join([posixpath.join(*CMR_HOST),*CMR_KEYS])
    logging.info('CMR request={0}'.format(cmr_query_url))
    #-- output list of granule names and urls
    producer_granule_ids = []
    granule_urls = []
    cmr_scroll_id = None
    while True:
        req = urllib2.Request(cmr_query_url)
        if cmr_scroll_id:
            req.add_header('cmr-scroll-id', cmr_scroll_id)
        response = urllib2.urlopen(req)
        #-- get scroll id for next iteration
        if not cmr_scroll_id:
            headers = {k.lower():v for k,v in dict(response.info()).items()}
            cmr_scroll_id = headers['cmr-scroll-id']
        #-- read the CMR search as JSON
        search_page = json.loads(response.read().decode('utf-8'))
        ids,urls = cmr_filter_json(search_page)
        if not urls:
            break
        #-- extend lists
        producer_granule_ids.extend(ids)
        granule_urls.extend(urls)
    #-- return the list of granule ids and urls
    return (producer_granule_ids, granule_urls)
