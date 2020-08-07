#!/usr/bin/env python
u"""
nsidc_icesat2_associated.py
Written by Tyler Sutterley (05/2020)

Program to acquire ICESat-2 datafiles from NSIDC server that is
    associated with an input file

https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-
    https-earthdata-login-enabled
http://www.voidspace.org.uk/python/articles/authentication.shtml#base64

Register with NASA Earthdata Login system:
https://urs.earthdata.nasa.gov

Add NSIDC_DATAPOOL_OPS to NASA Earthdata Applications
https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA

CALLING SEQUENCE:
    python nsidc_icesat2_associated.py --user=<username> --release=001
        --product=ATL03 ATL06_files
    where <username> is your NASA Earthdata username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user=X: username for NASA Earthdata Login
    -N X, --netrc=X: path to .netrc file for alternative authentication
    -D X, --directory: working data directory
    -P X, --product: associated product to download
        ATL03: Global Geolocated Photon Data
        ATL04: Normalized Relative Backscatter
        ATL06: Land Ice Height
        ATL07: Sea Ice Height
        ATL08: Land and Vegetation Height
        ATL09: Atmospheric Layer Characteristics
        ATL10: Sea Ice Freeboard
        ATL12: Ocean Surface Height
        ATL13: Inland Water Surface Height
    --auxiliary: Download ICESat-2 auxiliary files for each HDF5 file
    -M X, --mode=X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml

UPDATE HISTORY:
    Updated 05/2020: added option netrc to use alternative authentication
        adjust regular expression to allow syncing of ATL07 sea ice products
        adjust regular expression for auxiliary products
    Updated 09/2019: added ssl context to urlopen headers
    Written 08/2019
"""
from __future__ import print_function

import sys
import os
import re
import ssl
import netrc
import getopt
import shutil
import base64
import getpass
import builtins
import posixpath
import lxml.etree
import numpy as np
import calendar, time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    import urllib2
else:
    from http.cookiejar import CookieJar
    import urllib.request as urllib2

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to https host for NSIDC
    try:
        HOST = 'https://n5eil01u.ecs.nsidc.org'
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: download the ICESat-2 elevation data from NSIDC matching an file
def nsidc_icesat2_associated(file_list, PRODUCT, USER='', PASSWORD='',
    DIRECTORY=None, AUXILIARY=False, MODE=0o775):

    #-- https://docs.python.org/3/howto/urllib2.html#id5
    #-- create a password manager
    password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    #-- Add the username and password for NASA Earthdata Login system
    password_mgr.add_password(None, 'https://urs.earthdata.nasa.gov',
        USER, PASSWORD)
    #-- Encode username/password for request authorization headers
    base64_string = base64.b64encode('{0}:{1}'.format(USER,PASSWORD).encode())
    #-- compile HTML parser for lxml
    parser = lxml.etree.HTMLParser()
    #-- Create cookie jar for storing cookies. This is used to store and return
    #-- the session cookie given to use by the data server (otherwise will just
    #-- keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(
        urllib2.HTTPBasicAuthHandler(password_mgr),
        urllib2.HTTPSHandler(context=ssl.SSLContext()),
        urllib2.HTTPCookieProcessor(cookie_jar))
    #-- add Authorization header to opener
    authorization_header = "Basic {0}".format(base64_string.decode())
    opener.addheaders = [("Authorization", authorization_header)]
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.

    #-- remote https server for ICESat-2 Data
    HOST = 'https://n5eil01u.ecs.nsidc.org'
    #-- regular expression operator for extracting information from files
    rx = re.compile(r'(processed_)?(ATL\d{2})(-\d{2})?_(\d{4})(\d{2})(\d{2})'
        r'(\d{2})(\d{2})(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    #-- regular expression pattern for finding specific files
    regex_suffix = '(.*?)' if AUXILIARY else '(h5)'
    remote_regex_pattern = (r'{0}(-\d{{2}})?_(\d{{4}})(\d{{2}})(\d{{2}})'
        r'(\d{{2}})(\d{{2}})(\d{{2}})_({1})({2})({3})_({4})_({5})(.*?).{5}$')

    #-- for each input file
    for f in file_list:
        #-- print filename
        print(os.path.expanduser(f))

        #-- output data directory
        if DIRECTORY is None:
            #-- if not setting directory: use directory of file
            local_dir = os.path.dirname(os.path.expanduser(f))
        else:
            local_dir = os.path.expanduser(DIRECTORY)

        #-- extract parameters from ICESat-2 ATLAS HDF5 file name
        SUB,PRD,HEM,YY,MM,DD,HH,MN,SS,TRK,CYC,GRN,RL,VRS,AUX=rx.findall(f).pop()
        #-- set subdirectories for full remote directory (* splat operator)
        sd=['ATLAS','{0}.{1}'.format(PRODUCT,RL),'{0}.{1}.{2}'.format(YY,MM,DD)]
        d = posixpath.join(HOST,*sd)

        #-- compile regular expression operator for file parameters
        args = (PRODUCT,TRK,CYC,GRN,RL,VRS,regex_suffix)
        R1 = re.compile(remote_regex_pattern.format(*args), re.VERBOSE)

        #-- find ICESat-2 data file
        req = urllib2.Request(url=posixpath.join(d))
        #-- read and parse request for remote files (columns and dates)
        tree = lxml.etree.parse(urllib2.urlopen(req), parser)
        colnames = tree.xpath('//td[@class="indexcolname"]//a/@href')
        collastmod = tree.xpath('//td[@class="indexcollastmod"]/text()')
        #-- find matching files (for granule, release, version, track)
        remote_file_lines=[i for i,f in enumerate(colnames) if R1.match(f)]
        #-- sync each ICESat-2 data file
        for i in remote_file_lines:
            #-- remote and local versions of the file
            remote_file = posixpath.join(d,colnames[i])
            local_file = os.path.join(local_dir,colnames[i])
            #-- get last modified date and convert into unix time
            LMD = time.strptime(collastmod[i].rstrip(),'%Y-%m-%d %H:%M')
            remote_mtime = calendar.timegm(LMD)
            #-- sync ICESat-2 files with NSIDC server
            http_pull_file(remote_file, remote_mtime, local_file, MODE)
        #-- close request
        req = None

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(remote_file,remote_mtime,local_file,MODE):
    #-- Printing files transferred
    print('{0} -->\n\t{1}\n'.format(remote_file,local_file))
    #-- Create and submit request. There are a wide range of exceptions
    #-- that can be thrown here, including HTTPError and URLError.
    request = urllib2.Request(remote_file)
    response = urllib2.urlopen(request)
    #-- chunked transfer encoding size
    CHUNK = 16 * 1024
    #-- copy contents to local file using chunked transfer encoding
    #-- transfer should work properly with ascii and binary data formats
    with open(local_file, 'wb') as f:
        shutil.copyfileobj(response, f, CHUNK)
    #-- keep remote modification time of file and local access time
    os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
    os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -U X, --user=X\t\tUsername for NASA Earthdata Login')
    print(' -N X, --netrc=X\t\tPath to .netrc file for authentication')
    print(' -D X, --directory=X\tOutput data directory')
    print(' -P X, --product=X\tICESat-2 product to download')
    print(' --auxiliary\t\tSync ICESat-2 auxiliary files for each HDF5 file')
    print(' -M X, --mode=X\t\tPermission mode of files downloaded\n')

#-- Main program that calls nsidc_icesat2_associated()
def main():
    #-- Read the system arguments listed after the program
    long_options=['help','user=','netrc=','directory=','product=','auxiliary',
        'mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hU:N:D:P:AM:',long_options)

    #-- command line parameters
    USER = ''
    NETRC = None
    DIRECTORY = None
    PRODUCT = None
    AUXILIARY = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-U","--user"):
            USER = arg
        elif opt in ("-N","--netrc"):
            NETRC = os.path.expanduser(arg)
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--product"):
            PRODUCT = arg
        elif opt in ("--auxiliary"):
            AUXILIARY = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- Pre-ICESat-2 and IceBridge Products
    PROD = {}
    PROD['ATL03'] = 'Global Geolocated Photon Data'
    PROD['ATL04'] = 'Normalized Relative Backscatter'
    PROD['ATL06'] = 'Land Ice Height'
    PROD['ATL07'] = 'Sea Ice Height'
    PROD['ATL08'] = 'Land and Vegetation Height'
    PROD['ATL09'] = 'Atmospheric Layer Characteristics'
    PROD['ATL10'] = 'Sea Ice Freeboard'
    PROD['ATL12'] = 'Ocean Surface Height'
    PROD['ATL13'] = 'Inland Water Surface Height'

    #-- check that each data product entered was correctly typed
    keys = ','.join(sorted([key for key in PROD.keys()]))
    if PRODUCT not in PROD.keys():
        raise IOError('Incorrect Data Product Entered ({0})'.format(keys))

    #-- NASA Earthdata hostname
    HOST = 'urs.earthdata.nasa.gov'
    #-- get authentication
    if not USER and not NETRC:
        #-- check that NASA Earthdata credentials were entered
        USER = builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))
    elif NETRC:
        USER,LOGIN,PASSWORD = netrc.netrc(NETRC).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))

    #-- check internet connection before attempting to run program
    if check_connection():
        nsidc_icesat2_associated(arglist, PRODUCT, USER=USER, PASSWORD=PASSWORD,
            DIRECTORY=DIRECTORY, AUXILIARY=AUXILIARY, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
