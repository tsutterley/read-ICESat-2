#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (09/2020)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    lxml: processing XML and HTML in Python (https://pypi.python.org/pypi/lxml)

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_calendar_decimal.py: converts from calendar dates into decimal years
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 08/2020: added NASA Earthdata routines for downloading from CDDIS
    Written 07/2020
"""
import os
import re
import netrc
import datetime
import numpy as np
import icesat2_toolkit.convert_julian
import icesat2_toolkit.convert_calendar_decimal
import icesat2_toolkit.utilities

#-- PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0):
    """
    Convert delta time from seconds since epoch1 to time since epoch2

    Arguments
    ---------
    delta_time: seconds since epoch1

    Keyword arguments
    -----------------
    epoch1: epoch for input delta_time
    epoch2: epoch for output delta_time
    scale: scaling factor for converting time to output units
    """
    epoch1 = datetime.datetime(*epoch1)
    epoch2 = datetime.datetime(*epoch2)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

#-- PURPOSE: calculate the delta time from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0,
    epoch=(1992,1,1,0,0,0)):
    """
    Calculate the time in days since epoch from calendar dates

    Arguments
    ---------
    year: calendar month
    month: month of the year
    day: day of the month

    Keyword arguments
    -----------------
    hour: hour of the day
    minute: minute of the hour
    second: second of the minute
    epoch: epoch for output delta_time

    Returns
    -------
    delta_time: days since epoch
    """
    #-- calculate date in Modified Julian Days (MJD) from calendar date
    #-- MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    epoch1 = datetime.datetime(1858,11,17,0,0,0)
    epoch2 = datetime.datetime(*epoch)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- return the date in days since epoch
    return np.array(MJD - delta_time_epochs/86400.0,dtype=np.float)

#-- PURPOSE: Count number of leap seconds that have passed for each GPS time
def count_leap_seconds(GPS_Time):
    """
    Counts the number of leap seconds between a given GPS time and UTC

    Arguments
    ---------
    GPS_Time: seconds since January 6, 1980 at 00:00:00

    Returns
    -------
    n_leaps: number of elapsed leap seconds
    """
    #-- get the valid leap seconds
    leaps = get_leap_seconds()
    #-- number of leap seconds prior to GPS_Time
    n_leaps = np.zeros_like(GPS_Time,dtype=np.float)
    for i,leap in enumerate(leaps):
        count = np.count_nonzero(GPS_Time >= leap)
        if (count > 0):
            indices, = np.nonzero(GPS_Time >= leap)
            n_leaps[indices] += 1.0
    #-- return the number of leap seconds for converting to UTC
    return n_leaps

#-- PURPOSE: Define GPS leap seconds
def get_leap_seconds():
    """
    Gets a list of GPS times for when leap seconds occurred

    Returns
    -------
    GPS time (seconds since 1980-01-06T00:00:00) of leap seconds
    """
    FILE = icesat2_toolkit.utilities.get_data_path(['data','leap-seconds.list'])
    #-- find line with file expiration as delta time
    with open(FILE,'r') as fid:
        secs, = [re.findall(r'\d+',i).pop() for i in fid.read().splitlines()
            if re.match(r'^(?=#@)',i)]
    #-- check that leap seconds file is still valid
    expiry = datetime.datetime(1900,1,1) + datetime.timedelta(seconds=int(secs))
    today = datetime.datetime.now()
    update_leap_seconds() if (expiry < today) else None
    #-- get leap seconds
    leap_UTC,TAI_UTC=np.loadtxt(icesat2_toolkit.utilities.get_data_path(FILE)).T
    #-- TAI time is ahead of GPS by 19 seconds
    TAI_GPS = 19.0
    #-- convert leap second epochs from NTP to GPS
    #-- convert from time of 2nd leap second to time of 1st leap second
    leap_GPS = convert_delta_time(leap_UTC+TAI_UTC-TAI_GPS-1,
        epoch1=(1900,1,1,0,0,0), epoch2=(1980,1,6,0,0,0))
    #-- return the GPS times of leap second occurance
    return leap_GPS[leap_GPS >= 0].astype(np.float)

#-- PURPOSE: connects to servers and downloads leap esconds files
def update_leap_seconds(verbose=False, mode=0o775):
    """
    Connects to servers to download leap-seconds.list files from NIST servers
    https://www.nist.gov/pml/time-and-frequency-division/leap-seconds-faqs

    Servers and Mirrors
    ===================
    ftp://ftp.nist.gov/pub/time/leap-seconds.list
    https://www.ietf.org/timezones/data/leap-seconds.list

    Keyword arguments
    -----------------
    verbose: print file information about output file
    mode: permissions mode of output file
    """
    #-- local version of file
    FILE = 'leap-seconds.list'
    LOCAL = icesat2_toolkit.utilities.get_data_path(['data',FILE])
    HASH = icesat2_toolkit.utilities.get_hash(LOCAL)

    #-- try downloading from NIST ftp servers
    HOST = ['ftp.nist.gov','pub','time','iers',FILE]
    try:
        icesat2_toolkit.utilities.from_ftp(HOST, timeout=20, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except:
        pass
    else:
        return

    #-- try downloading from Internet Engineering Task Force (IETF) mirror
    REMOTE = ['https://www.ietf.org','timezones','data',FILE]
    try:
        icesat2_toolkit.utilities.from_http(REMOTE, timeout=5, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except:
        pass
    else:
        return

    #-- return exception that no server could be connected
    raise Exception('All Server Connection Error')
