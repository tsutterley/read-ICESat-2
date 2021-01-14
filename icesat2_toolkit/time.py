#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (01/2021)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 01/2021: added ftp connection checks
        merged with convert_julian and convert_calendar_decimal
        added calendar_days routine to get number of days per month
    Updated 08/2020: added NASA Earthdata routines for downloading from NSIDC
    Written 07/2020
"""
import os
import re
import datetime
import numpy as np
import icesat2_toolkit.utilities

#-- PURPOSE: gets the number of days per month for a given year
def calendar_days(year):
    """
    Calculates the number of days per month for a given year

    Arguments
    ---------
    year: calendar year

    Returns
    -------
    dpm: number of days for each month
    """
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31],dtype=np.float)
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float)
    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    #-- find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return dpm_leap
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return dpm_stnd

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
    epoch=(1992,1,1,0,0,0), scale=1.0):
    """
    Calculate the time in time units since epoch from calendar dates

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
    scale: scaling factor for converting time to output units

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
    return scale*np.array(MJD - delta_time_epochs/86400.0,dtype=np.float)

#-- PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(year, month, day=None, hour=None, minute=None,
    second=None, DofY=None):
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Dershowitz, N. and E.M. Reingold. 2008.  Calendrical Calculations.
        Cambridge: Cambridge University Press.

    Arguments
    ---------
    year: calendar year
    month: calendar month

    Keyword arguments
    -----------------
    day: day of the month
    hour: hour of the day
    minute: minute of the hour
    second: second of the minute
    DofY: day of the year (January 1 = 1)

    Returns
    -------
    t_date: date in decimal-year format
    """

    #-- number of dates
    n_dates = len(np.atleast_1d(year))

    #-- create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    #-- day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    #-- remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    #-- create output date variable
    t_date = np.zeros((n_dates))

    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap=np.array([31,29,31,30,31,30,31,31,30,31,30,31], dtype=np.float)
    dpm_stnd=np.array([31,28,31,30,31,30,31,31,30,31,30,31], dtype=np.float)

    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    #-- find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    #-- calculate the day of the year
    if DofY is not None:
        #-- if entered directly as an input
        #-- remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        #-- use calendar month and day of the month to calculate day of the year
        #-- month minus 1: January = 0, February = 1, etc (indice of month)
        #-- in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=np.int) - 1

        #-- day of month
        if day is not None:
            #-- remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            #-- if not entering days as an input
            #-- will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        #-- create matrix with the lower half = 1
        #-- this matrix will be used in a matrix multiplication
        #-- to calculate the total number of days for prior months
        #-- the -1 will make the diagonal == 0
        #-- i.e. first row == all zeros and the
        #-- last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        #-- using a dot product to calculate total number of days
        #-- for the months before the input date
        #-- basically is sum(i*dpm)
        #-- where i is 1 for all months < the month of interest
        #-- and i is 0 for all months >= the month of interest
        #-- month of interest is zero as the exact days will be
        #-- used to calculate the date

        #-- calculate the day of the year for leap and standard
        #-- use total days of all months before date
        #-- and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    #-- hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    #-- minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    #-- second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    #-- calculate decimal date
    #-- convert hours, minutes and seconds into days
    #-- convert calculated fractional days into decimal fractions of the year
    #-- Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    #-- Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

#-- PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD, ASTYPE=None, FORMAT='dict'):
    """
    Converts from Julian day to calendar date and time

    Translated from caldat in "Numerical Recipes in C", by William H. Press,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, 1988 (second printing).
    Hatcher, D. A., "Simple Formulae for Julian Day Numbers and Calendar Dates",
        Quarterly Journal of the Royal Astronomical Society, 25(1), 1984.


    Arguments
    ---------
    JD: Julian Day (days since 01-01-4713 BCE at 12:00:00)

    Keyword arguments
    -----------------
    ASTYPE: convert output to variable type
    FORMAT: format of output variables
        'dict': dictionary with variable keys
        'tuple': tuple with variable order YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
        'zip': aggregated variable sets

    Returns
    -------
    year: calendar year
    month: calendar month
    day: day of the month
    hour: hour of the day
    minute: minute of the hour
    second: second of the minute
    """

    #-- convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        SINGLE_VALUE = True
    else:
        SINGLE_VALUE = False

    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    #-- calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    #-- calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    #-- calculate day, month, year and hour
    DAY = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    MONTH = F - 1.0 - 12.0*np.floor(F/14.0)
    YEAR = D - 4715.0 - np.floor((7.0+MONTH)/10.0)
    HOUR = np.floor(24.0*(JD + 0.5 - JDO))
    #-- calculate minute and second
    G = (JD + 0.5 - JDO) - HOUR/24.0
    MINUTE = np.floor(G*1440.0)
    SECOND = (G - MINUTE/1440.0) * 86400.0

    #-- convert all variables to output type (from float)
    if ASTYPE is not None:
        YEAR = YEAR.astype(ASTYPE)
        MONTH = MONTH.astype(ASTYPE)
        DAY = DAY.astype(ASTYPE)
        HOUR = HOUR.astype(ASTYPE)
        MINUTE = MINUTE.astype(ASTYPE)
        SECOND = SECOND.astype(ASTYPE)

    #-- if only a single value was imported initially: remove singleton dims
    if SINGLE_VALUE:
        YEAR = YEAR.item(0)
        MONTH = MONTH.item(0)
        DAY = DAY.item(0)
        HOUR = HOUR.item(0)
        MINUTE = MINUTE.item(0)
        SECOND = SECOND.item(0)

    #-- return date variables in output format (default python dictionary)
    if (FORMAT == 'dict'):
        return dict(year=YEAR, month=MONTH, day=DAY,
            hour=HOUR, minute=MINUTE, second=SECOND)
    elif (FORMAT == 'tuple'):
        return (YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
    elif (FORMAT == 'zip'):
        return zip(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)

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
            indices = np.nonzero(GPS_Time >= leap)
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

#-- PURPOSE: connects to servers and downloads leap second files
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
        icesat2_toolkit.utilities.check_ftp_connection(HOST[0])
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
