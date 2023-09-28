#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (05/2023)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2023: allow epoch arguments to be numpy datetime64 or strings
        function to convert a string with time zone information to datetime
    Updated 04/2023: using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: output variables for some standard epochs
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 10/2022: added more time parsing for longer periods
        added encoding for reading leap seconds ascii files
    Updated 08/2022: output variables to unit conversion to seconds
        and the number of days per month for both leap and standard years
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 04/2021: updated NIST ftp server url for leap-seconds.list
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: parse date strings "time-units since yyyy-mm-dd hh:mm:ss"
        replaced numpy int/float to prevent deprecation warnings
    Updated 01/2021: added ftp connection checks
        merged with convert_julian and convert_calendar_decimal
        added calendar_days routine to get number of days per month
    Updated 08/2020: added NASA Earthdata routines for downloading from NSIDC
    Written 07/2020
"""
from __future__ import annotations

import re
import copy
import logging
import warnings
import datetime
import traceback
import numpy as np
import dateutil.parser
import icesat2_toolkit.utilities

# conversion factors between time units and seconds
_to_sec = {'microseconds': 1e-6, 'microsecond': 1e-6,
           'microsec': 1e-6, 'microsecs': 1e-6,
           'milliseconds': 1e-3, 'millisecond': 1e-3,
           'millisec': 1e-3, 'millisecs': 1e-3,
           'msec': 1e-3, 'msecs': 1e-3, 'ms': 1e-3,
           'seconds': 1.0, 'second': 1.0, 'sec': 1.0,
           'secs': 1.0, 's': 1.0,
           'minutes': 60.0, 'minute': 60.0,
           'min': 60.0, 'mins': 60.0,
           'hours': 3600.0, 'hour': 3600.0,
           'hr': 3600.0, 'hrs': 3600.0, 'h': 3600.0,
           'day': 86400.0, 'days': 86400.0, 'd': 86400.0}
# approximate conversions for longer periods
_to_sec['mon'] = 30.0 * 86400.0
_to_sec['month'] = 30.0 * 86400.0
_to_sec['months'] = 30.0 * 86400.0
_to_sec['common_year'] = 365.0 * 86400.0
_to_sec['common_years'] = 365.0 * 86400.0
_to_sec['year'] = 365.25 * 86400.0
_to_sec['years'] = 365.25 * 86400.0

# standard epochs
_mjd_epoch = (1858, 11, 17, 0, 0, 0)
_ntp_epoch = (1900, 1, 1, 0, 0, 0)
_unix_epoch = (1970, 1, 1, 0, 0, 0)
_gps_epoch = (1980, 1, 6, 0, 0, 0)
_tide_epoch = (1992, 1, 1, 0, 0, 0)
_j2000_epoch = (2000, 1, 1, 12, 0, 0)
_atlas_sdp_epoch = (2018, 1, 1, 0, 0, 0)

# PURPOSE: parse a date string and convert to a datetime object in UTC
def parse(date_string: str):
    """
    Parse a date string and convert to a naive ``datetime`` object in UTC

    Parameters
    ----------
    date_string: str
        formatted time string

    Returns
    -------
    date: obj
        output ``datetime`` object
    """
    # parse the date string
    date = dateutil.parser.parse(date_string)
    # convert to UTC if containing time-zone information
    # then drop the timezone information to prevent unsupported errors
    if date.tzinfo:
        date = date.astimezone(dateutil.tz.UTC).replace(tzinfo=None)
    # return the datetime object
    return date

# PURPOSE: parse a date string into epoch and units scale
def parse_date_string(date_string: str):
    """
    Parse a date string of the form

    - time-units since ``yyyy-mm-dd hh:mm:ss``
    - ``yyyy-mm-dd hh:mm:ss`` for exact calendar dates

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss

    Returns
    -------
    epoch: list
        epoch of ``delta_time``
    conversion_factor: float
        multiplication factor to convert to seconds
    """
    # try parsing the original date string as a date
    try:
        epoch = parse(date_string)
    except ValueError:
        pass
    else:
        # return the epoch (as list)
        return (datetime_to_list(epoch), 0.0)
    # split the date string into units and epoch
    units,epoch = split_date_string(date_string)
    if units not in _to_sec.keys():
        raise ValueError(f'Invalid units: {units}')
    # return the epoch (as list) and the time unit conversion factors
    return (datetime_to_list(epoch), _to_sec[units])

# PURPOSE: split a date string into units and epoch
def split_date_string(date_string: str):
    """
    Split a date string into units and epoch

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss
    """
    try:
        units,_,epoch = date_string.split(None, 2)
    except ValueError:
        raise ValueError(f'Invalid format: {date_string}')
    else:
        return (units.lower(), parse(epoch))

# PURPOSE: convert a datetime object into a list
def datetime_to_list(date):
    """
    Convert a ``datetime`` object into a list

    Parameters
    ----------
    date: obj
        Input ``datetime`` object to convert

    Returns
    -------
    date: list
        [year,month,day,hour,minute,second]
    """
    return [date.year, date.month, date.day,
            date.hour, date.minute, date.second]

# days per month in a leap and a standard year
# only difference is February (29 vs. 28)
_dpm_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_dpm_stnd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# PURPOSE: gets the number of days per month for a given year
def calendar_days(year: int | float | np.ndarray) -> np.ndarray:
    """
    Calculates the number of days per month for a given year

    Parameters
    ----------
    year: int, float or np.ndarray
        calendar year

    Returns
    -------
    dpm: list
        number of days for each month
    """
    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    # find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return np.array(_dpm_leap, dtype=np.float64)
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return np.array(_dpm_stnd, dtype=np.float64)

# PURPOSE: convert a numpy datetime array to delta times since an epoch
def convert_datetime(
        date: float | np.ndarray,
        epoch: str | tuple | list | np.datetime64 = _unix_epoch
    ):
    """
    Convert a ``numpy`` ``datetime`` array to seconds since ``epoch``

    Parameters
    ----------
    date: np.ndarray
        numpy datetime array
    epoch: str, tuple, list, np.ndarray, default (1970,1,1,0,0,0)
        epoch for output ``delta_time``

    Returns
    -------
    delta_time: float
        seconds since epoch
    """
    # convert epoch to datetime variables
    if isinstance(epoch, (tuple, list)):
        epoch = np.datetime64(datetime.datetime(*epoch))
    elif isinstance(epoch, str):
        epoch = np.datetime64(parse(epoch))
    # convert to delta time
    return (date - epoch) / np.timedelta64(1, 's')

# PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(
        delta_time: np.ndarray,
        epoch1: str | tuple | list | np.datetime64 | None = None,
        epoch2: str | tuple | list | np.datetime64 | None = None,
        scale: float = 1.0
    ):
    """
    Convert delta time from seconds since ``epoch1`` to time since ``epoch2``

    Parameters
    ----------
    delta_time: np.ndarray
        seconds since epoch1
    epoch1: tuple or NoneType, default None
        epoch for input ``delta_time``
    epoch2: tuple or NoneType, default None
        epoch for output ``delta_time``
    scale: float, default 1.0
        scaling factor for converting time to output units
    """
    # convert epochs to datetime variables
    if isinstance(epoch1, (tuple, list)):
        epoch1 = np.datetime64(datetime.datetime(*epoch1))
    elif isinstance(epoch1, str):
        epoch1 = np.datetime64(parse(epoch1))
    if isinstance(epoch2, (tuple, list)):
        epoch2 = np.datetime64(datetime.datetime(*epoch2))
    elif isinstance(epoch2, str):
        epoch2 = np.datetime64(parse(epoch2))
    # calculate the total difference in time in seconds
    delta_time_epochs = (epoch2 - epoch1) / np.timedelta64(1, 's')
    # subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

# PURPOSE: calculate the delta time from calendar date
# http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(
        year: np.ndarray,
        month: np.ndarray,
        day: np.ndarray,
        hour: np.ndarray | float = 0.0,
        minute: np.ndarray | float = 0.0,
        second: np.ndarray | float = 0.0,
        epoch: tuple | list | np.datetime64 = _tide_epoch,
        scale: float = 1.0
    ) -> np.ndarray:
    """
    Calculate the time in time units since ``epoch`` from calendar dates

    Parameters
    ----------
    year: np.ndarray
        calendar year
    month: np.ndarray
        month of the year
    day: np.ndarray
        day of the month
    hour: np.ndarray or float, default 0.0
        hour of the day
    minute: np.ndarray or float, default 0.0
        minute of the hour
    second: np.ndarray or float, default 0.0
        second of the minute
    epoch: tuple or list, default icesat2_toolkit.time._tide_epoch
        epoch for output ``delta_time``
    scale: float, default 1.0
        scaling factor for converting time to output units

    Returns
    -------
    delta_time: np.ndarray
        time since epoch
    """
    # calculate date in Modified Julian Days (MJD) from calendar date
    # MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    # convert epochs to datetime variables
    epoch1 = np.datetime64(datetime.datetime(*_mjd_epoch))
    if isinstance(epoch, (tuple, list)):
        epoch = np.datetime64(datetime.datetime(*epoch))
    elif isinstance(epoch, str):
        epoch = np.datetime64(parse(epoch))
    # calculate the total difference in time in days
    delta_time_epochs = (epoch - epoch1) / np.timedelta64(1, 'D')
    # return the date in units (default days) since epoch
    return scale*np.array(MJD - delta_time_epochs, dtype=np.float64)

# PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(
        year: np.ndarray,
        month: np.ndarray,
        day: np.ndarray,
        hour: np.ndarray | float | None = None,
        minute: np.ndarray | float | None = None,
        second: np.ndarray | float | None = None,
        DofY: np.ndarray | float | None = None,
    ) -> np.ndarray:
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Parameters
    ----------
    year: np.ndarray
        calendar year
    month: np.ndarray
        calendar month
    day: np.ndarray, float or NoneType, default None
        day of the month
    hour: np.ndarray, float or NoneType, default None
        hour of the day
    minute: np.ndarray, float or NoneType, default None
        minute of the hour
    second: np.ndarray, float or NoneType, default None
        second of the minute
    DofY: np.ndarray, float or NoneType, default None
        day of the year

    Returns
    -------
    t_date: np.ndarray
        date in decimal-year format

    References
    ----------
    .. [1] N. Dershowitz, and E. M. Reingold.
        *Calendrical Calculations*,
        Cambridge: Cambridge University Press, (2008).
    """

    # number of dates
    n_dates = len(np.atleast_1d(year))

    # create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    # day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    # remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    # create output date variable
    t_date = np.zeros((n_dates))

    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array(_dpm_leap, dtype=np.float64)
    dpm_stnd = np.array(_dpm_stnd, dtype=np.float64)

    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    # find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    # calculate the day of the year
    if DofY is not None:
        # if entered directly as an input
        # remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        # use calendar month and day of the month to calculate day of the year
        # month minus 1: January = 0, February = 1, etc (indice of month)
        # in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=int) - 1

        # day of month
        if day is not None:
            # remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            # if not entering days as an input
            # will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        # create matrix with the lower half = 1
        # this matrix will be used in a matrix multiplication
        # to calculate the total number of days for prior months
        # the -1 will make the diagonal == 0
        # i.e. first row == all zeros and the
        # last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        # using a dot product to calculate total number of days
        # for the months before the input date
        # basically is sum(i*dpm)
        # where i is 1 for all months < the month of interest
        # and i is 0 for all months >= the month of interest
        # month of interest is zero as the exact days will be
        # used to calculate the date

        # calculate the day of the year for leap and standard
        # use total days of all months before date
        # and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    # hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    # minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    # second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    # calculate decimal date
    # convert hours, minutes and seconds into days
    # convert calculated fractional days into decimal fractions of the year
    # Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    # Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

# PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD: np.ndarray, **kwargs):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: np.ndarray
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    astype: str or NoneType, default None
        convert output to variable type
    format: str, default 'dict'
        format of output variables

            - ``'dict'``: dictionary with variable keys
            - ``'tuple'``: tuple in most-to-least-significant order
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    year: np.ndarray
        calendar year
    month: np.ndarray
        calendar month
    day: np.ndarray
        day of the month
    hour: np.ndarray
        hour of the day
    minute: np.ndarray
        minute of the hour
    second: np.ndarray
        second of the minute

    References
    ----------
    .. [1] W. H. Press, *Numerical Recipes in C*,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, (1988).
    .. [2] D. A. Hatcher, "Simple Formulae for Julian Day Numbers and
        Calendar Dates", *Quarterly Journal of the Royal Astronomical
        Society*, 25(1), 1984.
    """
    # set default keyword arguments
    kwargs.setdefault('astype', None)
    kwargs.setdefault('format', 'dict')
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(ASTYPE='astype', FORMAT='format')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        single_value = True
    else:
        single_value = False

    # verify julian day
    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    # calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    # calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    # calculate day, month, year and hour
    day = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    month = F - 1.0 - 12.0*np.floor(F/14.0)
    year = D - 4715.0 - np.floor((7.0 + month)/10.0)
    hour = np.floor(24.0*(JD + 0.5 - JDO))
    # calculate minute and second
    G = (JD + 0.5 - JDO) - hour/24.0
    minute = np.floor(G*1440.0)
    second = (G - minute/1440.0) * 86400.0

    # convert all variables to output type (from float)
    if kwargs['astype'] is not None:
        year = year.astype(kwargs['astype'])
        month = month.astype(kwargs['astype'])
        day = day.astype(kwargs['astype'])
        hour = hour.astype(kwargs['astype'])
        minute = minute.astype(kwargs['astype'])
        second = second.astype(kwargs['astype'])

    # if only a single value was imported initially: remove singleton dims
    if single_value:
        year = year.item(0)
        month = month.item(0)
        day = day.item(0)
        hour = hour.item(0)
        minute = minute.item(0)
        second = second.item(0)

    # return date variables in output format
    if (kwargs['format'] == 'dict'):
        return dict(year=year, month=month, day=day,
            hour=hour, minute=minute, second=second)
    elif (kwargs['format'] == 'tuple'):
        return (year, month, day, hour, minute, second)
    elif (kwargs['format'] == 'zip'):
        return zip(year, month, day, hour, minute, second)

# PURPOSE: Count number of leap seconds that have passed for each GPS time
def count_leap_seconds(
        GPS_Time: np.ndarray | float,
        truncate: bool = True
    ):
    """
    Counts the number of leap seconds between a given GPS time and UTC

    Parameters
    ----------
    GPS_Time: np.ndarray or float
        seconds since January 6, 1980 at 00:00:00
    truncate: bool, default True
        Reduce list of leap seconds to positive GPS times

    Returns
    -------
    n_leaps: float
        number of elapsed leap seconds
    """
    # get the valid leap seconds
    leaps = get_leap_seconds(truncate=truncate)
    # number of leap seconds prior to GPS_Time
    n_leaps = np.zeros_like(GPS_Time,dtype=np.float64)
    for i,leap in enumerate(leaps):
        count = np.count_nonzero(GPS_Time >= leap)
        if (count > 0):
            indices = np.nonzero(GPS_Time >= leap)
            n_leaps[indices] += 1.0
    # return the number of leap seconds for converting to UTC
    return n_leaps

# PURPOSE: Define GPS leap seconds
def get_leap_seconds(truncate: bool = True):
    """
    Gets a list of GPS times for when leap seconds occurred

    Parameters
    ----------
    truncate: bool, default True
        Reduce list of leap seconds to positive GPS times

    Returns
    -------
    GPS time: float
        GPS seconds when leap seconds occurred
    """
    leap_secs = icesat2_toolkit.utilities.get_data_path(['data','leap-seconds.list'])
    # find line with file expiration as delta time
    with leap_secs.open(mode='r', encoding='utf8') as fid:
        secs, = [re.findall(r'\d+',i).pop() for i in fid.read().splitlines()
            if re.match(r'^(?=#@)',i)]
    # check that leap seconds file is still valid
    expiry = datetime.datetime(*_ntp_epoch) + datetime.timedelta(seconds=int(secs))
    today = datetime.datetime.utcnow()
    update_leap_seconds() if (expiry < today) else None
    # get leap seconds
    leap_UTC,TAI_UTC = np.loadtxt(leap_secs).T
    # TAI time is ahead of GPS by 19 seconds
    TAI_GPS = 19.0
    # convert leap second epochs from NTP to GPS
    # convert from time of 2nd leap second to time of 1st leap second
    leap_GPS = convert_delta_time(leap_UTC + TAI_UTC - TAI_GPS - 1,
        epoch1=_ntp_epoch, epoch2=_gps_epoch)
    # return the GPS times of leap second occurance
    if truncate:
        return leap_GPS[leap_GPS >= 0].astype(np.float64)
    else:
        return leap_GPS.astype(np.float64)

# PURPOSE: connects to servers and downloads leap second files
def update_leap_seconds(
        timeout: int | None = 20,
        verbose: bool = False,
        mode: oct = 0o775
    ):
    """
    Connects to servers to download leap-seconds.list files from NIST servers

    - https://www.nist.gov/pml/time-and-frequency-division/leap-seconds-faqs

    Servers and Mirrors

    - ftp://ftp.nist.gov/pub/time/leap-seconds.list
    - https://www.ietf.org/timezones/data/leap-seconds.list

    Parameters
    ----------
    timeout: int or None, default 20
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    # local version of file
    FILE = 'leap-seconds.list'
    LOCAL = icesat2_toolkit.utilities.get_data_path(['data',FILE])
    HASH = icesat2_toolkit.utilities.get_hash(LOCAL)

    # try downloading from NIST ftp servers
    HOST = ['ftp.nist.gov','pub','time',FILE]
    try:
        icesat2_toolkit.utilities.check_ftp_connection(HOST[0])
        icesat2_toolkit.utilities.from_ftp(HOST, timeout=timeout, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as exc:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

    # try downloading from Internet Engineering Task Force (IETF) mirror
    REMOTE = ['https://www.ietf.org','timezones','data',FILE]
    try:
        icesat2_toolkit.utilities.from_http(REMOTE, timeout=timeout, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as exc:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

    # return exception that no server could be connected
    raise Exception('All Server Connection Error')

