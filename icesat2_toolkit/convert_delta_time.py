#!/usr/bin/env python
u"""
convert_delta_time.py (08/2020)
Converts time from delta seconds into Julian and year-decimal

INPUTS:
    delta_time: seconds since gps_epoch

OPTIONS:
    gps_epoch: seconds between delta_time and GPS epoch (1980-01-06T00:00:00)
        ICESat-2 atlas_sdp_gps_epoch: 1198800018.0

OUTPUTS:
    julian: time in Julian days
    decimal: time in year-decimal

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

PROGRAM DEPENDENCIES:
    convert_julian.py: converts from Julian date to calendar date and time
    convert_calendar_decimal.py: converts from calendar date to decimal year
    time.py: Utilities for calculating time operations

UPDATE HISTORY:
    Updated 08/2020: time utilities for counting leap seconds and JD conversion
    Updated 07/2020: added function docstrings
    Written 05/2020
"""
import numpy as np
import icesat2_toolkit.time
from icesat2_toolkit.convert_julian import convert_julian
from icesat2_toolkit.convert_calendar_decimal import convert_calendar_decimal

#-- PURPOSE: convert time from delta seconds into Julian and year-decimal
def convert_delta_time(delta_time, gps_epoch=1198800018.0):
    """
    converts ICESat-2 delta_times into into Julian and year-decimal

    Arguments
    ---------
    delta_time: seconds since gps_epoch

    Keyword arguments
    -----------------
    gps_epoch: seconds between delta_time and GPS epoch (1980-01-06T00:00:00)

    Returns
    -------
    julian: time in Julian days
    decimal: time in year-decimal
    """
    #-- convert to array if single value
    if (np.ndim(delta_time) == 0):
        delta_time = np.array([delta_time])
    #-- calculate gps time from delta_time
    gps_seconds = gps_epoch + delta_time
    time_leaps = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
    #-- calculate Julian time (UTC) by converting to MJD and then adding offset
    time_julian = 2400000.5 + icesat2_toolkit.time.convert_delta_time(
        gps_seconds - time_leaps, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date with convert_julian.py
    Y,M,D,h,m,s = convert_julian(time_julian,FORMAT='tuple')
    #-- calculate year-decimal time (UTC)
    time_decimal = convert_calendar_decimal(Y,M,DAY=D,HOUR=h,MINUTE=m,SECOND=s)
    #-- return both the Julian and year-decimal formatted dates
    return dict(julian=np.squeeze(time_julian),decimal=np.squeeze(time_decimal))
