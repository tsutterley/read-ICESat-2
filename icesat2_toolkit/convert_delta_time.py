#!/usr/bin/env python
u"""
convert_delta_time.py (04/2022)
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
    time.py: Utilities for calculating time operations

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 08/2020: time utilities for counting leap seconds and JD conversion
    Updated 07/2020: added function docstrings
    Written 05/2020
"""
import numpy as np
import icesat2_toolkit.time

#-- PURPOSE: convert time from delta seconds into Julian and year-decimal
def convert_delta_time(delta_time, gps_epoch=1198800018.0):
    """
    converts ICESat-2 delta_times into into Julian and year-decimal

    Parameters
    ----------
    delta_time: float
        seconds since gps_epoch
    gps_epoch: float, default 1198800018.0
        seconds between delta_time and GPS epoch (1980-01-06T00:00:00)

    Returns
    -------
    julian: float
        time in Julian days
    decimal: float
        time in year-decimal
    """
    #-- convert to array if single value
    delta_time = np.atleast_1d(delta_time)
    #-- calculate gps time from delta_time
    gps_seconds = gps_epoch + delta_time
    time_leaps = icesat2_toolkit.time.count_leap_seconds(gps_seconds)
    #-- calculate Julian time (UTC) by converting to MJD and then adding offset
    time_julian = 2400000.5 + icesat2_toolkit.time.convert_delta_time(
        gps_seconds - time_leaps, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date
    Y,M,D,h,m,s = icesat2_toolkit.time.convert_julian(time_julian,format='tuple')
    #-- calculate year-decimal time (UTC)
    time_decimal = icesat2_toolkit.time.convert_calendar_decimal(Y,M,day=D,
        hour=h,minute=m,second=s)
    #-- return both the Julian and year-decimal formatted dates
    return dict(julian=np.squeeze(time_julian),decimal=np.squeeze(time_decimal))
