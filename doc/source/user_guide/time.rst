=======
time.py
=======

Utilities for calculating time operations

 - Can convert delta time from seconds since an epoch to time since a different epoch
 - Can calculate the time in days since epoch from calendar dates
 - Can count the number of leap seconds between a given GPS time and UTC
 - Syncs leap second files with NIST servers

Calling Sequence
================

Count the number of leap seconds between a GPS time and UTC

.. code-block:: python

    import icesat2_toolkit.time
    leap_seconds = icesat2_toolkit.time.count_leap_seconds(gps_seconds)

Convert a time from seconds since 1980-01-06T00:00:00 to Modified Julian Days (MJD)

.. code-block:: python

    import icesat2_toolkit.time
    MJD = icesat2_toolkit.time.convert_delta_time(delta_time, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)

Convert a calendar date into Modified Julian Days

.. code-block:: python

    import icesat2_toolkit.time
    MJD = icesat2_toolkit.time.convert_calendar_dates(YEAR,MONTH,DAY,hour=HOUR,
        minute=MINUTE,second=SECOND,epoch=(1858,11,17,0,0,0))

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/master/icesat2_toolkit/time.py


General Methods
===============


.. method:: icesat2_toolkit.time.convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0)

    Convert delta time from seconds since epoch1 to time since epoch2

    Arguments:

        `delta_time`: seconds since epoch1

    Keyword arguments:

        `epoch1`: epoch for input delta_time

        `epoch2`: epoch for output delta_time

        `scale`: scaling factor for converting time to output units


.. method:: icesat2_toolkit.time.convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0, epoch=None)

    Calculate the time in days since epoch from calendar dates

    Arguments:

        `year`: calendar month

        `month`: month of the year

        `day`: day of the month

    Keyword arguments:

        `hour`: hour of the day

        `minute`: minute of the hour

        `second`: second of the minute

        `epoch`: epoch for output delta_time


.. method:: icesat2_toolkit.time.count_leap_seconds(GPS_Time)

    Counts the number of leap seconds between a given GPS time and UTC

    Arguments:

        `GPS_Time`: seconds since January 6, 1980 at 00:00:00


.. method:: icesat2_toolkit.time.get_leap_seconds()

    Gets a list of GPS times for when leap seconds occurred


.. method:: icesat2_toolkit.time.update_leap_seconds(verbose=False, mode=0o775)

    Connects to servers to download leap-seconds.list files from `NIST servers`__

.. __: ftp://ftp.nist.gov/pub/time/leap-seconds.list

    Keyword arguments:

        `verbose`: print file information about output file

        `mode`: permissions mode of output file
