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

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/time.py


General Methods
===============

.. autofunction:: icesat2_toolkit.time.parse_date_string

.. autofunction:: icesat2_toolkit.time.split_date_string

.. autofunction:: icesat2_toolkit.time.datetime_to_list

.. autofunction:: icesat2_toolkit.time.calendar_days

.. autofunction:: icesat2_toolkit.time.convert_datetime

.. autofunction:: icesat2_toolkit.time.convert_delta_time

.. autofunction:: icesat2_toolkit.time.convert_calendar_dates

.. autofunction:: icesat2_toolkit.time.convert_calendar_decimal

.. autofunction:: icesat2_toolkit.time.convert_julian

.. autofunction:: icesat2_toolkit.time.count_leap_seconds

.. autofunction:: icesat2_toolkit.time.get_leap_seconds

.. autofunction:: icesat2_toolkit.time.update_leap_seconds
