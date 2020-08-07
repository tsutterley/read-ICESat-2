#!/usr/bin/env python
u"""
count_leap_seconds.py (08/2020)
Count number of leap seconds that have passed for each GPS time
"""
import warnings
import icesat2_toolkit.time

def count_leap_seconds(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use icesat2_toolkit.time instead",
        DeprecationWarning)
    # call renamed version to not break workflows
    return icesat2_toolkit.time.count_leap_seconds(*args)
