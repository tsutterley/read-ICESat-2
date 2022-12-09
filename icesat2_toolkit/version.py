#!/usr/bin/env python
u"""
version.py (04/2021)
Gets version number of a package
"""
from pkg_resources import get_distribution

# get version
version = get_distribution("icesat2_toolkit").version
# append "v" before the version
full_version = f"v{version}"
# get project name
project_name = get_distribution('icesat2_toolkit').project_name
