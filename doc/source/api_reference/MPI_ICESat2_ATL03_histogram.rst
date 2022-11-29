==============================
MPI_ICESat2_ATL03_histogram.py
==============================

- Read ICESat-2 ATL03 and ATL09 data files to calculate average segment surfaces

    * ATL03 datasets: Global Geolocated Photons
    * ATL09 datasets: Atmospheric Characteristics
- Alternative algorithm that uses gaussian/generalized gaussian decomposition to extract possibly multiple height surfaces from a histogram of photon events

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_ICESat2_ATL03_histogram.py

Calling Sequence
################

.. argparse::
    :filename: MPI_ICESat2_ATL03_histogram.py
    :func: arguments
    :prog: MPI_ICESat2_ATL03_histogram.py
    :nodescription:
    :nodefault: