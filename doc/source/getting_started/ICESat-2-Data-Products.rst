ICESat-2 Data Products
======================

.. _table1:
.. include:: Table.rst

ICESat-2 Granules
#################

Each orbit of ICESat-2 data is broken up into 14 granule regions.
The granule boundaries limit the size of each ATL03 file and simplify the formation of higher level data products.

.. image:: https://raw.githubusercontent.com/tsutterley/read-ICESat-2/main/icesat2_toolkit/data/ICESat-2_granules_global.png
  :alt: ICESat-2 global granules

.. image:: https://raw.githubusercontent.com/tsutterley/read-ICESat-2/main/icesat2_toolkit/data/ICESat-2_granules_polar.png
  :alt: ICESat-2 polar granules

ATLAS Beam Strength
###################

The Advanced Topographic Laser Altimeter System (ATLAS) is a photon-counting laser altimeter and the primary instrumentation onboard the ICESat-2 observatory.
ATLAS sends and receives data for 6 individual beams that are separated into three beam pairs.
The three beam pairs are separated by 3 kilometers, and the two paired beams are separated on the ground by 90 meters.
Each beam pair consists of a weak beam and a strong beam, with the strong beam approximately four times brighter than weak.
The spots illuminated on the ground from the ATLAS strong beams are numbered 1, 3, and 5, and the spots illuminated on the ground from the ATLAS weak beams are numbered 2, 4, and 6.

The ICESat-2 observatory can be oriented in one of two positions with respect to the direction of travel.
In the forward orientation, the weak beams lead the strong beams and a weak beam is on the left edge of the beam pattern (gt1l).
This is reversed in the backward orientation, and the strong beams lead the weak beams with a strong beam on the left edge of the beam pattern (:numref:`table2`).
The beam strength, spot number, atmospheric profile and spacecraft orientation can all be found as HDF5 attributes of each beam group.

.. _table2:
.. table:: ATLAS beam mapping when in the forward and backward orientations
    :align: center

    +----------------+------------------------+------------------------+
    | | Ground Track | | Forward Orientation  | | Backward Orientation |
    | | Designation  | | Beam Strength        | | Beam Strength        |
    +================+========================+========================+
    |      gt1l      |          Weak          |         Strong         |
    +----------------+------------------------+------------------------+
    |      gt1r      |         Strong         |          Weak          |
    +----------------+------------------------+------------------------+
    |      gt2l      |          Weak          |         Strong         |
    +----------------+------------------------+------------------------+
    |      gt2r      |         Strong         |          Weak          |
    +----------------+------------------------+------------------------+
    |      gt3l      |          Weak          |         Strong         |
    +----------------+------------------------+------------------------+
    |      gt3r      |         Strong         |          Weak          |
    +----------------+------------------------+------------------------+

