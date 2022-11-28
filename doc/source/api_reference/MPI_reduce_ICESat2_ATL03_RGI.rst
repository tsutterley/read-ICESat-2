===============================
MPI_reduce_ICESat2_ATL03_RGI.py
===============================

- Create masks for reducing ICESat-2 ATL03 data to the `Randolph Glacier Inventory <https://www.glims.org/RGI/rgi60_dl.html>`_

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_reduce_ICESat2_ATL03_RGI.py

Calling Sequence
################

.. argparse::
    :filename: MPI_reduce_ICESat2_ATL03_RGI.py
    :func: arguments
    :prog: MPI_reduce_ICESat2_ATL03_RGI.py
    :nodescription:
    :nodefault:

    --region -r : @replace
        Possible choices: 1--19

        Region of Randolph Glacier Inventory to run

        1. Alaska
        2. Western Canada and USA
        3. Arctic Canada North
        4. Arctic Canada South
        5. Greenland Periphery
        6. Iceland
        7. Svalbard
        8. Scandinavia
        9. Russian Arctic
        10. North Asia
        11. Central Europe
        12. Caucasus, Middle East
        13. Central Asia
        14. South Asia West
        15. South Asia East
        16. Low Latitudes
        17. Southern Andes
        18. New Zealand
        19. Antarctic, Subantarctic
