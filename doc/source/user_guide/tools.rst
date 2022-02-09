========
tools.py
========

Plotting tools and utilities

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/tools.py

General Attributes and Methods
==============================

.. function:: from_cpt(filename, use_extremes=True)

    Reads GMT color palette table files and registers the colormap to be recognizable by ``plt.cm.get_cmap()``

.. function:: custom_colormap(N, map_name)

    Calculates a custom colormap and registers it to be recognizable by ``plt.cm.get_cmap()``

    Arguments:

        - ``N``: number of slices in initial HSV color map

        - ``map_name``: name of color map

            * ``'Joughin'``: [Joughin2018]_ standard velocity colormap

            * ``'Rignot'``: [Rignot2011]_ standard velocity colormap

            * ``'Seroussi'``:  [Seroussi2011]_ velocity divergence colormap


References
##########

.. [Joughin2018] I. Joughin, B. E. Smith, and I. Howat, "Greenland Ice Mapping Project: ice flow velocity variation at sub-monthly to decadal timescales", *The Cryosphere*, 12, 2211--2227, (2018). `doi: 10.5194/tc-12-2211-2018 <https://doi.org/10.5194/tc-12-2211-2018>`_

.. [Rignot2011] E. Rignot J. Mouginot, and B. Scheuchl, "Ice Flow of the Antarctic Ice Sheet", *Science*, 333(6048), 1427--1430, (2011). `doi: 10.1126/science.1208336 <https://doi.org/10.1126/science.1208336>`_

.. [Seroussi2011] H. Seroussi, M. Morlighem, E. Rignot, E. Larour, D. Aubry, H. Ben Dhia, and S. S. Kristensen, "Ice flux divergence anomalies on 79north Glacier, Greenland", *Geophysical Research Letters*, 38(L09501), (2011). `doi: 10.1029/2011GL047338 <https://doi.org/10.1029/2011GL047338>`_
