#!/usr/bin/env python
u"""
tools.py
Written by Tyler Sutterley (04/2022)
Plotting tools and utilities

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2021: added custom colormap function for some common scales
    Written 09/2021
"""
import os
import re
import copy
import colorsys
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

def from_cpt(filename, use_extremes=True, **kwargs):
    """
    Reads GMT color palette table files and registers the
    colormap to be recognizable by ``plt.cm.get_cmap()``

    Can import HSV (hue-saturation-value) or RGB values

    Parameters
    ----------
    filename: str
        color palette table file
    use_extremes: bool, default True
        use the under, over and bad values from the cpt file
    **kwargs: dict
        optional arguments for LinearSegmentedColormap
    """

    # read the cpt file and get contents
    with open(filename,'r') as f:
        file_contents = f.read().splitlines()
    # extract basename from cpt filename
    name = re.sub(r'\.cpt','',os.path.basename(filename),flags=re.I)

    # compile regular expression operator to find numerical instances
    rx = re.compile(r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?')

    # create list objects for x, r, g, b
    x,r,g,b = ([],[],[],[])
    # assume RGB color model
    colorModel = "RGB"
    # back, forward and no data flags
    flags = dict(B=None,F=None,N=None)
    for line in file_contents:
        # find back, forward and no-data flags
        model = re.search(r'COLOR_MODEL.*(HSV|RGB)',line,re.I)
        BFN = re.match(r'[BFN]',line,re.I)
        # parse non-color data lines
        if model:
            # find color model
            colorModel = model.group(1)
            continue
        elif BFN:
            flags[BFN.group(0)] = [float(i) for i in rx.findall(line)]
            continue
        elif re.search(r"#",line):
            # skip over commented header text
            continue
        # find numerical instances within line
        x1,r1,g1,b1,x2,r2,g2,b2 = rx.findall(line)
        # append colors and locations to lists
        x.append(float(x1))
        r.append(float(r1))
        g.append(float(g1))
        b.append(float(b1))
    # append end colors and locations to lists
    x.append(float(x2))
    r.append(float(r2))
    g.append(float(g2))
    b.append(float(b2))

    # convert input colormap to output
    xNorm = [None]*len(x)
    if (colorModel == "HSV"):
        # convert HSV (hue-saturation-value) to RGB
        # calculate normalized locations (0:1)
        for i,xi in enumerate(x):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
            xNorm[i] = (xi - x[0])/(x[-1] - x[0])
    elif (colorModel == "RGB"):
        # normalize hexadecimal RGB triple from (0:255) to (0:1)
        # calculate normalized locations (0:1)
        for i,xi in enumerate(x):
            r[i] /= 255.0
            g[i] /= 255.0
            b[i] /= 255.0
            xNorm[i] = (xi - x[0])/(x[-1] - x[0])

    # output RGB lists containing normalized location and colors
    cdict = dict(red=[None]*len(x),green=[None]*len(x),blue=[None]*len(x))
    for i,xi in enumerate(x):
        cdict['red'][i] = [xNorm[i],r[i],r[i]]
        cdict['green'][i] = [xNorm[i],g[i],g[i]]
        cdict['blue'][i] = [xNorm[i],b[i],b[i]]

    # create colormap for use in matplotlib
    cmap = colors.LinearSegmentedColormap(name, cdict, **kwargs)
    # set flags for under, over and bad values
    extremes = dict(under=None,over=None,bad=None)
    for key,attr in zip(['B','F','N'],['under','over','bad']):
        if flags[key] is not None:
            r,g,b = flags[key]
            if (colorModel == "HSV"):
                # convert HSV (hue-saturation-value) to RGB
                r,g,b = colorsys.hsv_to_rgb(r/360.,g,b)
            elif (colorModel == 'RGB'):
                # normalize hexadecimal RGB triple from (0:255) to (0:1)
                r,g,b = (r/255.0,g/255.0,b/255.0)
            # set attribute for under, over and bad values
            extremes[attr] = (r,g,b)
    # create copy of colormap with extremes
    if use_extremes:
        cmap = cmap.with_extremes(**extremes)
    # register colormap to be recognizable by cm.get_cmap()
    cm.register_cmap(name=name, cmap=cmap)
    # return the colormap
    return cmap

def custom_colormap(N, map_name, **kwargs):
    """
    Calculates a custom colormap and registers it
    to be recognizable by ``plt.cm.get_cmap()``

    Parameters
    ----------
    N: int
        number of slices in initial HSV color map
    map_name: str
        name of color map

            - ``'Joughin'``: [Joughin2018]_ standard velocity colormap
            - ``'Rignot'``: [Rignot2011]_ standard velocity colormap
            - ``'Seroussi'``: [Seroussi2011]_ velocity divergence colormap
    **kwargs: dict
        optional arguments for LinearSegmentedColormap

    References
    ----------
    .. [Joughin2018] I. Joughin, B. E. Smith, and I. Howat,
        "Greenland Ice Mapping Project: ice flow velocity variation
        at sub-monthly to decadal timescales", *The Cryosphere*,
        12, 2211--2227, (2018).
        `doi: 10.5194/tc-12-2211-2018 <https://doi.org/10.5194/tc-12-2211-2018>`_
    .. [Rignot2011] E. Rignot J. Mouginot, and B. Scheuchl,
        "Ice Flow of the Antarctic Ice Sheet",
        *Science*, 333(6048), 1427--1430, (2011).
        `doi: 10.1126/science.1208336 <https://doi.org/10.1126/science.1208336>`_
    .. [Seroussi2011] H. Seroussi, M. Morlighem, E. Rignot, E. Larour,
        D. Aubry, H. Ben Dhia, and S. S. Kristensen,
        "Ice flux divergence anomalies on 79north Glacier, Greenland",
        *Geophysical Research Letters*, 38(L09501), (2011).
        `doi: 10.1029/2011GL047338 <https://doi.org/10.1029/2011GL047338>`_
    """

    # make sure map_name is properly formatted
    map_name = map_name.capitalize()
    if (map_name == 'Joughin'):
        # calculate initial HSV for Ian Joughin's color map
        h = np.linspace(0.1,1,N)
        s = np.ones((N))
        v = np.ones((N))
        # calculate RGB color map from HSV
        color_map = np.zeros((N,3))
        for i in range(N):
            color_map[i,:] = colorsys.hsv_to_rgb(h[i],s[i],v[i])
    elif (map_name == 'Seroussi'):
        # calculate initial HSV for Helene Seroussi's color map
        h = np.linspace(0,1,N)
        s = np.ones((N))
        v = np.ones((N))
        # calculate RGB color map from HSV
        RGB = np.zeros((N,3))
        for i in range(N):
            RGB[i,:] = colorsys.hsv_to_rgb(h[i],s[i],v[i])
        # reverse color order and trim to range
        RGB = RGB[::-1,:]
        RGB = RGB[1:np.floor(0.7*N).astype('i'),:]
        # calculate HSV color map from RGB
        HSV = np.zeros_like(RGB)
        for i,val in enumerate(RGB):
            HSV[i,:] = colorsys.rgb_to_hsv(val[0],val[1],val[2])
        # calculate saturation as a function of hue
        HSV[:,1] = np.clip(0.1 + HSV[:,0], 0, 1)
        # calculate RGB color map from HSV
        color_map = np.zeros_like(HSV)
        for i,val in enumerate(HSV):
            color_map[i,:] = colorsys.hsv_to_rgb(val[0],val[1],val[2])
    elif (map_name == 'Rignot'):
        # calculate initial HSV for Eric Rignot's color map
        h = np.linspace(0,1,N)
        s = np.clip(0.1 + h, 0, 1)
        v = np.ones((N))
        # calculate RGB color map from HSV
        color_map = np.zeros((N,3))
        for i in range(N):
            color_map[i,:] = colorsys.hsv_to_rgb(h[i],s[i],v[i])
    else:
        raise ValueError('Incorrect color map specified ({0})'.format(map_name))

    # output RGB lists containing normalized location and colors
    Xnorm = len(color_map) - 1.0
    cdict = dict(red=[None]*len(color_map),
        green=[None]*len(color_map),
        blue=[None]*len(color_map))
    for i,rgb in enumerate(color_map):
        cdict['red'][i] = [float(i)/Xnorm,rgb[0],rgb[0]]
        cdict['green'][i] = [float(i)/Xnorm,rgb[1],rgb[1]]
        cdict['blue'][i] = [float(i)/Xnorm,rgb[2],rgb[2]]

    # create colormap for use in matplotlib
    cmap = colors.LinearSegmentedColormap(map_name, cdict, **kwargs)
    # register colormap to be recognizable by cm.get_cmap()
    cm.register_cmap(name=map_name, cmap=cmap)
    # return the colormap
    return cmap
