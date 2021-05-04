======
fit.py
======

Utilities for calculating average fits from ATL03 Geolocated Photon Data

Calling Sequence
================

Count the number of leap seconds between a GPS time and UTC

.. code-block:: python

    import icesat2_toolkit.fit
    pe_weights = icesat2_toolkit.fit.classify_photons(x, h, h_win_width,
        indices, K=5, MIN_PH=5, MIN_XSPREAD=1.0, MIN_HSPREAD=0.01)

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/fit.py


General Methods
===============

.. method:: icesat2_toolkit.fit.compress_list(i,n)

    Compress complete list of valid indices into a set of ranges

    Arguments:

        ``i``: indices to compress

        ``n``: largest gap between indices to accept for range


.. method:: icesat2_toolkit.fit.windowed_manhattan(u, v, window=[], w=None)

    Create a windowed manhattan distance metric

    Arguments:

        ``u``: Input array

        ``v``: Input array for distance

    Keyword arguments:

        ``window``: distance window for reducing neighbors

        ``w``: weights for each value


.. method:: icesat2_toolkit.fit.classify_photons(x, h, h_win_width, indices, K=5, MIN_PH=5, MIN_XSPREAD=1.0, MIN_HSPREAD=0.01, METHOD='ball_tree')

    Use the NASA GSFC YAPC k-nearest neighbors algorithm to determine weights for each photon event within an ATL03 major frame

    Arguments:

        ``x``: along-track x coordinates for photon events for 3 major frames

        ``h``: photon event heights for 3 major frames

        ``h_win_width``: height of (possibly 2) telemetry bands

        ``indices``: indices of photon events in ATL03 major frame

    Keyword arguments:

        ``K``: number of values for KNN algorithm

        ``MIN_PH``: minimum number of photons for a major frame to be valid

        ``MIN_XSPREAD``: minimum along-track spread of photon events

        ``MIN_HSPREAD``: minimum window of heights for photon events

        ``METHOD``: algorithm for computing photon event weights

            ``'ball_tree'``: use scikit.learn.BallTree with custom distance metric

            ``'brute'``: use a brute-force approach


.. method:: icesat2_toolkit.fit.extract_tep_histogram(tep_hist_time,tep_hist,tep_range_prim)

    Centers the transmit-echo-path histogram reported by ATL03 using an iterative edit to distinguish between signal and noise

    Arguments:

        ``tep_hist_time``

        ``tep_hist``

        ``tep_range_prim``


.. method:: icesat2_toolkit.fit.filter_elevation(r0)

    Calculates the interquartile range [Pritchard2009]_ and robust dispersion estimator [Smith2017]_ of the model residuals

    Arguments:

        ``r0``: height residuals


.. method:: icesat2_toolkit.fit.try_surface_fit(x, y, z, confidence_mask, dist_along, SURF_TYPE='linear', ITERATE=25, CONFIDENCE=[4,3,2,1,0])

    Try fitting a surface to the signal photons with progressively less confidence if no valid surface is found

    Arguments:

        ``x``: along-track x-coordinates

        ``y``: along-track y-coordinates

        ``z``: along-track photon heights

        ``confidence_mask``: confidence level of each photon event

        ``dist_along``: center of segment in along-track x-coordinates

    Keyword arguments:

        ``SURF_TYPE``: surface polynomial to fit to photon heights

            ``'linear'``

            ``'quadratic'``

        ``ITERATE``: maximum number of iterations to use in fit

        ``CONFIDENCE``: minimum photon confidence levels to use in fit


.. method:: icesat2_toolkit.fit.reduce_surface_fit(x, y, z, centroid, ind, SURF_TYPE='linear', ITERATE=25)

    Iteratively fit a polynomial surface to the elevation data to reduce to within a valid surface window [Smith2019]_

    Arguments:

        ``x``: along-track x-coordinates

        ``y``: along-track y-coordinates

        ``z``: along-track photon heights

        ``centroid``: segment center for referencing along-track coordinates

        ``ind``: indices of photon events for confidence level to use in fit

    Keyword arguments:

        ``SURF_TYPE``: surface polynomial to fit to photon heights

            ``'linear'``

            ``'quadratic'``

        ``ITERATE``: maximum number of iterations to use in fit

.. method:: icesat2_toolkit.fit.fit_surface(x, y, z, centroid, SURF_TYPE='linear')

    Fit a polynomial surface to the elevation data

    Arguments:

        ``x``: along-track x-coordinates

        ``y``: along-track y-coordinates

        ``z``: along-track photon heights

        ``centroid``: segment center for referencing along-track coordinates

    Keyword arguments:

        ``SURF_TYPE``: surface polynomial to fit to photon heights

            ``'linear'``

            ``'quadratic'``

.. method:: icesat2_toolkit.fit.try_histogram_fit(x, y, z, confidence_mask, dist_along, dt, FIT_TYPE='gaussian', ITERATE=25, BACKGROUND=0, CONFIDENCE=[2,1,0])

    Try fitting a function to the signal photon histograms with progressively less confidence if no valid fit is found

    Arguments:

        ``x``: along-track x-coordinates

        ``y``: along-track y-coordinates

        ``z``: along-track photon heights

        ``confidence_mask``: confidence level of each photon event

        ``dist_along``: center of segment in along-track x-coordinates

        ``dt``: histogram bin size in seconds

    Keyword arguments:

        ``FIT_TYPE``: decomposition function to fit to photon height histograms

            ``'gaussian'``

            ``'general'``

        ``ITERATE``: maximum number of iterations to use in fit

        ``BACKGROUND``: vertical noise-photon density for segment

        ``CONFIDENCE``: minimum photon confidence levels to use in fit


.. method:: icesat2_toolkit.fit.reduce_histogram_fit(x, y, z, ind, dt, FIT_TYPE='gaussian', ITERATE=25, PEAKS=2, BACKGROUND=0)

    Iteratively use decomposition fitting to the elevation data to reduce to within a valid surface window

    Arguments:

        ``x``: along-track x-coordinates

        ``y``: along-track y-coordinates

        ``z``: along-track photon heights

        ``confidence_mask``: confidence level of each photon event

        ``dist_along``: center of segment in along-track x-coordinates

        ``dt``: histogram bin size in seconds

    Keyword arguments:

        ``FIT_TYPE``: decomposition function to fit to photon height histograms

            ``'gaussian'``

            ``'general'``

        ``ITERATE``: maximum number of iterations to use in fit

        ``PEAKS``: estimated number of signal peaks in the segment histogram

        ``BACKGROUND``: vertical noise-photon density for segment


.. method:: icesat2_toolkit.fit.fit_histogram(z, hist, priors, lower_bound, upper_bound, FIT_TYPE=None)

    Optimially fit a function to the photon event histogram with Levenberg-Marquardt algorithm

    Arguments:

        ``z``: photon height histogram bins

        ``hist``: photon height histogram

        ``priors``: mean estimate for each histogram fit parameter

        ``lower_bound``: lower-bound estimate for each histogram fit parameter

        ``upper_bound``: upper-bound estimate for each histogram fit parameter

    Keyword arguments:

        ``FIT_TYPE``: decomposition function to fit to photon height histograms

            ``'gaussian'``

            ``'general'``


.. method:: icesat2_toolkit.fit.fit_geolocation(var, distance_along_X, X_atc)

    Calculate the average of photon event variables by fitting with respect to the center of the along-track coordinates

    Arguments:

        ``var``: photon event variable to compute average

        ``distance_along_X``: along-track x-coordinates

        ``X_atc``: segment center in along-track x-coordinates


.. method:: icesat2_toolkit.fit.calc_first_photon_bias(temporal_residuals, n_pulses, n_pixels, dead_time, dt, METHOD='direct', ITERATE=20)

    Estimate mean and median first photon bias corrections using segment fit residuals [Smith2019]_

    Arguments:

        ``temporal_residuals``: photon height residuals in seconds

        ``n_pulses``: estimated number of laser pulses in segment

        ``n_pixels``: number of pixels for beam

        ``dead_time``: estimated dead time

        ``dt``: histogram bin size in seconds

    Keyword arguments:

        ``METHOD``: method for computing first photon bias

            ``'direct'``

            ``'logarithmic'``

        ``ITERATE``: maximum number of iterations to use in ``'logarithmic'`` method


.. method:: icesat2_toolkit.fit.histogram_first_photon_bias(t_full, hist, n_pulses, n_pixels, dead_time, dt, METHOD='direct', ITERATE=20)

    Estimate mean and median first photon bias corrections using histogram fit residuals

    Arguments:

        ``t_full``: histogram bins in seconds

        ``hist``: photon height residuals histogram

        ``n_pulses``: estimated number of laser pulses in segment

        ``n_pixels``: number of pixels for beam

        ``dead_time``: estimated dead time

        ``dt``: histogram bin size in seconds

    Keyword arguments:

        ``METHOD``: method for computing first photon bias

            ``'direct'``

            ``'logarithmic'``

        ``ITERATE``: maximum number of iterations to use in ``'logarithmic'`` method


.. method:: icesat2_toolkit.fit.calc_transmit_pulse_shape(t_TX,p_TX,W_TX,W_RX,dt_W,SNR,ITERATE=50)

    Estimate the transmit-pulse-shape correction needed for segment averages [Smith2019]_

    Arguments:

        ``t_TX``: windowed TEP histogram time with respect to histogram centroid

        ``p_TX``: windowed TEP histogram power with noise estimate removed

        ``W_TX``: Robust Dispersion Estimate (RDE) of windowed transmit pulse

        ``W_RX``: Robust Dispersion Estimate (RDE) of segment fit residuals

        ``dt_W``: Segment fit window

        ``SNR``: Estimated signal-to-noise ratio of segment photons

    Keyword arguments:

        ``ITERATE``: maximum number of iterations to use

References
##########

.. [Pritchard2009] H. D. Pritchard et al., "Extensive dynamic thinning on the margins of the Greenland and Antarctic ice sheets", *Nature*, 461(7266), 971--975, (2009). `doi:10.1038/nature08471 <https://doi.org/10.1038/nature08471>`_
.. [Smith2017] B. E. Smith el al., "Connected subglacial lake drainage beneath Thwaites Glacier, West Antarctica", *The Cryosphere*, 11(1), 451--467, (2017). `doi:10.5194/tc-11-451-2017 <https://doi.org/10.5194/tc-11-451-2017>`_
.. [Smith2019] B. E. Smith el al., "Land ice height-retrieval algorithm for NASA's ICESat-2 photon-counting laser altimeter", *Remote Sensing of Environment*, 233, 111352, (2019). `doi:10.1016/j.rse.2019.111352 <https://doi.org/10.1016/j.rse.2019.111352>`_
