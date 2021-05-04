======
fit.py
======

Utilities for calculating average fits from ATL03 Geolocated Photon Data

Calling Sequence
================

`Source code`__

.. __: https://github.com/tsutterley/read-ICESat-2/blob/main/icesat2_toolkit/fit.py


General Methods
===============

.. method:: icesat2_toolkit.fit.compress_list(i,n)

.. method:: icesat2_toolkit.fit.classify_photons(x, h, h_win_width, indices, K=5, MIN_PH=5, MIN_XSPREAD=1.0, MIN_HSPREAD=0.01)

.. method:: icesat2_toolkit.fit.extract_tep_histogram(tep_hist_time,tep_hist,tep_range_prim)

.. method:: icesat2_toolkit.fit.filter_elevation(r0)

.. method:: icesat2_toolkit.fit.try_surface_fit(x, y, z, confidence_mask, dist_along, SURF_TYPE='linear', ITERATE=25, CONFIDENCE=[4,3,2,1,0])

.. method:: icesat2_toolkit.fit.reduce_surface_fit(x, y, z, centroid, ind, SURF_TYPE='linear', ITERATE=25)

.. method:: icesat2_toolkit.fit.fit_surface(x, y, z, centroid, SURF_TYPE='linear')

.. method:: icesat2_toolkit.fit.try_histogram_fit(x, y, z, confidence_mask, dist_along, dt, FIT_TYPE='gaussian', ITERATE=25, BACKGROUND=0, CONFIDENCE=[2,1,0])

.. method:: icesat2_toolkit.fit.reduce_histogram_fit(x, y, z, ind, dt, FIT_TYPE='gaussian', ITERATE=25, PEAKS=2, BACKGROUND=0)

.. method:: icesat2_toolkit.fit.fit_histogram(z, hist, priors, lower_bound, upper_bound, FIT_TYPE=None)

.. method:: icesat2_toolkit.fit.fit_geolocation(var, distance_along_X, X_atc)

.. method:: icesat2_toolkit.fit.calc_first_photon_bias(temporal_residuals, n_pulses, n_pixels, dead_time, dt, METHOD='direct', ITERATE=20)

.. method:: icesat2_toolkit.fit.histogram_first_photon_bias(t_full, hist, n_pulses, n_pixels, dead_time, dt, METHOD='direct', ITERATE=20)

.. method:: icesat2_toolkit.fit.calc_transmit_pulse_shape(t_TX,p_TX,W_TX,W_RX,dt_W,SNR,ITERATE=50)
