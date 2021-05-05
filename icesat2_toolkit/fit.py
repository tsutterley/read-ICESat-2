#!/usr/bin/env python
u"""
fit.py
Written by Tyler Sutterley (05/2021)
Utilities for calculating average fits from ATL03 Geolocated Photon Data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    scikit-learn: Machine Learning in Python
        http://scikit-learn.org/stable/index.html
        https://github.com/scikit-learn/scikit-learn

UPDATE HISTORY:
    Written 05/2021
"""
import operator
import itertools
import numpy as np
import scipy.stats
import scipy.signal
import scipy.optimize
import sklearn.neighbors
import icesat2_toolkit._fit as _fit

# PURPOSE: compress complete list of valid indices into a set of ranges
def compress_list(i,n):
    """
    Compress complete list of valid indices into a set of ranges

    Arguments
    ---------
    i: indices to compress
    n: largest gap between indices to accept for range
    """
    for a,b in itertools.groupby(enumerate(i), lambda v: ((v[1]-v[0])//n)*n):
        group = list(map(operator.itemgetter(1),b))
        yield (group[0], group[-1])

# PURPOSE: create distance metric for windowed classifier
def windowed_manhattan(u, v, window=[], w=None):
    """
    Create a windowed Manhattan distance metric

    Arguments
    ---------
    u: Input array
    v: Input array for distance

    Keyword arguments
    -----------------
    window: distance window for reducing neighbors
    w: weights for each value
    """
    # verify dimensions
    u = np.atleast_1d(u)
    v = np.atleast_1d(v)
    # calculate Manhattan (rectilinear) distances
    l1_diff = np.abs(u - v)
    # weight differences
    if w is not None:
        w = np.atleast_1d(w)
        l1_diff = w * l1_diff
    # broadcast window to dimensions if using a square window
    window = np.broadcast_to(np.atleast_1d(window),l1_diff.shape)
    for d,wnd in enumerate(window):
        if (l1_diff[d] >= wnd):
            l1_diff[d] = np.inf
    return l1_diff.sum()

# PURPOSE: calculate distances between points as matrices
def distance_matrix(u, v, p=1, window=[]):
    """
    Calculate distances between points as matrices

    Arguments
    ---------
    u: Input array
    v: Input array for distance

    Keyword arguments
    -----------------
    p: power for calculating distance
        1: Manhattan distances
        2: Euclidean distances
    window: distance window for reducing neighbors
    """
    M,s = np.shape(u)
    N,s = np.shape(v)
    # allocate for output distance matrix
    D = np.zeros((M,N))
    # broadcast window to dimensions if using a square window
    window = np.broadcast_to(np.atleast_1d(window),(s,))
    for d in range(s):
        ii, = np.dot(d,np.ones((1,N))).astype(np.int)
        jj, = np.dot(d,np.ones((1,M))).astype(np.int)
        dx = np.abs(u[:,ii] - v[:,jj].T)
        # window differences for dimension
        dx[dx > window[d]] = np.inf
        # add differences to total distance matrix
        D += np.power(dx,p)
    # convert distances to output units
    return np.power(D,1.0/p)

# PURPOSE: use the GSFC YAPC k-nearest neighbors algorithm to determine
# weights for each photon event within an ATL03 major frame
def classify_photons(x, h, h_win_width, indices, K=5, MIN_PH=5,
    MIN_XSPREAD=1.0, MIN_HSPREAD=0.01, METHOD='linear'):
    """
    Use the NASA GSFC YAPC k-nearest neighbors algorithm to determine
    weights for each photon event within an ATL03 major frame

    Arguments
    ---------
    x: along-track x coordinates for photon events for 3 major frames
    h: photon event heights for 3 major frames
    h_win_width: height of (possibly 2) telemetry bands
    indices: indices of photon events in ATL03 major frame

    Keyword arguments
    -----------------
    K: number of values for KNN algorithm
    MIN_PH: minimum number of photons for a major frame to be valid
    MIN_XSPREAD: minimum along-track spread of photon events
    MIN_HSPREAD: minimum window of heights for photon events
    METHOD: algorithm for computing photon event weights
        `'ball_tree'`: use scikit.learn.BallTree with custom distance metric
        `'linear'`: use a brute-force approach with linear algebra
        `'brute'`: use a brute-force approach
    """
    # number of photon events in a major frame
    n_pe = len(h[indices])
    # output photon weights for major frame
    pe_weights = np.zeros((n_pe))
    # check that number of photons is greater than criteria
    # number of points but be greater than or equal to k
    if (n_pe < MIN_PH) | (n_pe < (K+1)):
        return pe_weights
    # along-track spread of photon events
    xspread = np.max(x[indices]) - np.min(x[indices])
    # height spread of photon events
    hspread = np.max(h[indices]) - np.min(h[indices])
    # check that spread widths are greater than criteria
    if (xspread < MIN_XSPREAD) | (hspread < MIN_HSPREAD):
        return pe_weights
    # photon density
    density = n_pe/(xspread*h_win_width)
    # minimum area to contain minimum number of photon events
    area_min = MIN_PH/density
    # calculate horizontal and vertical window sizes
    win_x = 0.75*MIN_PH*np.sqrt(density)
    win_h = 0.25*MIN_PH*np.sqrt(density)
    # reduce to a buffered window around major frame
    xmin = np.min(x[indices]) - win_x
    xmax = np.max(x[indices]) + win_x
    hmin = np.min(h[indices]) - win_h
    hmax = np.max(h[indices]) + win_h
    iwin, = np.nonzero((x >= xmin) & (x <= xmax) & (h >= hmin) & (h <= hmax))
    # method of calculating photon event weights
    if (METHOD == 'ball_tree'):
        # use BallTree with custom metric to calculate photon event weights
        # window for nearest neighbors
        window = np.array([win_x/2.0,win_h/2.0])
        # create ball tree with photon events in the 3 major frames
        # using a cythonized callable distance metric
        tree = sklearn.neighbors.BallTree(np.c_[x[iwin],h[iwin]],
            metric=_fit.windowed_manhattan, window=window)
        # K nearest neighbors with windowed manhattan metric
        # use K+1 to remove identity distances (d=0)
        dist,_ = tree.query(np.c_[x[indices],h[indices]], k=K+1,
            return_distance=True)
        # calculate photon event weights and normalize by window size
        valid = np.all(np.isfinite(dist),axis=1)
        inv_dist = np.sum(win_x/2.0 + win_h/2.0 - dist[:,1:],axis=1)
        pe_weights[valid] = inv_dist[valid]/(win_x*win_h)
    elif (METHOD == 'linear'):
        # use brute force with linear algebra to calculate photon event weights
        # window for nearest neighbors
        window = np.array([win_x/2.0,win_h/2.0])
        # calculate distance matrix between points
        dist = distance_matrix(np.c_[x[indices],h[indices]],
            np.c_[x[iwin],h[iwin]], p=1, window=window)
        # sort distances and get K nearest neighbors
        # use K+1 to remove identity distances (d=0)
        dist_sort = np.sort(dist, axis=1)[:,1:K+1]
        # calculate inverse distance of photon events in window
        inv_dist = win_x/2.0 + win_h/2.0 - dist_sort
        # calculate photon event weights and normalize by window size
        valid = np.all(np.isfinite(dist_sort),axis=1)
        pe_weights[valid] = np.sum(inv_dist[valid,:],axis=1)/(win_x*win_h)
    elif (METHOD == 'brute'):
        # use brute force approach to calculate photon event weights
        # for each photon in the major frame
        for j,i in enumerate(indices):
            # all photon events in buffer excluding source photon
            ii = sorted(set(iwin) - set([i]))
            # distance of photon events to source photon
            dx = np.abs(x[ii] - x[i])
            dh = np.abs(h[ii] - h[i])
            # indices of photons within window
            n, = np.nonzero((dx < (win_x/2.0)) & (dh < (win_h/2.0)))
            # skip iteration if there are less than K within window
            if (len(n) < K):
                continue
            # calculate inverse distance of photon events in window
            inv_dist = win_x/2.0 - dx[n] + win_h/2.0 - dh[n]
            # sort distances and get K nearest neighbors
            k_sort = np.argsort(dx[n] + dh[n])[:K]
            # sum of the K largest weights (normalized by the window size)
            pe_weights[j] = np.sum(inv_dist[k_sort])/(win_x*win_h)
    # return the weights for the major frame
    return pe_weights

# PURPOSE: centers the transmit-echo-path histogram reported by ATL03
# using an iterative edit to distinguish between signal and noise
def extract_tep_histogram(tep_hist_time,tep_hist,tep_range_prim):
    """
    Centers the transmit-echo-path histogram reported by ATL03
    using an iterative edit to distinguish between signal and noise
    """
    # ATL03 recommends subset between 15-30 ns to avoid secondary
    # using primary histogram range values from ATL03 tep attributes
    i, = np.nonzero((tep_hist_time >= tep_range_prim[0]) &
        (tep_hist_time < tep_range_prim[1]))
    t_tx = np.copy(tep_hist_time[i])
    n_tx = len(t_tx)
    # noise samples of tep_hist (first 5ns and last 10 ns)
    ns,ne = (tep_range_prim[0]+5e-9,tep_range_prim[1]-10e-9)
    noise, = np.nonzero((t_tx <= ns) | (t_tx >= ne))
    noise_p1 = []
    # signal samples of tep_hist
    signal = sorted(set(np.arange(n_tx)) - set(noise))
    # number of iterations
    n_iter = 0
    while (set(noise) != set(noise_p1)) & (n_iter < 10):
        # value of noise in tep histogram
        tep_noise_value = np.sqrt(np.sum(tep_hist[i][noise]**2)/n_tx)
        p_tx = np.abs(np.copy(tep_hist[i]) - tep_noise_value)
        # calculate centroid of tep_hist
        t0_tx = np.sum(t_tx[signal]*p_tx[signal])/np.sum(p_tx[signal])
        # calculate cumulative distribution function
        TX_cpdf = np.cumsum(p_tx[signal]/np.sum(p_tx[signal]))
        # linearly interpolate to 16th and 84th percentile for RDE
        TX16,TX84 = np.interp([0.16,0.84],TX_cpdf,t_tx[signal]-t0_tx)
        # calculate width of transmitted pulse (RDE)
        W_TX = 0.5*(TX84 - TX16)
        # recalculate noise
        noise_p1 = np.copy(noise)
        ns,ne = (t0_tx-6.0*W_TX,t0_tx+6.0*W_TX)
        noise, = np.nonzero((t_tx <= ns) | (t_tx >= ne))
        signal = sorted(set(np.arange(n_tx)) - set(noise))
        # add 1 to counter
        n_iter += 1
    # valid primary TEP return has full-width at half max < 3 ns
    mx = np.argmax(p_tx[signal])
    halfmax = np.max(p_tx[signal])/2.0
    H1 = np.interp(halfmax,p_tx[signal][:mx],t_tx[signal][:mx])
    H2 = np.interp(halfmax,p_tx[signal][:mx:-1],t_tx[signal][:mx:-1])
    FWHM = H2 - H1
    # return values
    return (t_tx[signal]-t0_tx,p_tx[signal],W_TX,FWHM,ns,ne)

# PURPOSE: calculate the interquartile range (Pritchard et al, 2009) and
# robust dispersion estimator (Smith et al, 2017) of the model residuals
def filter_elevation(r0):
    """
    Calculates the interquartile range (Pritchard et al, 2009) and
    robust dispersion estimator (Smith et al, 2017) of the model residuals

    Arguments
    ---------
    r0: height residuals
    """
    # calculate percentiles for IQR and RDE
    # IQR: first and third quartiles (25th and 75th percentiles)
    # RDE: 16th and 84th percentiles
    # median: 50th percentile
    Q1,Q3,P16,P84,MEDIAN = np.percentile(r0,[25,75,16,84,50])
    # calculate interquartile range
    IQR = Q3 - Q1
    # calculate robust dispersion estimator (RDE)
    RDE = P84 - P16
    # IQR pass: residual-(median value) is within 75% of IQR
    # RDE pass: residual-(median value) is within 50% of P84-P16
    return (0.75*IQR,0.5*RDE,MEDIAN)

# PURPOSE: try fitting a surface to the signal photons with progressively
# less confidence if no valid surface is found
def try_surface_fit(x, y, z, confidence_mask, dist_along, SURF_TYPE='linear',
    ITERATE=25, CONFIDENCE=[4,3,2,1,0]):
    """
    Try fitting a surface to the signal photons with progressively
    less confidence if no valid surface is found
    """
    # try with progressively less confidence
    for i,conf in enumerate(CONFIDENCE):
        ind, = np.nonzero(confidence_mask >= conf)
        centroid = dict(x=dist_along, y=np.mean(y[ind]))
        try:
            surf = reduce_surface_fit(x[ind], y[ind], z[ind], centroid, ind,
                SURF_TYPE=SURF_TYPE, ITERATE=ITERATE)
        except (ValueError, np.linalg.linalg.LinAlgError):
            pass
        else:
            return (i+1,surf,centroid)
    # if still no values found: return infinite values
    # will need to attempt a backup algorithm
    surf = dict(error=np.full(1,np.inf))
    centroid = None
    return (None,surf,centroid)

# PURPOSE: iteratively fit a polynomial surface to the elevation data to
# reduce to within a valid window
def reduce_surface_fit(x, y, z, centroid, ind, SURF_TYPE='linear', ITERATE=25):
    """
    Iteratively fit a polynomial surface to the elevation data to reduce to
    within a valid surface window
    """
    # calculate x and y relative to centroid point
    rel_x = x - centroid['x']
    # Constant Term
    Z0 = np.ones_like((z))
    if (SURF_TYPE == 'linear'):# linear fit
        SURFMAT = np.transpose([Z0,rel_x])
    elif (SURF_TYPE == 'quadratic'):# quadratic fit
        SURFMAT = np.transpose([Z0,rel_x,rel_x**2])
    # number of points for fit and number of terms in fit
    n_max,n_terms = np.shape(SURFMAT)
    # run only if number of points is above number of terms
    FLAG1 = ((n_max - n_terms) > 10)
    # maximum allowable window size
    H_win_max = 20.0
    # minimum allowable window size
    H_win_min = 3.0
    # set initial window to the full z range
    window = z.max() - z.min()
    window_p1 = np.copy(window)
    # initial indices for reducing to window
    filt = np.arange(n_max)
    filt_p1 = np.copy(filt)
    filt_p2 = np.copy(filt_p1)
    if FLAG1:
        # save initial indices for fitting all photons for confidence level
        indices = ind.copy()
        # run fit program for polynomial type
        fit = fit_surface(x, y, z, centroid, SURF_TYPE=SURF_TYPE)
        # number of iterations performed
        n_iter = 1
        # save beta coefficients
        beta_mat = np.copy(fit['beta'])
        error_mat = np.copy(fit['error'])
        # residuals of model fit
        resid = z - np.dot(SURFMAT,beta_mat)
        # standard deviation of the residuals
        resid_std = np.std(resid)
        # save MSE and DOF for error analysis
        MSE = np.copy(fit['MSE'])
        DOF = np.copy(fit['DOF'])
        # Root mean square error
        RMSE = np.sqrt(fit['MSE'])
        # Normalized root mean square error
        NRMSE = RMSE/(np.max(z)-np.min(z))
        # IQR pass: residual-(median value) is within 75% of IQR
        # RDE pass: residual-(median value) is within 50% of P84-P16
        IQR,RDE,MEDIAN = filter_elevation(resid)
        # checking if any residuals are outside of the window
        window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
        filt, = np.nonzero(np.abs(resid-MEDIAN) <= (window/2.0))
        # save iteration of window
        window_p1 = np.copy(window)
        # run only if number of points is above number of terms
        n_rem = np.count_nonzero(np.abs(resid-MEDIAN) <= (window/2.0))
        FLAG1 = ((n_rem - n_terms) > 10)
        # maximum number of iterations to prevent infinite loops
        FLAG2 = (n_iter <= ITERATE)
        # compare indices over two iterations to prevent false stoppages
        FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))
        # iterate until there are no additional removed photons
        while FLAG1 & FLAG2 & FLAG3:
            # fit selected photons for window
            x_filt,y_filt,z_filt,indices = (x[filt],y[filt],z[filt],ind[filt])
            # run fit program for polynomial type
            fit = fit_surface(x_filt,y_filt,z_filt,centroid,SURF_TYPE=SURF_TYPE)
            # add to number of iterations performed
            n_iter += 1
            # save model coefficients
            beta_mat = np.copy(fit['beta'])
            error_mat = np.copy(fit['error'])
            # save MSE and DOF for error analysis
            MSE = np.copy(fit['MSE'])
            DOF = np.copy(fit['DOF'])
            # Root mean square error
            RMSE = np.sqrt(fit['MSE'])
            # Normalized root mean square error
            NRMSE = RMSE/(np.max(z_filt)-np.min(z_filt))
            # save number of points
            n_max = len(z_filt)
            # residuals of model fit
            resid = z - np.dot(SURFMAT,beta_mat)
            # standard deviation of the residuals
            resid_std = np.std(resid)
            # IQR pass: residual-(median value) is within 75% of IQR
            # RDE pass: residual-(median value) is within 50% of P84-P16
            IQR,RDE,MEDIAN = filter_elevation(resid)
            # checking if any residuals are outside of the window
            window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
            # filter out using median statistics and refit
            filt_p2 = np.copy(filt_p1)
            filt_p1 = np.copy(filt)
            filt, = np.nonzero(np.abs(resid-MEDIAN) <= (window/2.0))
            # save iteration of window
            window_p1 = np.copy(window)
            # run only if number of points is above number of terms
            n_rem = np.count_nonzero(np.abs(resid-MEDIAN) <= (window/2.0))
            FLAG1 = ((n_rem - n_terms) > 10)
            # maximum number of iterations to prevent infinite loops
            FLAG2 = (n_iter <= ITERATE)
            # compare indices over two iterations to prevent false stoppages
            FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))

    # return reduced model fit
    FLAG3 = (set(filt) == set(filt_p1))
    if FLAG1 & FLAG3 & (window <= H_win_max):
        return {'beta':beta_mat, 'error':error_mat, 'MSE':MSE, 'NRMSE':NRMSE,
            'DOF':DOF, 'count':n_max, 'indices':indices, 'iterations':n_iter,
            'window':window, 'RDE':RDE}
    else:
        raise ValueError('No valid data points found')

# PURPOSE: fit a polynomial surface to the elevation data
def fit_surface(x, y, z, centroid, SURF_TYPE='linear'):
    """
    Fit a polynomial surface to the elevation data
    """
    # calculate x and y relative to centroid point
    rel_x = x - centroid['x']
    # Constant Term
    Z0 = np.ones_like((z))
    # Surface design matrix
    if (SURF_TYPE == 'linear'):# linear fit
        SURFMAT = np.transpose([Z0,rel_x])
    elif (SURF_TYPE == 'quadratic'):# quadratic fit
        SURFMAT = np.transpose([Z0,rel_x,rel_x**2])
    # number of points for fit and number of terms in fit
    n_max,n_terms = np.shape(SURFMAT)
    # Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(SURFMAT,z,rcond=-1)[0]
    # modelled surface elevation
    model = np.dot(SURFMAT,beta_mat)
    # residual of fit
    res = z - model
    # nu = Degrees of Freedom = number of measurements-number of parameters
    nu = n_max - n_terms

    # Mean square error
    # MSE = (1/nu)*sum((Y-X*B)**2)
    MSE = np.dot(np.transpose(z - model),(z - model))/nu

    # elevation surface error analysis
    Hinv = np.linalg.inv(np.dot(np.transpose(SURFMAT),SURFMAT))
    # Taking the diagonal components of the cov matrix
    hdiag = np.diag(Hinv)
    # Default is 95% confidence interval
    alpha = 1.0 - (0.95)
    # Student T-Distribution with D.O.F. nu
    # t.ppf parallels tinv in matlab
    tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)
    # beta_err = t(nu,1-alpha/2)*standard error
    std_error = np.sqrt(MSE*hdiag)
    model_error = np.dot(SURFMAT,tstar*std_error)

    return {'beta':beta_mat, 'error':tstar*std_error, 'model':model,
        'model_error': model_error, 'residuals':res, 'MSE':MSE, 'DOF':nu}

# PURPOSE: try fitting a function to the signal photon histograms
# with progressively less confidence if no valid fit is found
def try_histogram_fit(x, y, z, confidence_mask, dist_along, dt,
    FIT_TYPE='gaussian', ITERATE=25, BACKGROUND=0, CONFIDENCE=[2,1,0]):
    """
    Try fitting a function to the signal photon histograms with
    progressively less confidence if no valid fit is found
    """
    # try with progressively less confidence
    for i,conf in enumerate(CONFIDENCE):
        ind, = np.nonzero(confidence_mask >= conf)
        centroid = dict(x=dist_along, y=np.mean(y[ind]))
        try:
            surf = reduce_histogram_fit(x[ind], y[ind], z[ind], ind,
                dt, FIT_TYPE=FIT_TYPE, ITERATE=ITERATE, PEAKS=2,
                BACKGROUND=BACKGROUND)
        except (ValueError, RuntimeError, SyntaxError):
            pass
        else:
            return (i+1,surf,centroid)
    # if still no values found: return infinite values
    # will need to attempt a backup algorithm
    surf = dict(error=np.full(1,np.inf))
    centroid = None
    return (None,surf,centroid)

# PURPOSE: iteratively use decomposition fitting to the elevation data to
# reduce to within a valid window
def reduce_histogram_fit(x, y, z, ind, dt, FIT_TYPE='gaussian',
    ITERATE=25, PEAKS=2, BACKGROUND=0):
    """
    Iteratively use decomposition fitting to the elevation data to reduce
    to within a valid surface window
    """
    # speed of light
    c = 299792458.0
    # use same delta time as calculating first photon bias
    # so that the residuals will be the same
    dz = dt*c
    # number of background photons in each bin
    N_BG = dz*BACKGROUND
    # create a histogram of the heights
    zmin,zmax = (z.min(),z.max())
    z_full = np.arange(zmin,zmax+dz,dz)
    nz = len(z_full)
    # maximum allowable window size
    H_win_max = 20.0
    # minimum allowable window size
    H_win_min = 3.0
    # set initial window to the full z range
    window = zmax - zmin
    window_p1 = np.copy(window)

    # number of data points
    n_max = len(z)
    # number of terms in fit
    if (FIT_TYPE == 'gaussian'):# gaussian fit
        n_terms = 3
    elif (FIT_TYPE == 'general'):# generalized gaussian fit
        n_terms = 4
    # run only if number of histogram points is above number of terms
    FLAG1 = ((nz - n_terms) > 10)

    # using kernel density functions from scikit-learn neighbors
    # gaussian kernels will reflect more accurate distributions of the data
    # with less sensitivity to sampling width than histograms (tophat kernels)
    kde = sklearn.neighbors.KernelDensity(bandwidth=dz,kernel='gaussian')
    kde.fit(z[:,None])
    # kde score_samples outputs are normalized log density functions
    hist = np.exp(kde.score_samples(z_full[:,None]) + np.log(n_max*dz))
    # smooth histogram before determining differentials
    gw = scipy.signal.gaussian(nz,4)
    hist_smooth = scipy.signal.convolve(hist, gw/gw.sum(), mode='same')
    # First differentials to find zero crossings
    # histogram 1st differential
    dhist = np.zeros((nz))
    # forward differentiation for starting point
    dhist[0] = hist_smooth[1] - hist_smooth[0]
    # backward differentiation for end point
    dhist[-1] = hist_smooth[-1] - hist_smooth[-2]
    # centered differentiation for all others
    dhist[1:-1] = (hist_smooth[2:] - hist_smooth[0:-2])/2.0

    # find positive peaks above amplitude threshold (percent of max)
    # by calculating the histogram differentials
    # signal amplitude threshold greater than 10% of max or 5.5xbackground rate
    AmpThreshold = 0.10
    HistThreshold = np.max([5.5*N_BG, AmpThreshold*np.max(hist_smooth)])
    n_peaks = np.count_nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
        ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
    n_peaks = np.min([n_peaks,PEAKS])
    peak_index, = np.nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
        ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
    # initial indices for reducing to window
    filt = np.arange(n_max)
    filt_p1 = np.copy(filt)
    filt_p2 = np.copy(filt_p1)
    if FLAG1 and (n_peaks > 0):
        # save initial indices for fitting all photons for confidence level
        indices = ind.copy()
        # sort peak index by amplitude of peaks (descending from max to min)
        # and truncate to a finite number of peaks
        sorted_peaks = np.argsort(hist[peak_index])[::-1]
        peak_index = peak_index[sorted_peaks][:n_peaks]
        # amplitude of the maximum peak
        max_amp = hist[peak_index][0]
        # cumulative probability distribution function of initial histogram
        hist_cpdf = np.cumsum(hist/np.sum(hist))
        # IQR: first and third quartiles (25th and 75th percentiles)
        # RDE: 16th and 84th percentiles
        Q1,Q3,P16,P84 = np.interp([0.25,0.75,0.16,0.84],hist_cpdf,z_full)
        # create priors list
        priors = []
        lower_bound = []
        upper_bound = []
        for i,p in enumerate(peak_index):
            if (FIT_TYPE == 'gaussian'):
                # Fit Gaussian functions to photon event histogram
                # a*: amplitude of waveform
                # r*: range from differential index
                # w*: width as 0.75*IQR
                priors.append([hist[p],z_full[p],0.75*(Q3-Q1)])
                # bounds of each parameter
                # amplitude: 0 to histogram max+5.5xbackground rate
                # range: zmin to zmax
                # width: sz to half width of z
                lower_bound.extend([0,zmin,dz])
                upper_bound.extend([max_amp+5.5*N_BG,zmax,(zmax-zmin)/2.0])
            elif (FIT_TYPE == 'general'):
                # Fit Generalized Gaussian functions to photon event histogram
                # a*: amplitude of waveform
                # r*: range from differential index
                # w*: width as 0.75*IQR
                # p*: shape parameter = gaussian sqrt(2)
                priors.append([hist[p],z_full[p],0.75*(Q3-Q1),np.sqrt(2)])
                # bounds of each parameter
                # amplitude: 0 to histogram max+5.5xbackground rate
                # range: zmin to zmax
                # width: sz to half width of z
                # shape: positive
                lower_bound.extend([0,zmin,dz,0])
                upper_bound.extend([max_amp+5.5*N_BG,zmax,(zmax-zmin)/2.0,np.inf])
        # run optimized curve fit with Levenberg-Marquardt algorithm
        fit = fit_histogram(z_full,hist,priors,lower_bound,upper_bound,
            FIT_TYPE=FIT_TYPE)
        # number of iterations performed
        n_iter = 1
        # height fits and height fit errors
        height = fit['height'].copy()
        amplitude = fit['amplitude'].copy()
        height_errors = fit['error'].copy()
        # minimum and maximum heights
        min_peak = np.min(fit['height'])
        max_peak = np.max(fit['height'])
        # save MSE and DOF for error analysis
        MSE = np.copy(fit['MSE'])
        DOF = np.copy(fit['DOF'])
        # Root mean square error
        RMSE = np.sqrt(fit['MSE'])
        # Normalized root mean square error
        NRMSE = RMSE/(zmax-zmin)
        # histogram fit
        model = np.copy(fit['model'])
        # histogram fit residuals
        resid = np.copy(fit['residuals'])
        # cumulative probability distribution function of initial histogram
        cpdf = np.cumsum(fit['residuals']/np.sum(fit['residuals']))
        # interpolate residuals to percentiles of interest for statistics
        Q1,Q3,MEDIAN,P16,P84 = np.interp([0.25,0.75,0.5,0.16,0.84],cpdf,z_full)
        # IQR: first and third quartiles (25th and 75th percentiles)
        # RDE: 16th and 84th percentiles
        IQR = 0.75*(Q3-Q1)
        RDE = 0.50*(P84-P16)
        # checking if any residuals are outside of the window
        window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
        filt, = np.nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
        # run only if number of points is above number of terms
        n_rem = np.count_nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
        nz = (np.max(z[filt])-np.min(z[filt]))//dz + 1
        FLAG1 = ((nz - n_terms) > 10) & ((n_rem - n_terms) > 10)
        # maximum number of iterations to prevent infinite loops
        FLAG2 = (n_iter <= ITERATE)
        # compare indices over two iterations to prevent false stoppages
        FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))
        # iterate until there are no additional removed photons
        while FLAG1 & FLAG2 & FLAG3:
            # fit selected photons for window
            x_filt,y_filt,z_filt,indices = (x[filt],y[filt],z[filt],ind[filt])
            zmin,zmax = (z_filt.min(),z_filt.max())
            z_full = np.arange(zmin,zmax+dz,dz)
            nz = len(z_full)
            # using kernel density functions from scikit-learn neighbors
            # gaussian kernels will reflect more accurate distributions of the data
            # with less sensitivity to sampling width than histograms (tophat kernels)
            kde = sklearn.neighbors.KernelDensity(bandwidth=dz,kernel='gaussian')
            kde.fit(z_filt[:,None])
            # kde score_samples outputs are normalized log density functions
            hist = np.exp(kde.score_samples(z_full[:,None]) + np.log(nz*dz))
            # smooth histogram before determining differentials
            gw = scipy.signal.gaussian(nz,4)
            hist_smooth = scipy.signal.convolve(hist, gw/gw.sum(), mode='same')
            # First differentials to find zero crossings
            # histogram 1st differential
            dhist = np.zeros((nz))
            # forward differentiation for starting point
            dhist[0] = hist_smooth[1] - hist_smooth[0]
            # backward differentiation for end point
            dhist[-1] = hist_smooth[-1] - hist_smooth[-2]
            # centered differentiation for all others
            dhist[1:-1] = (hist_smooth[2:] - hist_smooth[0:-2])/2.0
            # find positive peaks above amplitude threshold (percent of max)
            # by calculating the histogram differentials
            # signal amplitude threshold greater than 10% of max or 5.5xbackground rate
            HistThreshold = np.max([5.5*N_BG, AmpThreshold*np.max(hist_smooth)])
            n_peaks = np.count_nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
                ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
            n_peaks = np.min([n_peaks,PEAKS])
            peak_index, = np.nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
                ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
            # sort peak index by amplitude of peaks (descending from max to min)
            # and truncate to a finite number of peaks
            sorted_peaks = np.argsort(hist[peak_index])[::-1]
            peak_index = peak_index[sorted_peaks][:n_peaks]
            # amplitude of the maximum peak
            max_amp = hist[peak_index][0]
            # cumulative probability distribution function of initial histogram
            hist_cpdf = np.cumsum(hist/np.sum(hist))
            # IQR: first and third quartiles (25th and 75th percentiles)
            # RDE: 16th and 84th percentiles
            Q1,Q3,P16,P84 = np.interp([0.25,0.75,0.16,0.84],hist_cpdf,z_full)
            # create priors list
            priors = []
            lower_bound = []
            upper_bound = []
            for i,p in enumerate(peak_index):
                if (FIT_TYPE == 'gaussian'):
                    # Fit Gaussian functions to photon event histogram
                    # a*: amplitude of waveform
                    # r*: range from differential index
                    # w*: width as 0.75*IQR
                    priors.append([hist[p],z_full[p],0.75*(Q3-Q1)])
                    # bounds of each parameter
                    # amplitude: 0 to histogram max+5.5xbackground rate
                    # range: zmin to zmax
                    # width: sz to half width of z
                    lower_bound.extend([0,zmin,dz])
                    upper_bound.extend([max_amp+5.5*N_BG,zmax,(zmax-zmin)/2.0])
                elif (FIT_TYPE == 'general'):
                    # Fit Generalized Gaussian functions to photon event histogram
                    # a*: amplitude of waveform
                    # r*: range from differential index
                    # w*: width as 0.75*IQR
                    # p*: shape parameter = gaussian sqrt(2)
                    priors.append([hist[p],z_full[p],0.75*(Q3-Q1),np.sqrt(2)])
                    # bounds of each parameter
                    # amplitude: 0 to histogram max+5.5xbackground rate
                    # range: zmin to zmax
                    # width: sz to half width of z
                    # shape: positive
                    lower_bound.extend([0,zmin,dz,0])
                    upper_bound.extend([max_amp+5.5*N_BG,zmax,(zmax-zmin)/2.0,np.inf])
            # run optimized curve fit with Levenberg-Marquardt algorithm
            fit = fit_histogram(z_full,hist,priors,lower_bound,upper_bound,
                FIT_TYPE=FIT_TYPE)
            # add to number of iterations performed
            n_iter += 1
            # height fits and height fit errors
            height = fit['height'].copy()
            amplitude = fit['amplitude'].copy()
            height_errors = fit['error'].copy()
            # minimum and maximum heights
            min_peak = np.min(fit['height'])
            max_peak = np.max(fit['height'])
            # save MSE and DOF for error analysis
            MSE = np.copy(fit['MSE'])
            DOF = np.copy(fit['DOF'])
            # Root mean square error
            RMSE = np.sqrt(fit['MSE'])
            # Normalized root mean square error
            NRMSE = RMSE/(zmax-zmin)
            # histogram fit
            model = np.copy(fit['model'])
            # histogram fit residuals
            resid = np.copy(fit['residuals'])
            # cumulative probability distribution function of initial histogram
            cpdf = np.cumsum(resid/np.sum(resid))
            # interpolate residuals to percentiles of interest for statistics
            Q1,Q3,MEDIAN,P16,P84 = np.interp([0.25,0.75,0.5,0.16,0.84],cpdf,z_full)
            # IQR: first and third quartiles (25th and 75th percentiles)
            # RDE: 16th and 84th percentiles
            IQR = 0.75*(Q3-Q1)
            RDE = 0.50*(P84-P16)
            # checking if any residuals are outside of the window
            window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
            # filter out using median statistics and refit
            filt_p2 = np.copy(filt_p1)
            filt_p1 = np.copy(filt)
            filt, = np.nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
            # save iteration of window
            window_p1 = np.copy(window)
            # run only if number of points is above number of terms
            n_rem = np.count_nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
            nz = (np.max(z[filt])-np.min(z[filt]))//dz + 1
            FLAG1 = ((nz - n_terms) > 10) & ((n_rem - n_terms) > 10)
            # maximum number of iterations to prevent infinite loops
            FLAG2 = (n_iter <= ITERATE)
            # compare indices over two iterations to prevent false stoppages
            FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))

    # return reduced model fit
    FLAG3 = (set(filt) == set(filt_p1))
    if FLAG1 & FLAG3 & (window <= H_win_max) & (n_peaks > 0):
        # calculate time with respect to mean of fit heights
        t_full = -2*(z_full-np.mean(height))/c
        # return values
        return {'height':height, 'error':height_errors, 'amplitude':amplitude,
            'MSE':MSE, 'NRMSE':NRMSE, 'residuals':resid, 'time': t_full,
            'model':model, 'DOF':DOF, 'count':n_max, 'indices':indices,
            'iterations':n_iter, 'window':window, 'RDE':RDE, 'peaks':n_peaks}
    else:
        raise ValueError('No valid fit found')

# PURPOSE: optimially fit a function to the photon event histogram
# with Levenberg-Marquardt algorithm
def fit_histogram(z, hist, priors, lower_bound, upper_bound, FIT_TYPE=None):
    """
    Optimially fit a function to the photon event histogram with
    Levenberg-Marquardt algorithm
    """
    # create lists for the initial parameters
    # parameters, and functions for each maximum
    plist = []
    flist = []
    n_peaks = len(priors)
    # function formatting string and parameter list for each fit type
    if (FIT_TYPE == 'gaussian'):
        # summation of gaussian functions with:
        # peak amplitudes a*
        # peak ranges r* (mean)
        # peak widths w* (standard deviation)
        # Gaussian function formatting string and parameters
        function = 'a{0:d}*np.exp(-(x-r{0:d})**2.0/(2*w{0:d}**2))'
        parameters = 'a{0:d}, r{0:d}, w{0:d}'
    elif (FIT_TYPE == 'general'):
        # summation of generalized gaussian functions with:
        # peak amplitudes a*
        # peak ranges r* (mean)
        # peak widths w* (standard deviation)
        # shape parameter p* (gaussian=sqrt(2))
        # Generalized Gaussian function formatting string and parameters
        function = 'a{0:d}*np.exp(-np.abs(x-r{0:d})**(p{0:d}**2.0)/(2*w{0:d}**2))'
        parameters = 'a{0:d}, r{0:d}, w{0:d}, p{0:d}'
    # fit decomposition functions to photon event histograms
    for n,p in enumerate(priors):
        # parameter list for peak n
        plist.append(parameters.format(n))
        # function definition list for peak n
        flist.append(function.format(n))
    # initial parameters for iteration n
    p0 = np.concatenate((priors),axis=0)
    # variables for iteration n
    lambda_parameters = ', '.join([p for p in plist])
    # full function for iteration n
    lambda_function = ' + '.join([f for f in flist])
    # tuple for parameter bounds (lower and upper)
    bounds = (lower_bound, upper_bound)
    # create lambda function for iteration n
    # lambda functions are inline definitions
    # with the parameters, variables and function definition
    fsum = eval('lambda x, {0}: {1}'.format(lambda_parameters, lambda_function))
    # optimized curve fit with Levenberg-Marquardt algorithm
    # with the initial guess parameters p0 and parameter bounds
    popt, pcov = scipy.optimize.curve_fit(fsum,z,hist,p0=p0,bounds=bounds)
    # modelled histogram fit
    model = fsum(z, *popt)
    # 1 standard deviation errors in parameters
    perr = np.sqrt(np.diag(pcov))
    # number of points for fit and number of terms in fit
    n_max = len(hist)
    n_terms = len(p0)
    # extract function outputs
    if (FIT_TYPE == 'gaussian'):
        # Gaussian function outputs
        n = np.arange(n_peaks)*3
        peak_amplitude = popt[n]
        peak_height = popt[n+1]
        peak_height_error = perr[n+1]
        peak_stdev = popt[n+2]
    elif (FIT_TYPE == 'general'):
        # Generalized Gaussian function outputs
        n = np.arange(n_peaks)*4
        peak_amplitude = popt[n]
        peak_height = popt[n+1]
        peak_height_error = perr[n+1]
        peak_stdev = popt[n+2]
    # residual of fit
    res = hist - model
    # nu = Degrees of Freedom = number of measurements-number of parameters
    nu = n_max - n_terms
    # Mean square error
    # MSE = (1/nu)*sum((Y-X*B)**2)
    MSE = np.dot(np.transpose(hist - model),(hist - model))/nu
    # Default is 95% confidence interval
    alpha = 1.0 - (0.95)
    # Student T-Distribution with D.O.F. nu
    # t.ppf parallels tinv in matlab
    tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)
    return {'height':peak_height, 'amplitude':peak_amplitude,
        'error':tstar*peak_height_error, 'stdev': peak_stdev,
        'model':model, 'residuals':np.abs(res), 'MSE':MSE, 'DOF':nu}

# PURPOSE: calculate delta_time, latitude and longitude of the segment center
def fit_geolocation(var, distance_along_X, X_atc):
    """
    Calculate the average of photon event variables by fitting with respect
    to the center of the along-track coordinates
    """
    # calculate x relative to centroid point
    rel_x = distance_along_X - X_atc
    # design matrix
    XMAT = np.transpose([np.ones_like((distance_along_X)),rel_x])
    # Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(XMAT,var,rcond=-1)[0]
    # return the fitted geolocation
    return beta_mat[0]

# PURPOSE: estimate mean and median first photon bias corrections
def calc_first_photon_bias(temporal_residuals,n_pulses,n_pixels,dead_time,dt,
    METHOD='direct',ITERATE=20):
    """
    Estimate mean and median first photon bias corrections
    """
    # create a histogram of the temporal residuals
    t_full = np.arange(temporal_residuals.min(),temporal_residuals.max()+dt,dt)
    nt = len(t_full)
    # number of input photon events
    cnt = len(temporal_residuals)

    # using kernel density functions from scikit-learn neighbors
    # gaussian kernels will reflect more accurate distributions of the data
    # with less sensitivity to sampling width than histograms (tophat kernels)
    kde = sklearn.neighbors.KernelDensity(bandwidth=dt,kernel='gaussian')
    kde.fit(temporal_residuals[:,None])
    # kde score_samples outputs are normalized log density functions
    hist = np.exp(kde.score_samples(t_full[:,None]) + np.log(cnt*dt))
    N0_full = hist/(n_pulses*n_pixels)
    # centroid of initial histogram
    hist_centroid = np.sum(t_full*hist)/np.sum(hist)
    # cumulative probability distribution function of initial histogram
    hist_cpdf = np.cumsum(hist/np.sum(hist))
    # linearly interpolate to 10th, 50th, and 90th percentiles
    H10,hist_median,H90 = np.interp([0.1,0.5,0.9],hist_cpdf,t_full)

    # calculate moving total of normalized histogram
    # average number of pixels in the detector that were inactive
    P_dead = np.zeros((nt))
    # dead time as a function of the number of bins
    n_dead = int(dead_time/dt)
    # calculate moving total of last n_dead bins
    kernel = np.triu(np.tri(nt,nt,0),k=-n_dead)
    P_dead[:] = np.dot(kernel,N0_full[:,None]).flatten()

    # estimate gain directly
    if (METHOD == 'direct'):
        # estimate gain
        G_est_full = 1.0 - P_dead
        # parameters for calculating first photon bias from calibration products
        width = np.abs(H90 - H10)
        strength = np.sum(N0_full)
        # calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full
        # calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        # calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        # calculate median first photon bias correction
        # linearly interpolate to 40th, 50th and 60th percentiles
        PE40,FPB_median_corr,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    elif (METHOD == 'logarithmic') and np.count_nonzero(P_dead > 0.01):
        # find indices above threshold for computing correction
        ii, = np.nonzero(P_dead > 0.01)
        # complete gain over entire histogram
        G_est_full = np.ones((nt))
        # segment indices (above threshold and +/- dead time)
        imin,imax = (np.min(ii)-n_dead,np.max(ii)+n_dead)
        # truncate values to range of segment
        N0 = N0_full[imin:imax+1]
        N_corr = np.copy(N0)
        nr = len(N0)
        # calculate gain for segment
        gain = np.ones((nr))
        gain_prev = np.zeros((nr))
        kernel = np.triu(np.tri(nr,nr,0),k=-n_dead)
        # counter for number of iterations for segment
        n_iter = 0
        # iterate until convergence or until reaching limit of iterations
        # using matrix algebra to avoid using a nested loop
        while np.any(np.abs(gain-gain_prev) > 0.001) & (n_iter <= ITERATE):
            gain_prev=np.copy(gain)
            gain=np.exp(np.dot(kernel,np.log(1.0-N_corr[:,None]))).flatten()
            N_corr=N0/gain
            n_iter += 1
        # add segment to complete gain array
        G_est_full[imin:imax+1] = gain[:]

        # calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full

        # calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        # calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        # calculate median first photon bias correction
        # linearly interpolate to 40th, 50th and 60th percentiles
        PE40,FPB_median_corr,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    else:
        # possible that no first photon bias correction is necessary
        FPB_mean_corr = 0.0
        FPB_mean_sigma = 0.0
        FPB_median_corr = 0.0
        FPB_median_sigma = 0.0
        N_PE = np.sum(hist)

    # return first photon bias corrections
    return {'mean':FPB_mean_corr, 'mean_sigma':np.abs(FPB_mean_sigma),
        'median':FPB_median_corr, 'median_sigma':np.abs(FPB_median_sigma),
        'width':width, 'strength':strength, 'count':N_PE}

# PURPOSE: estimate mean and median first photon bias corrections
def histogram_first_photon_bias(t_full,hist,n_pulses,n_pixels,dead_time,dt,
    METHOD='direct',ITERATE=20):
    """
    Estimate mean and median first photon bias corrections using
    histogram fit residuals
    """
    # number of time points
    nt = len(t_full)
    # normalize residual histogram by number of pulses and number of pixels
    N0_full = hist/(n_pulses*n_pixels)
    # centroid of initial histogram
    hist_centroid = np.sum(t_full*hist)/np.sum(hist)
    # cumulative probability distribution function of initial histogram
    hist_cpdf = np.cumsum(hist/np.sum(hist))
    # linearly interpolate to 10th, 50th, and 90th percentiles
    H10,hist_median,H90 = np.interp([0.1,0.5,0.9],hist_cpdf,t_full)

    # calculate moving total of normalized histogram
    # average number of pixels in the detector that were inactive
    P_dead = np.zeros((nt))
    # dead time as a function of the number of bins
    n_dead = int(dead_time/dt)
    # calculate moving total of last n_dead bins
    kernel = np.triu(np.tri(nt,nt,0),k=-n_dead)
    P_dead[:] = np.dot(kernel,N0_full[:,None]).flatten()

    # estimate gain directly
    if (METHOD == 'direct'):
        # estimate gain
        G_est_full = 1.0 - P_dead
        # parameters for calculating first photon bias from calibration products
        width = np.abs(H90 - H10)
        strength = np.sum(N0_full)
        # calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full
        # calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE - hist_centroid
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        # calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        # calculate median first photon bias correction
        # linearly interpolate to 40th, 50th and 60th percentiles
        PE40,PE50,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_corr = PE50 - hist_median
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    elif (METHOD == 'logarithmic') and np.count_nonzero(P_dead > 0.01):
        # find indices above threshold for computing correction
        ii, = np.nonzero(P_dead > 0.01)
        # complete gain over entire histogram
        G_est_full = np.ones((nt))
        # segment indices (above threshold and +/- dead time)
        imin,imax = (np.min(ii)-n_dead,np.max(ii)+n_dead)
        # truncate values to range of segment
        N0 = N0_full[imin:imax+1]
        N_corr = np.copy(N0)
        nr = len(N0)
        # calculate gain for segment
        gain = np.ones((nr))
        gain_prev = np.zeros((nr))
        kernel = np.triu(np.tri(nr,nr,0),k=-n_dead)
        # counter for number of iterations for segment
        n_iter = 0
        # iterate until convergence or until reaching limit of iterations
        # using matrix algebra to avoid using a nested loop
        while np.any(np.abs(gain-gain_prev) > 0.001) & (n_iter <= ITERATE):
            gain_prev=np.copy(gain)
            gain=np.exp(np.dot(kernel,np.log(1.0-N_corr[:,None]))).flatten()
            N_corr=N0/gain
            n_iter += 1
        # add segment to complete gain array
        G_est_full[imin:imax+1] = gain[:]

        # calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full

        # calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE - hist_centroid
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        # calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        # calculate median first photon bias correction
        # linearly interpolate to 40th, 50th and 60th percentiles
        PE40,PE50,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_corr = PE50 - hist_median
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    else:
        # possible that no first photon bias correction is necessary
        FPB_mean_corr = 0.0
        FPB_mean_sigma = 0.0
        FPB_median_corr = 0.0
        FPB_median_sigma = 0.0
        N_PE = np.sum(hist)

    # return first photon bias corrections
    return {'mean':FPB_mean_corr, 'mean_sigma':np.abs(FPB_mean_sigma),
        'median':FPB_median_corr, 'median_sigma':np.abs(FPB_median_sigma),
        'width':width, 'strength':strength, 'count':N_PE}

# PURPOSE: Estimate transmit-pulse-shape correction
def calc_transmit_pulse_shape(t_TX,p_TX,W_TX,W_RX,dt_W,SNR,ITERATE=50):
    """
    Estimate the transmit-pulse-shape correction needed for segment averages
    """
    # length of the transmit pulse
    nt = len(p_TX)
    # average time step of the transmit pulse
    dt = np.abs(t_TX[1]-t_TX[0])

    # calculate broadening of the received pulse
    W_spread = np.sqrt(np.max([W_RX**2 - W_TX**2,1e-22]))
    # create zero padded transmit and received pulses (by 4*W_spread samples)
    dw = np.ceil(W_spread/dt)
    wmn = -int(np.min([0,np.round((-t_TX[0])/dt)-4*dw]))
    wmx = int(np.max([nt,np.round((-t_TX[0])/dt)+4*dw])-nt)
    t_RX = np.arange(t_TX[0]-wmn*dt,t_TX[-1]+(wmx+1)*dt,dt)
    nr = len(t_RX)
    TX = np.zeros((nr))
    TX[wmn:wmn+nt] = np.copy(p_TX)

    # smooth the transmit pulse by the spread
    gw = scipy.signal.gaussian(nr, W_spread/dt)
    RX = scipy.signal.convolve(TX/TX.sum(), gw/gw.sum(), mode='same')
    # normalize and add a random noise estimate
    RX /= np.sum(RX)
    RX += (1.0-2.0*np.random.rand(nr))*(dt/dt_W)/SNR
    # verify that all values of the synthetic received pulse are positive
    RX = np.abs(RX)

    # calculate median estimate of the synthetic received pulse
    RX_cpdf = np.cumsum(RX/np.sum(RX))
    # linearly interpolate to 50th percentile to calculate median
    t_synthetic_med = np.interp(0.5,RX_cpdf,t_RX)
    # calculate centroid for mean of the synthetic received pulse
    t_synthetic_mean = np.sum(t_RX*RX)/np.sum(RX)

    # speed of light
    c = 299792458.0
    # number of iterations
    n_iter = 0
    # threshold for stopping iteration
    threshold = 2e-4/c
    # iterate until convergence of both mean and median
    FLAG1,FLAG2 = (False,False)
    while (FLAG1 | FLAG2) and (n_iter < ITERATE):
        # copy previous mean and median times
        tmd_prev = np.copy(t_synthetic_med)
        tmn_prev = np.copy(t_synthetic_mean)
        # truncate to within window
        i, = np.nonzero((t_RX >= (t_synthetic_mean-0.5*dt_W)) &
            (t_RX <= (t_synthetic_mean+0.5*dt_W)))
        # linearly interpolate to 50th percentile to calculate median
        t_synthetic_med = np.interp(0.5,np.cumsum(RX[i]/np.sum(RX[i])),t_RX[i])
        # calculate mean time for window
        t_synthetic_mean = np.sum(t_RX[i]*RX[i])/np.sum(RX[i])
        # add to iteration
        n_iter += 1
        # check iteration
        FLAG1 = (np.abs(t_synthetic_med - tmd_prev) > threshold)
        FLAG2 = (np.abs(t_synthetic_mean - tmn_prev) > threshold)

    # return estimated transmit pulse corrections corrections
    return {'mean':t_synthetic_mean,'median':t_synthetic_med,'spread':W_spread}
