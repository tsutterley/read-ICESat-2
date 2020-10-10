#!/usr/bin/env python
u"""
MPI_ICESat2_ATL03_histogram.py (10/2020)
Read ICESat-2 ATL03 and ATL09 data files to calculate average segment surfaces
    ATL03 datasets: Global Geolocated Photons
    ATL09 datasets: Atmospheric Characteristics
Alternative algorithm that uses gaussian/generalized gaussian decomposition to
    extract possibly multiple height surfaces from a histogram of photon events

CALLING SEQUENCE:
    mpiexec -np 6 python MPI_ICESat2_ATL03_histogram.py ATL03_file ATL09_file

COMMAND LINE OPTIONS:
    -O X, --output X: Name and path of output file
    -V, --verbose: Verbose output to track progress
    -M X, --mode X: Permission mode of files created

REQUIRES MPI PROGRAM
    MPI: standardized and portable message-passing system
        https://www.open-mpi.org/
        http://mpitutorial.com/

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    mpi4py: MPI for Python
        http://pythonhosted.org/mpi4py/
        http://mpi4py.readthedocs.org/en/stable/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://h5py.org
        http://docs.h5py.org/en/stable/mpi.html
    scikit-learn: Machine Learning in Python
        http://scikit-learn.org/stable/index.html
        https://github.com/scikit-learn/scikit-learn

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    convert_calendar_decimal.py: converts from calendar date to decimal year
    time.py: Utilities for calculating time operations
    utilities: download and management utilities for syncing files

REFERENCES:
    A Chauve, C Mallet, F Bretar, S Durrieu, M P Deseilligny and W Puech.
        "Processing full-waveform lidar data: Modelling raw signals"
        International Archives of Photogrammetry
        Remote Sensing and Spatial Information Sciences, 36, 102--107 (2007)
    C Mallet and F Bretar.
        "Full-waveform topographic lidar: State-of-the-art"
        ISPRS Journal of Photogrammetry and Remote Sensing, 64, 1--16 (2009)
    J B Blair, D L Rabine, and M A Hofton.
        "The Laser Vegetation Imaging Sensor: A Medium-Altitude, Digitisation-Only,
        Airborne Laser Altimeter for Mapping Vegetation and Topography"
        ISPRS Journal of Photogrammetry and Remote Sensing, 54, 115--122 (1999)
    M A Hofton, J B Blair, J B Minster, J R Ridgway, N P Williams, J L Bufton,
        and D L Rabine. "An Airborne Scanning Laser Altimetry Survey of Long Valley,
        California" International Journal of Remote Sensing, 21(12), 2413--2437 (2000)
    M A Hofton, J B Minster, and J B Blair.
        "Decomposition of Laser Altimeter Waveforms"
        IEEE Transactions on Geoscience and Remote Sensing, 38(4), 1899--1996 (2000)
    M A Hofton, J B Blair, S B Luthcke, and D L Rabine.
        "Assessing the Performance of 20-25 m Footprint Waveform Lidar Data
        Collected in ICESat Data Corridors in Greenland"
        Geophysical Research Letters, 35, L24501 (2008), doi:10.1029/2008GL035774
    J R Ridgway, J B Minster, N Williams, J L Bufton and W B Krabill.
        "Airborne laser altimeter survey of Long Valley, California"
        Geophysical Journal International (1997) 131, 267-280

UPDATE HISTORY:
    Updated 10/2020: using argparse to set parameters
    Updated 09/2020: using reference photon delta time to interpolate ATL09
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 07/2020: "re-tiding" is no longer unnecessary
    Updated 06/2020: reduce the maximum number of peaks to fit and reduce threshold
        verify that complementary beam pair is in list of beams
        set masks of output arrays after reading from HDF5
        save histogram fit amplitudes to output HDF5 file
        add additional beam check within heights groups
    Updated 05/2020: add mean median difference of histogram fit residuals
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: adding segment quality summary variable
    Updated 04/2019: updated backup algorithm for when the histogram fit fails
        estimate both mean and median first photon bias corrections
        estimate both mean and median transmit pulse shape corrections
    Updated 03/2019: extract a set of ATL09 parameters for each ATL03 segment_ID
    Written 05/2017
"""
from __future__ import print_function, division

import sys
import os
import re
import h5py
import argparse
import datetime
import operator
import itertools
import numpy as np
import scipy.stats
import scipy.signal
import scipy.optimize
import scipy.interpolate
import sklearn.neighbors
from mpi4py import MPI
from icesat2_toolkit.convert_julian import convert_julian
from icesat2_toolkit.convert_delta_time import convert_delta_time

#-- PURPOSE: keep track of MPI threads
def info(rank, size):
    print('Rank {0:d} of {1:d}'.format(rank+1,size))
    print('module name: {0}'.format(__name__))
    if hasattr(os, 'getppid'):
        print('parent process: {0:d}'.format(os.getppid()))
    print('process id: {0:d}'.format(os.getpid()))

#-- PURPOSE: try fitting a function to the signal photons with progressively
#-- less confidence if no valid histogram fit is found
def try_histogram_fit(x, y, z, confidence_mask, dist_along, dt,
    FIT_TYPE='gaussian', ITERATE=25, BACKGROUND=0, CONFIDENCE=[2,1,0]):
    #-- try with progressively less confidence
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
    #-- if still no values found: return infinite values
    #-- will need to attempt a backup algorithm
    surf = dict(error=np.full(1,np.inf))
    centroid = None
    return (None,surf,centroid)

#-- PURPOSE: iteratively use decomposition fitting to the elevation data to
#-- reduce to within a valid window
def reduce_histogram_fit(x, y, z, ind, dt, FIT_TYPE='gaussian',
    ITERATE=25, PEAKS=2, BACKGROUND=0):
    #-- speed of light
    c = 299792458.0
    #-- use same delta time as calculating first photon bias
    #-- so that the residuals will be the same
    dz = dt*c
    #-- number of background photons in each bin
    N_BG = dz*BACKGROUND
    #-- create a histogram of the heights
    zmin,zmax = (z.min(),z.max())
    z_full = np.arange(zmin,zmax+dz,dz)
    nz = len(z_full)
    #-- maximum allowable window size
    H_win_max = 20.0
    #-- minimum allowable window size
    H_win_min = 3.0
    #-- set initial window to the full z range
    window = zmax - zmin
    window_p1 = np.copy(window)

    #-- number of data points
    n_max = len(z)
    #-- number of terms in fit
    if (FIT_TYPE == 'gaussian'):#-- gaussian fit
        n_terms = 3
    elif (FIT_TYPE == 'general'):#-- generalized gaussian fit
        n_terms = 4
    #-- run only if number of histogram points is above number of terms
    FLAG1 = ((nz - n_terms) > 10)

    #-- using kernel density functions from scikit-learn neighbors
    #-- gaussian kernels will reflect more accurate distributions of the data
    #-- with less sensitivity to sampling width than histograms (tophat kernels)
    kde = sklearn.neighbors.KernelDensity(bandwidth=dz,kernel='gaussian')
    kde.fit(z[:,None])
    #-- kde score_samples outputs are normalized log density functions
    hist = np.exp(kde.score_samples(z_full[:,None]) + np.log(n_max*dz))
    #-- smooth histogram before determining differentials
    gw = scipy.signal.gaussian(nz,4)
    hist_smooth = scipy.signal.convolve(hist, gw/gw.sum(), mode='same')
    #-- First differentials to find zero crossings
    #-- histogram 1st differential
    dhist = np.zeros((nz))
    #-- forward differentiation for starting point
    dhist[0] = hist_smooth[1] - hist_smooth[0]
    #-- backward differentiation for end point
    dhist[-1] = hist_smooth[-1] - hist_smooth[-2]
    #-- centered differentiation for all others
    dhist[1:-1] = (hist_smooth[2:] - hist_smooth[0:-2])/2.0

    #-- find positive peaks above amplitude threshold (percent of max)
    #-- by calculating the histogram differentials
    #-- signal amplitude threshold greater than 10% of max or 5.5xbackground rate
    AmpThreshold = 0.10
    HistThreshold = np.max([5.5*N_BG, AmpThreshold*np.max(hist_smooth)])
    n_peaks = np.count_nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
        ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
    n_peaks = np.min([n_peaks,PEAKS])
    peak_index, = np.nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
        ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
    #-- initial indices for reducing to window
    filt = np.arange(n_max)
    filt_p1 = np.copy(filt)
    filt_p2 = np.copy(filt_p1)
    if FLAG1 and (n_peaks > 0):
        #-- save initial indices for fitting all photons for confidence level
        indices = ind.copy()
        #-- sort peak index by amplitude of peaks (descending from max to min)
        #-- and truncate to a finite number of peaks
        sorted_peaks = np.argsort(hist[peak_index])[::-1]
        peak_index = peak_index[sorted_peaks][:n_peaks]
        #-- amplitude of the maximum peak
        max_amp = hist[peak_index][0]
        #-- estimated mean and standard deviation of peaks
        hist_mean = np.sum(hist)/nz
        hist_stdev = np.sqrt(np.sum((hist-hist_mean)**2)/nz)
        #-- cumulative probability distribution function of initial histogram
        hist_cpdf = np.cumsum(hist/np.sum(hist))
        #-- IQR: first and third quartiles (25th and 75th percentiles)
        #-- RDE: 16th and 84th percentiles
        Q1,Q3,P16,P84 = np.interp([0.25,0.75,0.16,0.84],hist_cpdf,z_full)
        #-- create priors list
        priors = []
        lower_bound = []
        upper_bound = []
        for i,p in enumerate(peak_index):
            if (FIT_TYPE == 'gaussian'):
                #-- Fit Gaussian functions to photon event histogram
                #-- a*: amplitude of waveform
                #-- r*: range from differential index
                #-- w*: width as 0.75*IQR
                priors.append([hist[p],z_full[p],0.75*(Q3-Q1)])
                #-- bounds of each parameter
                #-- amplitude: 0 to histogram max+3std
                #-- range: zmin to zmax
                #-- width: sz to half width of z
                lower_bound.extend([0,zmin,dz])
                upper_bound.extend([max_amp+3*hist_stdev,zmax,(zmax-zmin)/2.0])
            elif (FIT_TYPE == 'general'):
                #-- Fit Generalized Gaussian functions to photon event histogram
                #-- a*: amplitude of waveform
                #-- r*: range from differential index
                #-- w*: width as 0.75*IQR
                #-- p*: shape parameter = gaussian sqrt(2)
                priors.append([hist[p],z_full[p],0.75*(Q3-Q1),np.sqrt(2)])
                #-- bounds of each parameter
                #-- amplitude: 0 to histogram max+3std
                #-- range: zmin to zmax
                #-- width: sz to half width of z
                #-- shape: positive
                lower_bound.extend([0,zmin,dz,0])
                upper_bound.extend([max_amp+3*hist_stdev,zmax,(zmax-zmin)/2.0,np.inf])
        #-- run optimized curve fit with Levenberg-Marquardt algorithm
        fit = fit_histogram(z_full,hist,priors,lower_bound,upper_bound,
            FIT_TYPE=FIT_TYPE)
        #-- number of iterations performed
        n_iter = 1
        #-- height fits and height fit errors
        height = fit['height'].copy()
        amplitude = fit['amplitude'].copy()
        height_errors = fit['error'].copy()
        #-- minimum and maximum heights
        min_peak = np.min(fit['height'])
        max_peak = np.max(fit['height'])
        #-- save MSE and DOF for error analysis
        MSE = np.copy(fit['MSE'])
        DOF = np.copy(fit['DOF'])
        #-- Root mean square error
        RMSE = np.sqrt(fit['MSE'])
        #-- Normalized root mean square error
        NRMSE = RMSE/(zmax-zmin)
        #-- histogram fit
        model = np.copy(fit['model'])
        #-- histogram fit residuals
        resid = np.copy(fit['residuals'])
        #-- cumulative probability distribution function of initial histogram
        cpdf = np.cumsum(fit['residuals']/np.sum(fit['residuals']))
        #-- interpolate residuals to percentiles of interest for statistics
        Q1,Q3,MEDIAN,P16,P84 = np.interp([0.25,0.75,0.5,0.16,0.84],cpdf,z_full)
        #-- IQR: first and third quartiles (25th and 75th percentiles)
        #-- RDE: 16th and 84th percentiles
        IQR = 0.75*(Q3-Q1)
        RDE = 0.50*(P84-P16)
        #-- checking if any residuals are outside of the window
        window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
        filt, = np.nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
        #-- run only if number of points is above number of terms
        n_rem = np.count_nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
        nz = (np.max(z[filt])-np.min(z[filt]))//dz + 1
        FLAG1 = ((nz - n_terms) > 10) & ((n_rem - n_terms) > 10)
        #-- maximum number of iterations to prevent infinite loops
        FLAG2 = (n_iter <= ITERATE)
        #-- compare indices over two iterations to prevent false stoppages
        FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))
        #-- iterate until there are no additional removed photons
        while FLAG1 & FLAG2 & FLAG3:
            #-- fit selected photons for window
            x_filt,y_filt,z_filt,indices = (x[filt],y[filt],z[filt],ind[filt])
            zmin,zmax = (z_filt.min(),z_filt.max())
            z_full = np.arange(zmin,zmax+dz,dz)
            nz = len(z_full)
            #-- using kernel density functions from scikit-learn neighbors
            #-- gaussian kernels will reflect more accurate distributions of the data
            #-- with less sensitivity to sampling width than histograms (tophat kernels)
            kde = sklearn.neighbors.KernelDensity(bandwidth=dz,kernel='gaussian')
            kde.fit(z_filt[:,None])
            #-- kde score_samples outputs are normalized log density functions
            hist = np.exp(kde.score_samples(z_full[:,None]) + np.log(nz*dz))
            #-- smooth histogram before determining differentials
            gw = scipy.signal.gaussian(nz,4)
            hist_smooth = scipy.signal.convolve(hist, gw/gw.sum(), mode='same')
            #-- First differentials to find zero crossings
            #-- histogram 1st differential
            dhist = np.zeros((nz))
            #-- forward differentiation for starting point
            dhist[0] = hist_smooth[1] - hist_smooth[0]
            #-- backward differentiation for end point
            dhist[-1] = hist_smooth[-1] - hist_smooth[-2]
            #-- centered differentiation for all others
            dhist[1:-1] = (hist_smooth[2:] - hist_smooth[0:-2])/2.0
            #-- find positive peaks above amplitude threshold (percent of max)
            #-- by calculating the histogram differentials
            #-- signal amplitude threshold greater than 10% of max or 5.5xbackground rate
            HistThreshold = np.max([5.5*N_BG, AmpThreshold*np.max(hist_smooth)])
            n_peaks = np.count_nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
                ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
            n_peaks = np.min([n_peaks,PEAKS])
            peak_index, = np.nonzero((np.sign(dhist[0:-1]) >= 0) & (np.sign(dhist[1:]) < 0) &
                ((hist_smooth[0:-1] > HistThreshold) | (hist_smooth[1:] > HistThreshold)))
            #-- sort peak index by amplitude of peaks (descending from max to min)
            #-- and truncate to a finite number of peaks
            sorted_peaks = np.argsort(hist[peak_index])[::-1]
            peak_index = peak_index[sorted_peaks][:n_peaks]
            #-- amplitude of the maximum peak
            max_amp = hist[peak_index][0]
            #-- estimated mean and standard deviation of peaks
            hist_mean = np.sum(hist)/nz
            hist_stdev = np.sqrt(np.sum((hist-hist_mean)**2)/nz)
            #-- cumulative probability distribution function of initial histogram
            hist_cpdf = np.cumsum(hist/np.sum(hist))
            #-- IQR: first and third quartiles (25th and 75th percentiles)
            #-- RDE: 16th and 84th percentiles
            Q1,Q3,P16,P84 = np.interp([0.25,0.75,0.16,0.84],hist_cpdf,z_full)
            #-- create priors list
            priors = []
            lower_bound = []
            upper_bound = []
            for i,p in enumerate(peak_index):
                if (FIT_TYPE == 'gaussian'):
                    #-- Fit Gaussian functions to photon event histogram
                    #-- a*: amplitude of waveform
                    #-- r*: range from differential index
                    #-- w*: width as 0.75*IQR
                    priors.append([hist[p],z_full[p],0.75*(Q3-Q1)])
                    #-- bounds of each parameter
                    #-- amplitude: 0 to histogram max+3std
                    #-- range: zmin to zmax
                    #-- width: sz to half width of z
                    lower_bound.extend([0,zmin,dz])
                    upper_bound.extend([max_amp+3*hist_stdev,zmax,(zmax-zmin)/2.0])
                elif (FIT_TYPE == 'general'):
                    #-- Fit Generalized Gaussian functions to photon event histogram
                    #-- a*: amplitude of waveform
                    #-- r*: range from differential index
                    #-- w*: width as 0.75*IQR
                    #-- p*: shape parameter = gaussian sqrt(2)
                    priors.append([hist[p],z_full[p],0.75*(Q3-Q1),np.sqrt(2)])
                    #-- bounds of each parameter
                    #-- amplitude: 0 to histogram max+3std
                    #-- range: zmin to zmax
                    #-- width: sz to half width of z
                    #-- shape: positive
                    lower_bound.extend([0,zmin,dz,0])
                    upper_bound.extend([max_amp+3*hist_stdev,zmax,(zmax-zmin)/2.0,np.inf])
            #-- run optimized curve fit with Levenberg-Marquardt algorithm
            fit = fit_histogram(z_full,hist,priors,lower_bound,upper_bound,
                FIT_TYPE=FIT_TYPE)
            #-- add to number of iterations performed
            n_iter += 1
            #-- height fits and height fit errors
            height = fit['height'].copy()
            amplitude = fit['amplitude'].copy()
            height_errors = fit['error'].copy()
            #-- minimum and maximum heights
            min_peak = np.min(fit['height'])
            max_peak = np.max(fit['height'])
            #-- save MSE and DOF for error analysis
            MSE = np.copy(fit['MSE'])
            DOF = np.copy(fit['DOF'])
            #-- Root mean square error
            RMSE = np.sqrt(fit['MSE'])
            #-- Normalized root mean square error
            NRMSE = RMSE/(zmax-zmin)
            #-- histogram fit
            model = np.copy(fit['model'])
            #-- histogram fit residuals
            resid = np.copy(fit['residuals'])
            #-- cumulative probability distribution function of initial histogram
            cpdf = np.cumsum(resid/np.sum(resid))
            #-- interpolate residuals to percentiles of interest for statistics
            Q1,Q3,MEDIAN,P16,P84 = np.interp([0.25,0.75,0.5,0.16,0.84],cpdf,z_full)
            #-- IQR: first and third quartiles (25th and 75th percentiles)
            #-- RDE: 16th and 84th percentiles
            IQR = 0.75*(Q3-Q1)
            RDE = 0.50*(P84-P16)
            #-- checking if any residuals are outside of the window
            window = np.max([H_win_min,6.0*RDE,0.5*window_p1])
            #-- filter out using median statistics and refit
            filt_p2 = np.copy(filt_p1)
            filt_p1 = np.copy(filt)
            filt, = np.nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
            #-- save iteration of window
            window_p1 = np.copy(window)
            #-- run only if number of points is above number of terms
            n_rem = np.count_nonzero((z > (min_peak-window/2.0)) & (z < (max_peak+window/2.0)))
            nz = (np.max(z[filt])-np.min(z[filt]))//dz + 1
            FLAG1 = ((nz - n_terms) > 10) & ((n_rem - n_terms) > 10)
            #-- maximum number of iterations to prevent infinite loops
            FLAG2 = (n_iter <= ITERATE)
            #-- compare indices over two iterations to prevent false stoppages
            FLAG3 = (set(filt) != set(filt_p1)) | (set(filt_p1) != set(filt_p2))

    #-- return reduced model fit
    FLAG3 = (set(filt) == set(filt_p1))
    if FLAG1 & FLAG3 & (window <= H_win_max) & (n_peaks > 0):
        #-- calculate time with respect to height of maximum amplitude
        iamp = np.argmax(amplitude)
        t_full = -2*(z_full-height[iamp])/c
        #-- return values
        return {'height':height, 'error':height_errors, 'amplitude':amplitude,
            'MSE':MSE, 'NRMSE':NRMSE, 'residuals':resid, 'time': t_full,
            'model':model, 'DOF':DOF, 'count':n_max, 'indices':indices,
            'iterations':n_iter, 'window':window, 'RDE':RDE, 'peaks':n_peaks}
    else:
        raise ValueError('No valid fit found')

#-- PURPOSE: optimially fit a function to the photon event histogram
#-- with Levenberg-Marquardt algorithm
def fit_histogram(z, hist, priors, lower_bound, upper_bound, FIT_TYPE=None):
    #-- create lists for the initial parameters
    #-- parameters, and functions for each maximum
    plist = []
    flist = []
    n_peaks = len(priors)
    #-- function formatting string and parameter list for each fit type
    if (FIT_TYPE == 'gaussian'):
        #-- summation of gaussian functions with:
        #-- pulse amplitudes a*
        #-- pulse ranges r* (mean)
        #-- pulse widths w* (standard deviation)
        #-- Gaussian function formatting string and parameters
        function = 'a{0:d}*np.exp(-(x-r{0:d})**2.0/(2*w{0:d}**2))'
        parameters = 'a{0:d}, r{0:d}, w{0:d}'
    elif (FIT_TYPE == 'general'):
        #-- summation of generalized gaussian functions with:
        #-- pulse amplitudes a*
        #-- pulse ranges r* (mean)
        #-- pulse widths w* (standard deviation)
        #-- shape parameter p* (gaussian=sqrt(2))
        #-- Generalized Gaussian function formatting string and parameters
        function = 'a{0:d}*np.exp(-np.abs(x-r{0:d})**(p{0:d}**2.0)/(2*w{0:d}**2))'
        parameters = 'a{0:d}, r{0:d}, w{0:d}, p{0:d}'
    #-- fit decomposition functions to photon event histograms
    for n,p in enumerate(priors):
        #-- parameter list for peak n
        plist.append(parameters.format(n))
        #-- function definition list for peak n
        flist.append(function.format(n))
    #-- initial parameters for iteration n
    p0 = np.concatenate((priors),axis=0)
    #-- variables for iteration n
    lambda_parameters = ', '.join([p for p in plist])
    #-- full function for iteration n
    lambda_function = ' + '.join([f for f in flist])
    #-- tuple for parameter bounds (lower and upper)
    bounds = (lower_bound, upper_bound)
    #-- create lambda function for iteration n
    #-- lambda functions are inline definitions
    #-- with the parameters, variables and function definition
    fsum = eval('lambda x, {0}: {1}'.format(lambda_parameters, lambda_function))
    #-- optimized curve fit with Levenberg-Marquardt algorithm
    #-- with the initial guess parameters p0 and parameter bounds
    popt, pcov = scipy.optimize.curve_fit(fsum,z,hist,p0=p0,bounds=bounds)
    #-- modelled histogram fit
    model = fsum(z, *popt)
    #-- 1 standard deviation errors in parameters
    perr = np.sqrt(np.diag(pcov))
    #-- number of points for fit and number of terms in fit
    n_max = len(hist)
    n_terms = len(p0)
    #-- extract function outputs
    if (FIT_TYPE == 'gaussian'):
        #-- Gaussian function outputs
        n = np.arange(n_peaks)*3
        peak_amplitude = popt[n]
        peak_height = popt[n+1]
        peak_height_error = perr[n+1]
        peak_stdev = popt[n+2]
    elif (FIT_TYPE == 'general'):
        #-- Generalized Gaussian function outputs
        n = np.arange(n_peaks)*4
        peak_amplitude = popt[n]
        peak_height = popt[n+1]
        peak_height_error = perr[n+1]
        peak_stdev = popt[n+2]
    #-- residual of fit
    res = hist - model
    #-- nu = Degrees of Freedom = number of measurements-number of parameters
    nu = n_max - n_terms
    #-- Mean square error
    #-- MSE = (1/nu)*sum((Y-X*B)**2)
    MSE = np.dot(np.transpose(hist - model),(hist - model))/nu
    #-- Default is 95% confidence interval
    alpha = 1.0 - (0.95)
    #-- Student T-Distribution with D.O.F. nu
    #-- t.ppf parallels tinv in matlab
    tstar = scipy.stats.t.ppf(1.0-(alpha/2.0),nu)
    return {'height':peak_height, 'amplitude':peak_amplitude,
        'error':tstar*peak_height_error, 'stdev': peak_stdev,
        'model':model, 'residuals':np.abs(res), 'MSE':MSE, 'DOF':nu}

#-- PURPOSE: calculate delta_time, latitude and longitude of the segment center
def fit_geolocation(var, distance_along_X, X_atc):
    #-- calculate x relative to centroid point
    rel_x = distance_along_X - X_atc
    #-- design matrix
    XMAT = np.transpose([np.ones_like((distance_along_X)),rel_x])
    #-- Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(XMAT,var,rcond=-1)[0]
    #-- return the fitted geolocation
    return beta_mat[0]

#-- PURPOSE: estimate mean and median first photon bias corrections
def calc_first_photon_bias(t_full,hist,n_pulses,n_pixels,dead_time,dt,
    METHOD='direct',ITERATE=20):
    #-- number of time points
    nt = len(t_full)
    #-- normalize residual histogram by number of pulses and number of pixels
    N0_full = hist/(n_pulses*n_pixels)
    #-- centroid of initial histogram
    hist_centroid = np.sum(t_full*hist)/np.sum(hist)
    #-- cumulative probability distribution function of initial histogram
    hist_cpdf = np.cumsum(hist/np.sum(hist))
    #-- linearly interpolate to 10th, 50th, and 90th percentiles
    H10,hist_median,H90 = np.interp([0.1,0.5,0.9],hist_cpdf,t_full)

    #-- calculate moving total of normalized histogram
    #-- average number of pixels in the detector that were inactive
    P_dead = np.zeros((nt))
    #-- dead time as a function of the number of bins
    n_dead = np.int(dead_time/dt)
    #-- calculate moving total of last n_dead bins
    kernel = np.triu(np.tri(nt,nt,0),k=-n_dead)
    P_dead[:] = np.dot(kernel,N0_full[:,None]).flatten()

    #-- estimate gain directly
    if (METHOD == 'direct'):
        #-- estimate gain
        G_est_full = 1.0 - P_dead
        #-- parameters for calculating first photon bias from calibration products
        width = np.abs(H90 - H10)
        strength = np.sum(N0_full)
        #-- calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full
        #-- calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        #-- calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        #-- calculate median first photon bias correction
        #-- linearly interpolate to 40th, 50th and 60th percentiles
        PE40,FPB_median_corr,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    elif (METHOD == 'logarithmic') and np.count_nonzero(P_dead > 0.01):
        #-- find indices above threshold for computing correction
        ii, = np.nonzero(P_dead > 0.01)
        #-- complete gain over entire histogram
        G_est_full = np.ones((nt))
        #-- segment indices (above threshold and +/- dead time)
        imin,imax = (np.min(ii)-n_dead,np.max(ii)+n_dead)
        #-- truncate values to range of segment
        N0 = N0_full[imin:imax+1]
        N_corr = np.copy(N0)
        nr = len(N0)
        #-- calculate gain for segment
        gain = np.ones((nr))
        gain_prev = np.zeros((nr))
        kernel = np.triu(np.tri(nr,nr,0),k=-n_dead)
        #-- counter for number of iterations for segment
        n_iter = 0
        #-- iterate until convergence or until reaching limit of iterations
        #-- using matrix algebra to avoid using a nested loop
        while np.any(np.abs(gain-gain_prev) > 0.001) & (n_iter <= ITERATE):
            gain_prev=np.copy(gain)
            gain=np.exp(np.dot(kernel,np.log(1.0-N_corr[:,None]))).flatten()
            N_corr=N0/gain
            n_iter += 1
        #-- add segment to complete gain array
        G_est_full[imin:imax+1] = gain[:]

        #-- calculate corrected histogram of photon events
        N_PEcorr = (n_pulses*n_pixels)*N0_full/G_est_full
        N_PE = np.sum(N_PEcorr)
        N_sigma = np.sqrt(n_pulses*n_pixels*N0_full)/G_est_full

        #-- calculate mean corrected estimate
        FPB_mean_corr = np.sum(t_full*N_PEcorr)/N_PE
        FPB_mean_sigma = np.sqrt(np.sum((N_sigma*(t_full-FPB_mean_corr)/N_PE)**2))
        #-- calculate median corrected estimate
        PEcorr_cpdf = np.cumsum(N_PEcorr/N_PE)
        sigma_cpdf = np.sqrt(np.cumsum(N_sigma**2))/N_PE
        #-- calculate median first photon bias correction
        #-- linearly interpolate to 40th, 50th and 60th percentiles
        PE40,FPB_median_corr,PE60 = np.interp([0.4,0.5,0.6],PEcorr_cpdf,t_full)
        FPB_median_sigma = (PE60-PE40)*np.interp(0.5,PEcorr_cpdf,sigma_cpdf)/0.2
    else:
        #-- possible that no first photon bias correction is necessary
        FPB_mean_corr = 0.0
        FPB_mean_sigma = 0.0
        FPB_median_corr = 0.0
        FPB_median_sigma = 0.0
        N_PE = np.sum(hist)

    #-- return first photon bias corrections
    return {'mean':FPB_mean_corr, 'mean_sigma':np.abs(FPB_mean_sigma),
        'median':FPB_median_corr, 'median_sigma':np.abs(FPB_median_sigma),
        'width':width, 'strength':strength, 'count':N_PE}

#-- PURPOSE: compress complete list of valid indices into a set of ranges
def compress_list(i,n):
    for a,b in itertools.groupby(enumerate(i), lambda v: ((v[1]-v[0])//n)*n):
        group = list(map(operator.itemgetter(1),b))
        yield (group[0], group[-1])

#-- PURPOSE: centers the transmit-echo-path histogram reported by ATL03
#-- using an iterative edit to distinguish between signal and noise
def extract_tep_histogram(tep_hist_time,tep_hist,tep_range_prim):
    #-- ATL03 recommends subset between 15-30 ns to avoid secondary
    #-- using primary histogram range values from ATL03 tep attributes
    i, = np.nonzero((tep_hist_time >= tep_range_prim[0]) &
        (tep_hist_time < tep_range_prim[1]))
    t_tx = np.copy(tep_hist_time[i])
    n_tx = len(t_tx)
    #-- noise samples of tep_hist (first 5ns and last 10 ns)
    ns,ne = (tep_range_prim[0]+5e-9,tep_range_prim[1]-10e-9)
    noise, = np.nonzero((t_tx <= ns) | (t_tx >= ne))
    noise_p1 = []
    #-- signal samples of tep_hist
    signal = sorted(set(np.arange(n_tx)) - set(noise))
    #-- number of iterations
    n_iter = 0
    while (set(noise) != set(noise_p1)) & (n_iter < 10):
        #-- value of noise in tep histogram
        tep_noise_value = np.sqrt(np.sum(tep_hist[i][noise]**2)/n_tx)
        p_tx = np.abs(np.copy(tep_hist[i]) - tep_noise_value)
        #-- calculate centroid of tep_hist
        t0_tx = np.sum(t_tx[signal]*p_tx[signal])/np.sum(p_tx[signal])
        #-- calculate cumulative distribution function
        TX_cpdf = np.cumsum(p_tx[signal]/np.sum(p_tx[signal]))
        #-- linearly interpolate to 16th and 84th percentile for RDE
        TX16,TX84 = np.interp([0.16,0.84],TX_cpdf,t_tx[signal]-t0_tx)
        #-- calculate width of transmitted pulse (RDE)
        W_TX = 0.5*(TX84 - TX16)
        #-- recalculate noise
        noise_p1 = np.copy(noise)
        ns,ne = (t0_tx-6.0*W_TX,t0_tx+6.0*W_TX)
        noise, = np.nonzero((t_tx <= ns) | (t_tx >= ne))
        signal = sorted(set(np.arange(n_tx)) - set(noise))
        #-- add 1 to counter
        n_iter += 1
    #-- valid primary TEP return has full-width at half max < 3 ns
    mx = np.argmax(p_tx[signal])
    halfmax = np.max(p_tx[signal])/2.0
    H1 = np.interp(halfmax,p_tx[signal][:mx],t_tx[signal][:mx])
    H2 = np.interp(halfmax,p_tx[signal][:mx:-1],t_tx[signal][:mx:-1])
    FWHM = H2 - H1
    #-- return values
    return (t_tx[signal]-t0_tx,p_tx[signal],W_TX,FWHM,ns,ne)

#-- PURPOSE: Estimate transmit-pulse-shape correction
def calc_transmit_pulse_shape(t_TX,p_TX,W_TX,W_RX,dt_W,SNR,ITERATE=50):
    #-- length of the transmit pulse
    nt = len(p_TX)
    #-- average time step of the transmit pulse
    dt = np.abs(t_TX[1]-t_TX[0])

    #-- calculate broadening of the received pulse
    W_spread = np.sqrt(np.max([W_RX**2 - W_TX**2,1e-22]))
    #-- create zero padded transmit and received pulses (by 4*W_spread samples)
    dw = np.ceil(W_spread/dt)
    wmn = -np.int(np.min([0,np.round((-t_TX[0])/dt)-4*dw]))
    wmx = np.int(np.max([nt,np.round((-t_TX[0])/dt)+4*dw])-nt)
    t_RX = np.arange(t_TX[0]-wmn*dt,t_TX[-1]+(wmx+1)*dt,dt)
    nr = len(t_RX)
    TX = np.zeros((nr))
    TX[wmn:wmn+nt] = np.copy(p_TX)

    #-- smooth the transmit pulse by the spread
    gw = scipy.signal.gaussian(nr, W_spread/dt)
    RX = scipy.signal.convolve(TX/TX.sum(), gw/gw.sum(), mode='same')
    #-- normalize and add a random noise estimate
    RX /= np.sum(RX)
    RX += (1.0-2.0*np.random.rand(nr))*(dt/dt_W)/SNR
    #-- verify that all values of the synthetic received pulse are positive
    RX = np.abs(RX)

    #-- calculate median estimate of the synthetic received pulse
    RX_cpdf = np.cumsum(RX/np.sum(RX))
    #-- linearly interpolate to 50th percentile to calculate median
    t_synthetic_med = np.interp(0.5,RX_cpdf,t_RX)
    #-- calculate centroid for mean of the synthetic received pulse
    t_synthetic_mean = np.sum(t_RX*RX)/np.sum(RX)

    #-- number of iterations
    n_iter = 0
    #-- threshold for stopping iteration
    threshold = 2e-4/299792458.0
    #-- iterate until convergence of both mean and median
    FLAG1,FLAG2 = (False,False)
    while (FLAG1 | FLAG2) and (n_iter < ITERATE):
        #-- copy previous mean and median times
        tmd_prev = np.copy(t_synthetic_med)
        tmn_prev = np.copy(t_synthetic_mean)
        #-- truncate to within window
        i, = np.nonzero((t_RX >= (t_synthetic_mean-0.5*dt_W)) &
            (t_RX <= (t_synthetic_mean+0.5*dt_W)))
        #-- linearly interpolate to 50th percentile to calculate median
        t_synthetic_med = np.interp(0.5,np.cumsum(RX[i]/np.sum(RX[i])),t_RX[i])
        #-- calculate mean time for window
        t_synthetic_mean = np.sum(t_RX[i]*RX[i])/np.sum(RX[i])
        #-- add to iteration
        n_iter += 1
        #-- check iteration
        FLAG1 = (np.abs(t_synthetic_med - tmd_prev) > threshold)
        FLAG2 = (np.abs(t_synthetic_mean - tmn_prev) > threshold)

    #-- return estimated transmit pulse corrections corrections
    return {'mean':t_synthetic_mean,'median':t_synthetic_med,'spread':W_spread}

#-- PURPOSE: reads ICESat-2 ATL03 and ATL09 HDF5 files
#-- and computes heights over segments using the decomposition of histograms
def main():
    #-- start MPI communicator
    comm = MPI.COMM_WORLD

    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Read ICESat-2 ATL03 and ATL09 data files to calculate
            average segment surfaces using gaussian/generalized gaussian
            decomposition to extract possibly multiple height surfaces from a
            histogram of photon events
            """
    )
    #-- command line parameters
    #-- first file listed contains the ATL03 file
    #-- second file listed is the associated ATL09 file
    parser.add_argument('ATL03',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='ICESat-2 ATL03 file to run')
    parser.add_argument('ATL09',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='ICESat-2 ATL09 file to run')
    #-- use default output file name
    parser.add_argument('--output','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Name and path of output file')
    #-- verbosity settings
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output files')
    args = parser.parse_args()

    #-- output module information for process
    if args.verbose:
        info(comm.rank,comm.size)
    if args.verbose and (comm.rank==0):
        print('{0} -->'.format(args.ATL03))
    #-- directory setup
    ATL03_dir = os.path.dirname(args.ATL03)

    #-- compile regular expression operator for extracting data from ATL03 files
    rx1 = re.compile(r'(processed)?(ATL\d+)_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    #-- universal variables
    #-- speed of light
    c = 299792458.0
    #-- associated beam pairs
    associated_beam_pair = dict(gt1l='gt1r',gt1r='gt1l',gt2l='gt2r',gt2r='gt2l',
        gt3l='gt3r',gt3r='gt3l')

    #-- read ICESat-2 ATL03 HDF5 files (extract base parameters)
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX=rx1.findall(args.ATL03).pop()

    #-- Open the HDF5 file for reading
    fileID = h5py.File(args.ATL03, 'r', driver='mpio', comm=comm)

    #-- read each input beam within the file
    IS2_atl03_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        #-- check if subsetted beam contains data
        #-- check in both the geolocation and heights groups
        try:
            fileID[gtx]['geolocation']['segment_id']
            fileID[gtx]['heights']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl03_beams.append(gtx)

    #-- number of GPS seconds between the GPS epoch
    #-- and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = fileID['ancillary_data']['atlas_sdp_gps_epoch'][:]
    #-- which TEP to use for a given spot (convert to 0-based index)
    tep_valid_spot = fileID['ancillary_data']['tep']['tep_valid_spot'][:] - 1
    tep_pce = ['pce1_spot1','pce2_spot3']
    #-- valid range of times for each TEP histogram
    tep_range_prim = fileID['ancillary_data']['tep']['tep_range_prim'][:]
    #-- save tep parameters for a given beam
    tep = {}

    #-- variables of interest for generating corrected elevation estimates
    Segment_ID = {}
    Segment_Index_begin = {}
    Segment_PE_count = {}
    Segment_Distance = {}
    Segment_Length = {}
    Segment_Background = {}
    #-- fit parameters
    Segment_delta_time = {}
    Segment_Height = {}
    Segment_Land_Ice = {}
    Segment_Minimum = {}
    Segment_Maximum = {}
    Segment_Amplitude = {}
    Segment_Minimum_Amplitude = {}
    Segment_Maximum_Amplitude = {}
    Segment_dH_along = {}
    Segment_dH_across = {}
    Segment_Height_Error = {}
    Segment_Land_Ice_Error = {}
    Segment_Minimum_Error = {}
    Segment_Maximum_Error = {}
    Segment_dH_along_Error = {}
    Segment_dH_across_Error = {}
    Segment_Mean_Median = {}
    Segment_X_atc = {}
    Segment_X_spread = {}
    Segment_Y_atc = {}
    Segment_sigma_geo = {}
    Segment_Longitude = {}
    Segment_Latitude = {}
    Segment_N_Fit = {}
    Segment_N_Peaks = {}
    Segment_Window = {}
    Segment_RDE = {}
    Segment_SNR = {}
    Segment_Summary = {}
    Segment_Iterations = {}
    Segment_Source = {}
    Segment_Pulses = {}
    #-- correction parameters
    FPB_mean_corr = {}
    FPB_mean_sigma = {}
    FPB_median_corr = {}
    FPB_median_sigma = {}
    mean_dead_time = {}
    FPB_n_corr = {}
    FPB_cal_corr = {}
    TPS_mean_corr = {}
    TPS_median_corr = {}

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl03_beams):
        print(gtx) if args.verbose and (comm.rank == 0) else None
        #-- beam type (weak versus strong) for time
        atlas_beam_type = fileID[gtx].attrs['atlas_beam_type'].decode('utf-8')
        n_pixels = 16.0 if (atlas_beam_type == "strong") else 4.0
        #-- ATL03 Segment ID
        Segment_ID[gtx] = fileID[gtx]['geolocation']['segment_id'][:]
        #-- number of valid overlapping ATL03 segments
        n_seg = len(Segment_ID[gtx]) - 1

        #-- first photon in the segment (convert to 0-based indexing)
        Segment_Index_begin[gtx] = fileID[gtx]['geolocation']['ph_index_beg'][:] - 1
        #-- number of photon events in the segment
        Segment_PE_count[gtx] = fileID[gtx]['geolocation']['segment_ph_cnt'][:]
        #-- along-track distance for each ATL03 segment
        Segment_Distance[gtx] = fileID[gtx]['geolocation']['segment_dist_x'][:]
        #-- along-track length for each ATL03 segment
        Segment_Length[gtx] = fileID[gtx]['geolocation']['segment_length'][:]
        #-- ocean tide
        fv = fileID[gtx]['geophys_corr']['tide_ocean'].attrs['_FillValue']
        tide_ocean = np.ma.array(fileID[gtx]['geophys_corr']['tide_ocean'][:],
            fill_value=fv)
        tide_ocean.mask = tide_ocean.data == tide_ocean.fill_value
        #-- interpolate background photon rate based on 50-shot summation
        background_delta_time = fileID[gtx]['bckgrd_atlas']['delta_time'][:]
        SPL = scipy.interpolate.UnivariateSpline(background_delta_time,
            fileID[gtx]['bckgrd_atlas']['bckgrd_rate'][:],k=3,s=0)
        Segment_Background[gtx] = SPL(fileID[gtx]['geolocation']['delta_time'][:])

        #-- ATLAS spot number for beam in current orientation
        spot = np.int(fileID[gtx].attrs['atlas_spot_number'])
        #-- get ATLAS impulse response variables for the transmitter echo path (TEP)
        tep1,tep2 = ('atlas_impulse_response','tep_histogram')
        #-- get appropriate transmitter-echo-path histogram for spot
        associated_pce = tep_valid_spot[spot-1]
        pce = tep_pce[associated_pce]
        #-- delta time of TEP histogram
        tep_tod, = fileID[tep1][pce][tep2]['tep_tod'][:]
        #-- truncate tep to primary histogram (reflection 43-50 ns)
        #-- and extract signal tep from noise tep.  calculate width of tep
        #-- ATL03 recommends subsetting between 15-30 ns to avoid secondary
        tep_hist_time = np.copy(fileID[tep1][pce][tep2]['tep_hist_time'][:])
        tep_hist = np.copy(fileID[tep1][pce][tep2]['tep_hist'][:])
        t_TX,p_TX,W_TX,FWHM,TXs,TXe = extract_tep_histogram(tep_hist_time,
            tep_hist, tep_range_prim)
        #-- save tep information and statistics
        tep[gtx] = {}
        tep[gtx]['pce'] = pce
        tep[gtx]['tep_tod'] = tep_tod
        tep[gtx]['tx_start'] = TXs
        tep[gtx]['tx_end'] = TXe
        tep[gtx]['tx_robust_sprd'] = W_TX
        tep[gtx]['sigma_tx'] = FWHM

        #-- channel dead time and first photon bias table for beam
        cal1,cal2 = ('ancillary_data','calibrations')
        channel_dead_time = fileID[cal1][cal2]['dead_time'][gtx]['dead_time'][:]
        mean_dead_time[gtx] = np.mean(channel_dead_time)
        fpb_dead_time = fileID[cal1][cal2]['first_photon_bias'][gtx]['dead_time'][:]
        fpb_strength = fileID[cal1][cal2]['first_photon_bias'][gtx]['strength'][:]
        fpb_width = fileID[cal1][cal2]['first_photon_bias'][gtx]['width'][:]
        fpb_corr = fileID[cal1][cal2]['first_photon_bias'][gtx]['ffb_corr'][:]
        #-- calculate first photon bias as a function of strength and width
        #-- for the calculated mean dead time of the beam
        ndt,ns,nw = np.shape(fpb_corr)
        fpb_corr_dead_time = np.zeros((ns,nw))
        for s in range(ns):
            for w in range(nw):
                SPL = scipy.interpolate.UnivariateSpline(fpb_dead_time/1e9,
                    fpb_corr[:,s,w],k=3,s=0)
                fpb_corr_dead_time[s,w] = SPL(mean_dead_time[gtx])
        #-- bivariate spline for estimating first-photon bias using CAL-19
        CAL19 = scipy.interpolate.RectBivariateSpline(fpb_strength[0,:],
            fpb_width[0,:]/1e9, fpb_corr_dead_time/1e12, kx=1, ky=1)

        #-- allocate for output segment fit data
        fill_value = fileID[gtx]['geolocation']['sigma_h'].attrs['_FillValue']
        #-- delta time of fit photons
        Distributed_delta_time = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_delta_time.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit heights
        Distributed_Height = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Height.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Minimum = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Minimum.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Maximum = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Maximum.mask = np.ones((n_seg),dtype=np.bool)
        #-- land ice height corrected for first photon bias and transmit-pulse shape
        Distributed_Land_Ice = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Land_Ice.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit amplitudes
        Distributed_Amplitude = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Amplitude.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Minimum_Amplitude = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Minimum_Amplitude.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Maximum_Amplitude = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Maximum_Amplitude.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit along-track slopes
        Distributed_dH_along = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_along.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit height errors
        Distributed_Height_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Height_Error.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Minimum_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Minimum_Error.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_Maximum_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Maximum_Error.mask = np.ones((n_seg),dtype=np.bool)
        #-- land ice height errors (max of fit or first photon bias uncertainties)
        Distributed_Land_Ice_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Land_Ice_Error.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit along-track slope errors
        Distributed_dH_along_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_along_Error.mask = np.ones((n_seg),dtype=np.bool)
        #-- difference between the mean and median of the residuals from fit height
        Distributed_Mean_Median = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Mean_Median.mask = np.ones((n_seg),dtype=np.bool)
        #-- along-track X coordinates of segment fit
        Distributed_X_atc = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_X_atc.mask = np.ones((n_seg),dtype=np.bool)
        #-- along-track X coordinate spread of points used in segment fit
        Distributed_X_spread = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_X_spread.mask = np.ones((n_seg),dtype=np.bool)
        #-- along-track Y coordinates of segment fit
        Distributed_Y_atc = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Y_atc.mask = np.ones((n_seg),dtype=np.bool)
        #-- longitude of fit photons
        Distributed_Longitude = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Longitude.mask = np.ones((n_seg),dtype=np.bool)
        #-- latitude of fit photons
        Distributed_Latitude = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Latitude.mask = np.ones((n_seg),dtype=np.bool)
        #-- number of photons in fit
        Distributed_N_Fit = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_N_Fit.mask = np.ones((n_seg),dtype=np.bool)
        #-- number of peaks in the final histogram fit
        Distributed_N_Peaks = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_N_Peaks.mask = np.ones((n_seg),dtype=np.bool)
        #-- size of the window used in the fit
        Distributed_Window = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Window.mask = np.ones((n_seg),dtype=np.bool)
        #-- robust dispersion estimator
        Distributed_RDE = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_RDE.mask = np.ones((n_seg),dtype=np.bool)
        #-- signal-to-noise ratio
        Distributed_SNR = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_SNR.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment quality summary
        Distributed_Summary = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_Summary.mask = np.ones((n_seg),dtype=np.bool)
        #-- number of iterations for fit
        Distributed_Iterations = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_Iterations.mask = np.ones((n_seg),dtype=np.bool)
        #-- signal source selection
        Distributed_Source = np.ma.zeros((n_seg),fill_value=4,dtype=np.int)
        Distributed_Source.mask = np.ones((n_seg),dtype=np.bool)
        #-- number of pulses in segment
        Distributed_Pulses = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_Pulses.mask = np.ones((n_seg),dtype=np.bool)
        #-- first photon bias estimates
        Distributed_FPB_mean_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_mean_corr.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_FPB_mean_sigma = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_mean_sigma.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_FPB_median_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_median_corr.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_FPB_median_sigma = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_median_sigma.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_FPB_n_corr = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Distributed_FPB_n_corr.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_FPB_cal_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_cal_corr.mask = np.ones((n_seg),dtype=np.bool)
        #-- transmit pulse shape bias estimates
        Distributed_TPS_mean_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_TPS_mean_corr.mask = np.ones((n_seg),dtype=np.bool)
        Distributed_TPS_median_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_TPS_median_corr.mask = np.ones((n_seg),dtype=np.bool)

        #-- iterate over valid ATL03 segments
        #-- in ATL03 1-based indexing: invalid == 0
        #-- here in 0-based indexing: invalid == -1
        segment_indices, = np.nonzero((Segment_Index_begin[gtx][:-1] >= 0) &
            (Segment_Index_begin[gtx][1:] >= 0))
        iteration_count = len(segment_indices)
        #-- run for each geoseg (distributed over comm.size # of processes)
        for iteration in range(comm.rank, iteration_count, comm.size):
            #-- indice for iteration (can run through a subset of segments)
            j = segment_indices[iteration]

            #-- iterate over valid ATL03 segments
            #-- in ATL03 1-based indexing: invalid == 0
            #-- here in 0-based indexing: invalid == -1
            if (Segment_Index_begin[gtx][j] >= 0):
                #-- index for segment j
                idx = Segment_Index_begin[gtx][j]
                #-- number of photons in segment (use 2 ATL03 segments)
                c1 = np.int(Segment_PE_count[gtx][j])
                c2 = np.int(Segment_PE_count[gtx][j+1])
                cnt = c1 + c2
                #-- time of each Photon event (PE)
                segment_times = np.copy(fileID[gtx]['heights']['delta_time'][idx:idx+cnt])
                #-- Photon event lat/lon and elevation (re-tided WGS84)
                segment_heights = np.copy(fileID[gtx]['heights']['h_ph'][idx:idx+cnt])
                #-- ATL03 pe heights no longer apply the ocean tide
                #-- and so "re-tiding" is no longer unnecessary
                # segment_heights[:c1] += tide_ocean[j]
                # segment_heights[c1:] += tide_ocean[j+1]
                segment_lats = np.copy(fileID[gtx]['heights']['lat_ph'][idx:idx+cnt])
                segment_lons = np.copy(fileID[gtx]['heights']['lon_ph'][idx:idx+cnt])
                #-- Photon event channel and identification
                ID_channel = np.copy(fileID[gtx]['heights']['ph_id_channel'][idx:idx+cnt])
                ID_pulse = np.copy(fileID[gtx]['heights']['ph_id_pulse'][idx:idx+cnt])
                n_pulses = np.unique(ID_pulse).__len__()
                frame_number = np.copy(fileID[gtx]['heights']['pce_mframe_cnt'][idx:idx+cnt])
                #-- vertical noise-photon density
                background_density = 2.0*n_pulses*Segment_Background[gtx][j]/c
                #-- along-track X and Y coordinates
                distance_along_X = np.copy(fileID[gtx]['heights']['dist_ph_along'][idx:idx+cnt])
                distance_along_X[:c1] += Segment_Distance[gtx][j]
                distance_along_X[c1:] += Segment_Distance[gtx][j+1]
                distance_along_Y = np.copy(fileID[gtx]['heights']['dist_ph_across'][idx:idx+cnt])
                #-- check the spread of photons along-track (must be > 20m)
                along_X_spread = distance_along_X.max() - distance_along_X.min()
                #-- check confidence level associated with each photon event
                #-- -2: TEP
                #-- -1: Events not associated with a specific surface type
                #--  0: noise
                #--  1: buffer but algorithm classifies as background
                #--  2: low
                #--  3: medium
                #--  4: high
                #-- Surface types for signal classification confidence
                #-- 0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater
                ice_sig_conf = np.copy(fileID[gtx]['heights']['signal_conf_ph'][idx:idx+cnt,3])
                ice_sig_low_count = np.count_nonzero(ice_sig_conf > 1)
                #-- check if segment has photon events classified for land ice
                #-- that are at or above low-confidence threshold
                #-- and that the spread of photons is greater than 20m
                if (ice_sig_low_count > 10) & (along_X_spread > 20):
                    #-- perform a histogram fit procedure
                    Segment_X = Segment_Distance[gtx][j] + Segment_Length[gtx][j]
                    #-- step-size for histograms (50 ps ~ 7.5mm height)
                    valid,fit,centroid = try_histogram_fit(distance_along_X,
                        distance_along_Y, segment_heights, ice_sig_conf,
                        Segment_X, 5e-11, FIT_TYPE='gaussian', ITERATE=20,
                        BACKGROUND=background_density, CONFIDENCE=[2,1,0])
                    #-- indices of points used in final iterated fit
                    ifit = fit['indices'] if valid else None
                    if bool(valid) & (np.abs(fit['error'][0]) < 20):
                        #-- maximum likely height from amplitudes
                        iamp = np.argmax(fit['amplitude'])
                        Distributed_Height.data[j] = fit['height'][iamp]
                        Distributed_Height.mask[j] = False
                        Distributed_Height_Error.data[j] = fit['error'][iamp]
                        Distributed_Height_Error.mask[j] = False
                        Distributed_Amplitude.data[j] = fit['amplitude'][iamp]
                        Distributed_Amplitude.mask[j] = False
                        #-- minimum decomposed height and error
                        imin = np.argmin(fit['height'])
                        Distributed_Minimum.data[j] = fit['height'][imin]
                        Distributed_Minimum.mask[j] = False
                        Distributed_Minimum_Error.data[j] = fit['error'][imin]
                        Distributed_Minimum_Error.mask[j] = False
                        Distributed_Minimum_Amplitude.data[j] = fit['amplitude'][imin]
                        Distributed_Minimum_Amplitude.mask[j] = False
                        #-- maximum decomposed height and error
                        imax = np.argmax(fit['height'])
                        Distributed_Maximum.data[j] = fit['height'][imax]
                        Distributed_Maximum.mask[j] = False
                        Distributed_Maximum_Error.data[j] = fit['error'][imax]
                        Distributed_Maximum_Error.mask[j] = False
                        Distributed_Maximum_Amplitude.data[j] = fit['amplitude'][imax]
                        Distributed_Maximum_Amplitude.mask[j] = False
                        #-- along-track and cross-track coordinates
                        Distributed_X_atc.data[j] = np.copy(centroid['x'])
                        Distributed_X_atc.mask[j] = False
                        Distributed_X_spread.data[j] = np.copy(along_X_spread)
                        Distributed_X_spread.mask[j] = False
                        Distributed_Y_atc.data[j] = np.copy(centroid['y'])
                        Distributed_Y_atc.mask[j] = False
                        #-- fit geolocation to the along-track distance of segment
                        Distributed_delta_time.data[j] = fit_geolocation(segment_times[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_delta_time.mask[j] = False
                        Distributed_Longitude.data[j] = fit_geolocation(segment_lons[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Longitude.mask[j] = False
                        Distributed_Latitude.data[j] = fit_geolocation(segment_lats[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Latitude.mask[j] = False
                        #-- number of photons used in fit
                        Distributed_N_Fit.data[j] = len(ifit)
                        Distributed_N_Fit.mask[j] = False
                        #-- number of peaks used in histogram fit
                        Distributed_N_Peaks.data[j] = np.copy(fit['peaks'])
                        Distributed_N_Peaks.mask[j] = False
                        #-- size of the final window
                        Distributed_Window.data[j] = np.copy(fit['window'])
                        Distributed_Window.mask[j] = False
                        #-- robust dispersion estimator
                        Distributed_RDE.data[j] = np.copy(fit['RDE'])
                        Distributed_RDE.mask[j] = False
                        #-- signal to noise ratio
                        N_BG = background_density*Distributed_Window.data[j]
                        Distributed_SNR.data[j] = Distributed_N_Fit.data[j]/N_BG
                        Distributed_SNR.mask[j] = False
                        #-- number of iterations used in fit
                        Distributed_Iterations.data[j] = np.copy(fit['iterations'])
                        Distributed_Iterations.mask[j] = False
                        Distributed_Source.data[j] = np.copy(valid)
                        Distributed_Source.mask[j] = False
                        Distributed_Pulses.data[j] = np.copy(n_pulses)
                        Distributed_Pulses.mask[j] = False
                        #-- calculate difference between the mean and the median from the fit
                        z_full = -c*fit['time']/2.0
                        residual_mean = np.sum(z_full*fit['residuals'])/np.sum(fit['residuals'])
                        residual_cpdf = np.cumsum(fit['residuals']/np.sum(fit['residuals']))
                        residual_median = np.interp(0.5,residual_cpdf,z_full)
                        Distributed_Mean_Median.data[j] = residual_mean - residual_median
                        Distributed_Mean_Median.mask[j] = False
                        #-- estimate first photon bias corrections
                        #-- step-size for histograms (50 ps ~ 7.5mm height)
                        FPB = calc_first_photon_bias(fit['time'],fit['residuals'],
                            n_pulses,n_pixels,mean_dead_time[gtx],5e-11,ITERATE=20)
                        Distributed_FPB_mean_corr.data[j] = -0.5*FPB['mean']*c
                        Distributed_FPB_mean_corr.mask[j] = False
                        Distributed_FPB_mean_sigma.data[j] = 0.5*FPB['mean_sigma']*c
                        Distributed_FPB_mean_sigma.mask[j] = False
                        Distributed_FPB_median_corr.data[j] = -0.5*FPB['median']*c
                        Distributed_FPB_median_corr.mask[j] = False
                        Distributed_FPB_median_sigma.data[j] = 0.5*FPB['median_sigma']*c
                        Distributed_FPB_median_sigma.mask[j] = False
                        Distributed_FPB_n_corr.data[j] = np.copy(FPB['count'])
                        Distributed_FPB_n_corr.mask[j] = False
                        #-- first photon bias correction from CAL-19
                        FPB_calibrated = CAL19.ev(FPB['strength'],FPB['width'])
                        Distributed_FPB_cal_corr.data[j] = -0.5*FPB_calibrated*c
                        Distributed_FPB_cal_corr.mask[j] = False
                        #-- estimate transmit pulse shape correction
                        W_RX = 2.0*Distributed_RDE.data[j]/c
                        dt_W = 2.0*Distributed_Window.data[j]/c
                        TPS = calc_transmit_pulse_shape(t_TX, p_TX, W_TX, W_RX,
                            dt_W, Distributed_SNR.data[j], ITERATE=50)
                        Distributed_TPS_mean_corr.data[j] = 0.5*TPS['mean']*c
                        Distributed_TPS_mean_corr.mask[j] = False
                        Distributed_TPS_median_corr.data[j] = 0.5*TPS['median']*c
                        Distributed_TPS_median_corr.mask[j] = False
                        #-- calculate flags for quality summary
                        VPD = Distributed_N_Fit.data[j]/Distributed_Window.data[j]
                        Distributed_Summary.data[j] = np.int(
                            (Distributed_RDE.data[j] >= 1) |
                            (Distributed_Height_Error.data[j] >= 1) |
                            (VPD <= (n_pixels/4.0)))
                        Distributed_Summary.mask[j] = False

            #-- some ATL03 segments will not result in a valid fit
            #-- backup algorithm uses 4 segments to find a valid surface
            if (j not in (0,n_seg-2,n_seg-1)) & Distributed_Height.mask[j] & \
                (Segment_Index_begin[gtx][j-1] > 0):
                #-- index for segment j
                idx = Segment_Index_begin[gtx][j-1]
                #-- number of photons in segment (use 4 ATL03 segments)
                c1 = Segment_PE_count[gtx][j-1].astype(np.int)
                c2 = Segment_PE_count[gtx][j].astype(np.int)
                c3 = Segment_PE_count[gtx][j+1].astype(np.int)
                c4 = Segment_PE_count[gtx][j+2].astype(np.int)
                cnt = c1 + c2 + c3 + c4
                #-- time of each Photon event (PE)
                segment_times = np.copy(fileID[gtx]['heights']['delta_time'][idx:idx+cnt])
                #-- Photon event lat/lon and elevation (re-tided WGS84)
                segment_heights = np.copy(fileID[gtx]['heights']['h_ph'][idx:idx+cnt])
                #-- ATL03 pe heights no longer apply the ocean tide
                #-- and so "re-tiding" is no longer unnecessary
                # segment_heights[:c1] += tide_ocean[j-1]
                # segment_heights[c1:c1+c2] += tide_ocean[j]
                # segment_heights[c1+c2:c1+c2+c3] += tide_ocean[j+1]
                # segment_heights[c1+c2+c3:] += tide_ocean[j+2]
                segment_lats = np.copy(fileID[gtx]['heights']['lat_ph'][idx:idx+cnt])
                segment_lons = np.copy(fileID[gtx]['heights']['lon_ph'][idx:idx+cnt])
                #-- Photon event channel and identification
                ID_channel = np.copy(fileID[gtx]['heights']['ph_id_channel'][idx:idx+cnt])
                ID_pulse = np.copy(fileID[gtx]['heights']['ph_id_pulse'][idx:idx+cnt])
                n_pulses = np.unique(ID_pulse).__len__()
                frame_number = np.copy(fileID[gtx]['heights']['pce_mframe_cnt'][idx:idx+cnt])
                #-- vertical noise-photon density
                background_density = 2.0*n_pulses*Segment_Background[gtx][j]/c
                #-- along-track X and Y coordinates
                distance_along_X = np.copy(fileID[gtx]['heights']['dist_ph_along'][idx:idx+cnt])
                distance_along_X[:c1] += Segment_Distance[gtx][j-1]
                distance_along_X[c1:c1+c2] += Segment_Distance[gtx][j]
                distance_along_X[c1+c2:c1+c2+c3] += Segment_Distance[gtx][j+1]
                distance_along_X[c1+c2+c3:] += Segment_Distance[gtx][j+2]
                distance_along_Y = np.copy(fileID[gtx]['heights']['dist_ph_across'][idx:idx+cnt])
                #-- check the spread of photons along-track (must be > 40m)
                along_X_spread = distance_along_X.max() - distance_along_X.min()
                #-- check confidence level associated with each photon event
                #-- -2: TEP
                #-- -1: Events not associated with a specific surface type
                #--  0: noise
                #--  1: buffer but algorithm classifies as background
                #--  2: low
                #--  3: medium
                #--  4: high
                #-- Surface types for signal classification confidence
                #-- 0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater
                ice_sig_conf = np.copy(fileID[gtx]['heights']['signal_conf_ph'][idx:idx+cnt,3])
                ice_sig_low_count = np.count_nonzero(ice_sig_conf > 1)
                #-- check if segment has photon events classified for land ice
                #-- that are at or above low-confidence threshold
                #-- and that the spread of photons is greater than 40m
                if (ice_sig_low_count > 10) & (along_X_spread > 40):
                    #-- perform a histogram fit procedure
                    Segment_X = Segment_Distance[gtx][j] + Segment_Length[gtx][j]
                    #-- step-size for histograms (50 ps ~ 7.5mm height)
                    valid,fit,centroid = try_histogram_fit(distance_along_X,
                        distance_along_Y, segment_heights, ice_sig_conf,
                        Segment_X, 5e-11, FIT_TYPE='gaussian', ITERATE=20,
                        BACKGROUND=background_density, CONFIDENCE=[1,0])
                    #-- indices of points used in final iterated fit
                    ifit = fit['indices'] if valid else None
                    if bool(valid) & (np.abs(fit['error'][0]) < 20):
                        #-- maximum likely height from amplitudes
                        iamp = np.argmax(fit['amplitude'])
                        Distributed_Height.data[j] = fit['height'][iamp]
                        Distributed_Height.mask[j] = False
                        Distributed_Height_Error.data[j] = fit['error'][iamp]
                        Distributed_Height_Error.mask[j] = False
                        Distributed_Amplitude.data[j] = fit['amplitude'][iamp]
                        Distributed_Amplitude.mask[j] = False
                        #-- minimum decomposed height and error
                        imin = np.argmin(fit['height'])
                        Distributed_Minimum.data[j] = fit['height'][imin]
                        Distributed_Minimum.mask[j] = False
                        Distributed_Minimum_Error.data[j] = fit['error'][imin]
                        Distributed_Minimum_Error.mask[j] = False
                        Distributed_Minimum_Amplitude.data[j] = fit['amplitude'][imin]
                        Distributed_Minimum_Amplitude.mask[j] = False
                        #-- maximum decomposed height and error
                        imax = np.argmax(fit['height'])
                        Distributed_Maximum.data[j] = fit['height'][imax]
                        Distributed_Maximum.mask[j] = False
                        Distributed_Maximum_Error.data[j] = fit['error'][imax]
                        Distributed_Maximum_Error.mask[j] = False
                        Distributed_Maximum_Amplitude.data[j] = fit['amplitude'][imax]
                        Distributed_Maximum_Amplitude.mask[j] = False
                        #-- along-track and cross-track coordinates
                        Distributed_X_atc.data[j] = np.copy(centroid['x'])
                        Distributed_X_atc.mask[j] = False
                        Distributed_X_spread.data[j] = np.copy(along_X_spread)
                        Distributed_X_spread.mask[j] = False
                        Distributed_Y_atc.data[j] = np.copy(centroid['y'])
                        Distributed_Y_atc.mask[j] = False
                        #-- fit geolocation to the along-track distance of segment
                        Distributed_delta_time.data[j] = fit_geolocation(segment_times[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Longitude.data[j] = fit_geolocation(segment_lons[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Longitude.mask[j] = False
                        Distributed_Latitude.data[j] = fit_geolocation(segment_lats[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        #-- number of photons used in fit
                        Distributed_N_Fit.data[j] = len(ifit)
                        Distributed_N_Fit.mask[j] = False
                        #-- number of peaks used in histogram fit
                        Distributed_N_Peaks.data[j] = np.copy(fit['peaks'])
                        Distributed_N_Peaks.mask[j] = False
                        #-- size of the final window
                        Distributed_Window.data[j] = np.copy(fit['window'])
                        Distributed_Window.mask[j] = False
                        #-- robust dispersion estimator
                        Distributed_RDE.data[j] = np.copy(fit['RDE'])
                        Distributed_RDE.mask[j] = False
                        #-- signal to noise ratio
                        N_BG = background_density*Distributed_Window.data[j]
                        Distributed_SNR.data[j] = Distributed_N_Fit.data[j]/N_BG
                        Distributed_SNR.mask[j] = False
                        #-- number of iterations used in fit
                        Distributed_Iterations.data[j] = np.copy(fit['iterations'])
                        Distributed_Iterations.mask[j] = False
                        Distributed_Source.data[j] = 2 + np.copy(valid)
                        Distributed_Source.mask[j] = False
                        Distributed_Pulses.data[j] = np.copy(n_pulses)
                        Distributed_Pulses.mask[j] = False
                        #-- calculate difference between the mean and the median from the fit
                        z_full = -c*fit['time']/2.0
                        residual_mean = np.sum(z_full*fit['residuals'])/np.sum(fit['residuals'])
                        residual_cpdf = np.cumsum(fit['residuals']/np.sum(fit['residuals']))
                        residual_median = np.interp(0.5,residual_cpdf,z_full)
                        Distributed_Mean_Median.data[j] = residual_mean - residual_median
                        Distributed_Mean_Median.mask[j] = False
                        #-- estimate first photon bias corrections
                        #-- step-size for histograms (50 ps ~ 7.5mm height)
                        FPB = calc_first_photon_bias(fit['time'],fit['residuals'],
                            n_pulses,n_pixels,mean_dead_time[gtx],5e-11,ITERATE=20)
                        Distributed_FPB_mean_corr.data[j] = -0.5*FPB['mean']*c
                        Distributed_FPB_mean_corr.mask[j] = False
                        Distributed_FPB_mean_sigma.data[j] = 0.5*FPB['mean_sigma']*c
                        Distributed_FPB_mean_sigma.mask[j] = False
                        Distributed_FPB_median_corr.data[j] = -0.5*FPB['median']*c
                        Distributed_FPB_median_corr.mask[j] = False
                        Distributed_FPB_median_sigma.data[j] = 0.5*FPB['median_sigma']*c
                        Distributed_FPB_median_sigma.mask[j] = False
                        Distributed_FPB_n_corr.data[j] = np.copy(FPB['count'])
                        Distributed_FPB_n_corr.mask[j] = False
                        #-- first photon bias correction from CAL-19
                        FPB_calibrated = CAL19.ev(FPB['strength'],FPB['width'])
                        Distributed_FPB_cal_corr.data[j] = -0.5*FPB_calibrated*c
                        Distributed_FPB_cal_corr.mask[j] = False
                        #-- estimate transmit pulse shape correction
                        W_RX = 2.0*Distributed_RDE.data[j]/c
                        dt_W = 2.0*Distributed_Window.data[j]/c
                        TPS = calc_transmit_pulse_shape(t_TX, p_TX, W_TX, W_RX,
                            dt_W, Distributed_SNR.data[j], ITERATE=50)
                        Distributed_TPS_mean_corr.data[j] = 0.5*TPS['mean']*c
                        Distributed_TPS_mean_corr.mask[j] = False
                        Distributed_TPS_median_corr.data[j] = 0.5*TPS['median']*c
                        Distributed_TPS_median_corr.mask[j] = False
                        #-- calculate flags for quality summary
                        VPD = Distributed_N_Fit.data[j]/Distributed_Window.data[j]
                        Distributed_Summary.data[j] = np.int(
                            (Distributed_RDE.data[j] >= 1) |
                            (Distributed_Height_Error.data[j] >= 1) |
                            (VPD <= (n_pixels/4.0)))
                        Distributed_Summary.mask[j] = False

            #-- if there is a valid land ice height
            if (~Distributed_Height.mask[j]):
                #-- land ice height corrected for first photon bias and transmit-pulse shape
                #-- segment heights have already been "re-tided"
                Distributed_Land_Ice.data[j] = Distributed_Height.data[j] + \
                    Distributed_FPB_median_corr.data[j] + Distributed_TPS_median_corr.data[j]
                Distributed_Land_Ice.mask[j] = False
                #-- land ice height errors (max of fit or first photon bias uncertainties)
                Distributed_Land_Ice_Error.data[j] = np.sqrt(np.max([
                    Distributed_Height_Error.data[j]**2,
                    Distributed_FPB_median_sigma.data[j]**2]))
                Distributed_Land_Ice_Error.mask[j] = False

        #-- communicate output MPI matrices between ranks
        #-- operations are element summations and logical "and" across elements

        #-- delta time of fit photons
        Segment_delta_time[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_delta_time[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_delta_time.data, MPI.DOUBLE], \
            recvbuf=[Segment_delta_time[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_delta_time.mask, MPI.BOOL], \
            recvbuf=[Segment_delta_time[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_delta_time = None
        #-- segment fit heights (maximum amplitude)
        Segment_Height[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Height[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Height.data, MPI.DOUBLE], \
            recvbuf=[Segment_Height[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Height.mask, MPI.BOOL], \
            recvbuf=[Segment_Height[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Height = None
        #-- segment fit height amplitudes (maximum amplitude)
        Segment_Amplitude[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Amplitude[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Amplitude.data, MPI.DOUBLE], \
            recvbuf=[Segment_Amplitude[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Amplitude.mask, MPI.BOOL], \
            recvbuf=[Segment_Amplitude[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Amplitude = None
        #-- land ice height corrected for first photon bias and transmit-pulse shape
        Segment_Land_Ice[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Land_Ice[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice.data, MPI.DOUBLE], \
            recvbuf=[Segment_Land_Ice[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice.mask, MPI.BOOL], \
            recvbuf=[Segment_Land_Ice[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Land_Ice = None
        #-- segment fit along-track slopes
        Segment_dH_along[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_along[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_dH_along.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_along[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_along.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_along[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_along = None
        #-- segment fit height errors
        Segment_Height_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Height_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Height_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Height_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Height_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Height_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Height_Error = None
        #-- land ice height errors (max of fit or first photon bias uncertainties)
        Segment_Land_Ice_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Land_Ice_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Land_Ice_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Land_Ice_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Land_Ice_Error = None
        #-- segment fit along-track slope errors
        Segment_dH_along_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_along_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_dH_along_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_along_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_along_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_along_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_along_Error = None
        #-- segment fit heights (minimum)
        Segment_Minimum[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Minimum[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Minimum.data, MPI.DOUBLE], \
            recvbuf=[Segment_Minimum[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Minimum.mask, MPI.BOOL], \
            recvbuf=[Segment_Minimum[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Minimum = None
        #-- segment fit heights (maximum)
        Segment_Maximum[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Maximum[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Maximum.data, MPI.DOUBLE], \
            recvbuf=[Segment_Maximum[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Maximum.mask, MPI.BOOL], \
            recvbuf=[Segment_Maximum[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Maximum = None
        #-- segment fit height errors (minimum)
        Segment_Minimum_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Minimum_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Minimum_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Minimum_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Minimum_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Minimum_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Minimum_Error = None
        #-- segment fit height errors (maximum)
        Segment_Maximum_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Maximum_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Maximum_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Maximum_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Maximum_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Maximum_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Maximum_Error = None
        #-- segment fit height amplitudes (minimum)
        Segment_Minimum_Amplitude[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Minimum_Amplitude[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Minimum_Amplitude.data, MPI.DOUBLE], \
            recvbuf=[Segment_Minimum_Amplitude[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Minimum_Amplitude.mask, MPI.BOOL], \
            recvbuf=[Segment_Minimum_Amplitude[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Minimum_Amplitude = None
        #-- segment fit height amplitudes (maximum)
        Segment_Maximum_Amplitude[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Maximum_Amplitude[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Maximum_Amplitude.data, MPI.DOUBLE], \
            recvbuf=[Segment_Maximum_Amplitude[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Maximum_Amplitude.mask, MPI.BOOL], \
            recvbuf=[Segment_Maximum_Amplitude[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Maximum_Amplitude = None
        #-- difference between the mean and median of the residuals from fit height
        Segment_Mean_Median[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Mean_Median[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Mean_Median.data, MPI.DOUBLE], \
            recvbuf=[Segment_Mean_Median[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Mean_Median.mask, MPI.BOOL], \
            recvbuf=[Segment_Mean_Median[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Mean_Median = None
        #-- along-track X coordinates of segment fit
        Segment_X_atc[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_X_atc[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_X_atc.data, MPI.DOUBLE], \
            recvbuf=[Segment_X_atc[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_X_atc.mask, MPI.BOOL], \
            recvbuf=[Segment_X_atc[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_X_atc = None
        #-- along-track X coordinate spread of points used in segment fit
        Segment_X_spread[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_X_spread[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_X_spread.data, MPI.DOUBLE], \
            recvbuf=[Segment_X_spread[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_X_spread.mask, MPI.BOOL], \
            recvbuf=[Segment_X_spread[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_X_spread = None
        #-- along-track Y coordinates of segment fit
        Segment_Y_atc[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Y_atc[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Y_atc.data, MPI.DOUBLE], \
            recvbuf=[Segment_Y_atc[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Y_atc.mask, MPI.BOOL], \
            recvbuf=[Segment_Y_atc[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Y_atc = None
        #-- longitude of fit photons
        Segment_Longitude[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Longitude[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Longitude.data, MPI.DOUBLE], \
            recvbuf=[Segment_Longitude[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Longitude.mask, MPI.BOOL], \
            recvbuf=[Segment_Longitude[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Longitude = None
        #-- latitude of fit photons
        Segment_Latitude[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Latitude[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Latitude.data, MPI.DOUBLE], \
            recvbuf=[Segment_Latitude[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Latitude.mask, MPI.BOOL], \
            recvbuf=[Segment_Latitude[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Latitude = None
        #-- number of photons in fit
        Segment_N_Fit[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Segment_N_Fit[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_N_Fit.data, MPI.INT], \
            recvbuf=[Segment_N_Fit[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_N_Fit.mask, MPI.BOOL], \
            recvbuf=[Segment_N_Fit[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_N_Fit = None
        #-- number of peaks used in histogram fit
        Segment_N_Peaks[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Segment_N_Peaks[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_N_Peaks.data, MPI.INT], \
            recvbuf=[Segment_N_Peaks[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_N_Peaks.mask, MPI.BOOL], \
            recvbuf=[Segment_N_Peaks[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_N_Peaks = None
        #-- size of the window used in the fit
        Segment_Window[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Window[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Window.data, MPI.DOUBLE], \
            recvbuf=[Segment_Window[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Window.mask, MPI.BOOL], \
            recvbuf=[Segment_Window[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Window = None
        #-- robust dispersion estimator
        Segment_RDE[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_RDE[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_RDE.data, MPI.DOUBLE], \
            recvbuf=[Segment_RDE[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_RDE.mask, MPI.BOOL], \
            recvbuf=[Segment_RDE[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_RDE = None
        #-- signal-to-noise ratio
        Segment_SNR[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_SNR[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_SNR.data, MPI.DOUBLE], \
            recvbuf=[Segment_SNR[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_SNR.mask, MPI.BOOL], \
            recvbuf=[Segment_SNR[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_SNR = None
        #-- segment quality summary
        Segment_Summary[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Segment_Summary[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Summary.data, MPI.INT], \
            recvbuf=[Segment_Summary[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Summary.mask, MPI.BOOL], \
            recvbuf=[Segment_Summary[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Summary = None
        #-- number of iterations for fit
        Segment_Iterations[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Segment_Iterations[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Iterations.data, MPI.INT], \
            recvbuf=[Segment_Iterations[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Iterations.mask, MPI.BOOL], \
            recvbuf=[Segment_Iterations[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Iterations = None
        #-- signal source selection
        Segment_Source[gtx] = np.ma.zeros((n_seg),fill_value=4,dtype=np.int)
        Segment_Source[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Source.data, MPI.INT], \
            recvbuf=[Segment_Source[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Source.mask, MPI.BOOL], \
            recvbuf=[Segment_Source[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Source = None
        #-- number of pulses in segment
        Segment_Pulses[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        Segment_Pulses[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_Pulses.data, MPI.INT], \
            recvbuf=[Segment_Pulses[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Pulses.mask, MPI.BOOL], \
            recvbuf=[Segment_Pulses[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Pulses = None
        #-- first photon bias estimates
        FPB_mean_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_mean_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_mean_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_mean_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_mean_corr = None
        FPB_mean_sigma[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_mean_sigma[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_sigma.data, MPI.DOUBLE], \
            recvbuf=[FPB_mean_sigma[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_sigma.mask, MPI.BOOL], \
            recvbuf=[FPB_mean_sigma[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_mean_sigma = None
        FPB_median_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_median_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_median_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_median_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_median_corr = None
        FPB_median_sigma[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_median_sigma[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_sigma.data, MPI.DOUBLE], \
            recvbuf=[FPB_median_sigma[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_sigma.mask, MPI.BOOL], \
            recvbuf=[FPB_median_sigma[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_median_sigma = None
        FPB_n_corr[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=np.int)
        FPB_n_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_n_corr.data, MPI.INT], \
            recvbuf=[FPB_n_corr[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_n_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_n_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_n_corr = None
        FPB_cal_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_cal_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_cal_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_cal_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_cal_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_cal_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_cal_corr = None
        #-- transmit pulse shape bias estimates
        TPS_mean_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        TPS_mean_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_TPS_mean_corr.data, MPI.DOUBLE], \
            recvbuf=[TPS_mean_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_TPS_mean_corr.mask, MPI.BOOL], \
            recvbuf=[TPS_mean_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_TPS_mean_corr = None
        TPS_median_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        TPS_median_corr[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_TPS_median_corr.data, MPI.DOUBLE], \
            recvbuf=[TPS_median_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_TPS_median_corr.mask, MPI.BOOL], \
            recvbuf=[TPS_median_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_TPS_median_corr = None
        #-- wait for all distributed processes to finish for beam
        comm.Barrier()

    #-- copy variables for outputting to HDF5 file
    IS2_atl03_fit = {}
    IS2_atl03_fill = {}
    IS2_atl03_attrs = {}

    #-- ICESat-2 spacecraft orientation at time
    IS2_atl03_fit['orbit_info'] = {}
    IS2_atl03_attrs['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl03_fit['orbit_info'][key] = val[:]
        #-- Getting attributes of group and included variables
        #-- Global Group Attributes
        for att_name,att_val in fileID['orbit_info'].attrs.items():
            IS2_atl03_attrs['orbit_info'][att_name] = att_val
        #-- Variable Attributes
        IS2_atl03_attrs['orbit_info'][key] = {}
        for att_name,att_val in val.attrs.items():
            IS2_atl03_attrs['orbit_info'][key][att_name] = att_val

    #-- information ancillary to the data product
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl03_fit['ancillary_data'] = {}
    IS2_atl03_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc','end_cycle',
        'end_geoseg','end_gpssow','end_gpsweek','end_orbit','end_region',
        'end_rgt','granule_end_utc','granule_start_utc','release','start_cycle',
        'start_geoseg','start_gpssow','start_gpsweek','start_orbit','start_region',
        'start_rgt','version']:
        #-- get each HDF5 variable
        IS2_atl03_fit['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        IS2_atl03_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl03_attrs['ancillary_data'][key][att_name] = att_val

    #-- for each output beam
    for gtx in sorted(IS2_atl03_beams):
        #-- atmospheric profile for beam gtx from ATL09 dataset
        pfl = fileID[gtx].attrs['atmosphere_profile']
        #-- complementary beam in pair
        cmp = associated_beam_pair[gtx]
        #-- extract and interpolate atmospheric parameters from ATL09
        dtime = fileID[gtx]['geolocation']['delta_time'][:]
        IS2_atl09_mds,IS2_atl09_attrs = read_HDF5_ATL09(args.ATL09, pfl,
            dtime, ATTRIBUTES=True, VERBOSE=args.verbose, COMM=comm)

        #-- segment fit across-track slopes
        Distributed_dH_across = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_across.mask = np.ones((n_seg),dtype=np.bool)
        #-- segment fit across-track slope errors
        Distributed_dH_across_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_across_Error.mask = np.ones((n_seg),dtype=np.bool)
        #-- contribution of geolocation uncertainty to height error
        Distributed_sigma_geo = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_sigma_geo.mask = np.ones((n_seg),dtype=np.bool)

        #-- iterate over valid ATL03 segments
        #-- in ATL03 1-based indexing: invalid == 0
        #-- here in 0-based indexing: invalid == -1
        segment_indices, = np.nonzero((Segment_Index_begin[gtx][:-1] >= 0) &
            (Segment_Index_begin[gtx][1:] >= 0))
        #-- verify that complementary beam pair is in list of beams
        iteration_count = len(segment_indices) if (cmp in IS2_atl03_beams) else 0
        #-- run for each geoseg (distributed over comm.size # of processes)
        for iteration in range(comm.rank, iteration_count, comm.size):
            #-- indice for iteration (can run through a subset of segments)
            j = segment_indices[iteration]
            #-- across track slopes for beam
            if ((~Segment_Height[gtx].mask[j]) & (~Segment_Height[cmp].mask[j])):
                #-- segment fit across-track slopes
                dY = (Segment_Y_atc[gtx].data[j] - Segment_Y_atc[cmp].data[j])
                Distributed_dH_across.data[j] = (Segment_Land_Ice[gtx].data[j] -
                    Segment_Land_Ice[cmp].data[j])/dY
                Distributed_dH_across.mask[j] = False
                #-- segment fit across-track slope errors
                Distributed_dH_across_Error.data[j] = np.sqrt(
                    Segment_Land_Ice_Error[gtx].data[j]**2 +
                    Segment_Land_Ice_Error[cmp].data[j]**2)/np.abs(dY)
                Distributed_dH_across_Error.mask[j] = False
                #-- geolocation uncertainty
                sigma_geo_across = fileID[gtx]['geolocation']['sigma_across'][j]
                sigma_geo_along = fileID[gtx]['geolocation']['sigma_along'][j]
                sigma_geo_h = fileID[gtx]['geolocation']['sigma_h'][j]
                #-- contribution of geolocation uncertainty to height errors
                Distributed_sigma_geo.data[j] = np.sqrt(sigma_geo_h**2 +
                    (sigma_geo_along*Segment_dH_along[gtx].data[j])**2 +
                    (sigma_geo_across*Distributed_dH_across.data[j])**2)
                Distributed_sigma_geo.mask[j] = False

        #-- segment fit across-track slopes
        Segment_dH_across[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_across[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_dH_across.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_across[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_across.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_across[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_across = None
        #-- segment fit across-track slope errors
        Segment_dH_across_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_across_Error[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_dH_across_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_across_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_across_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_across_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_across_Error = None
        #-- contribution of geolocation uncertainty to height errors
        Segment_sigma_geo[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_sigma_geo[gtx].mask = np.ones((n_seg),dtype=np.bool)
        comm.Allreduce(sendbuf=[Distributed_sigma_geo.data, MPI.DOUBLE], \
            recvbuf=[Segment_sigma_geo[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_sigma_geo.mask, MPI.BOOL], \
            recvbuf=[Segment_sigma_geo[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_sigma_geo = None
        #-- wait for all distributed processes to finish for beam
        comm.Barrier()

        #-- set values for invalid segments to fill_value of each variable
        Segment_delta_time[gtx].data[Segment_delta_time[gtx].mask] = Segment_delta_time[gtx].fill_value
        Segment_Height[gtx].data[Segment_Height[gtx].mask] = Segment_Height[gtx].fill_value
        Segment_Minimum[gtx].data[Segment_Minimum[gtx].mask] = Segment_Minimum[gtx].fill_value
        Segment_Maximum[gtx].data[Segment_Maximum[gtx].mask] = Segment_Maximum[gtx].fill_value
        Segment_Land_Ice[gtx].data[Segment_Land_Ice[gtx].mask] = Segment_Land_Ice[gtx].fill_value
        Segment_dH_along[gtx].data[Segment_dH_along[gtx].mask] = Segment_dH_along[gtx].fill_value
        Segment_dH_across[gtx].data[Segment_dH_across[gtx].mask] = Segment_dH_across[gtx].fill_value
        Segment_Height_Error[gtx].data[Segment_Height_Error[gtx].mask] = Segment_Height_Error[gtx].fill_value
        Segment_Minimum_Error[gtx].data[Segment_Minimum_Error[gtx].mask] = Segment_Minimum_Error[gtx].fill_value
        Segment_Maximum_Error[gtx].data[Segment_Maximum_Error[gtx].mask] = Segment_Maximum_Error[gtx].fill_value
        Segment_Land_Ice_Error[gtx].data[Segment_Land_Ice_Error[gtx].mask] = Segment_Land_Ice_Error[gtx].fill_value
        Segment_Amplitude[gtx].data[Segment_Amplitude[gtx].mask] = Segment_Amplitude[gtx].fill_value
        Segment_Minimum_Amplitude[gtx].data[Segment_Minimum_Amplitude[gtx].mask] = Segment_Minimum_Amplitude[gtx].fill_value
        Segment_Maximum_Amplitude[gtx].data[Segment_Maximum_Amplitude[gtx].mask] = Segment_Maximum_Amplitude[gtx].fill_value
        Segment_dH_along_Error[gtx].data[Segment_dH_along_Error[gtx].mask] = Segment_dH_along_Error[gtx].fill_value
        Segment_dH_across_Error[gtx].data[Segment_dH_across_Error[gtx].mask] = Segment_dH_across_Error[gtx].fill_value
        Segment_Mean_Median[gtx].data[Segment_Mean_Median[gtx].mask] = Segment_Mean_Median[gtx].fill_value
        Segment_X_atc[gtx].data[Segment_X_atc[gtx].mask] = Segment_X_atc[gtx].fill_value
        Segment_X_spread[gtx].data[Segment_X_spread[gtx].mask] = Segment_X_spread[gtx].fill_value
        Segment_Y_atc[gtx].data[Segment_Y_atc[gtx].mask] = Segment_Y_atc[gtx].fill_value
        Segment_sigma_geo[gtx].data[Segment_sigma_geo[gtx].mask] = Segment_sigma_geo[gtx].fill_value
        Segment_Longitude[gtx].data[Segment_Longitude[gtx].mask] = Segment_Longitude[gtx].fill_value
        Segment_Latitude[gtx].data[Segment_Latitude[gtx].mask] = Segment_Latitude[gtx].fill_value
        Segment_N_Fit[gtx].data[Segment_N_Fit[gtx].mask] = Segment_N_Fit[gtx].fill_value
        Segment_N_Peaks[gtx].data[Segment_N_Peaks[gtx].mask] = Segment_N_Peaks[gtx].fill_value
        Segment_Window[gtx].data[Segment_Window[gtx].mask] = Segment_Window[gtx].fill_value
        Segment_RDE[gtx].data[Segment_RDE[gtx].mask] = Segment_RDE[gtx].fill_value
        Segment_SNR[gtx].data[Segment_SNR[gtx].mask] = Segment_SNR[gtx].fill_value
        Segment_Summary[gtx].data[Segment_Summary[gtx].mask] = Segment_Summary[gtx].fill_value
        Segment_Iterations[gtx].data[Segment_Iterations[gtx].mask] = Segment_Iterations[gtx].fill_value
        Segment_Source[gtx].data[Segment_Source[gtx].mask] = Segment_Source[gtx].fill_value
        Segment_Pulses[gtx].data[Segment_Pulses[gtx].mask] = Segment_Pulses[gtx].fill_value
        FPB_mean_corr[gtx].data[FPB_mean_corr[gtx].mask] = FPB_mean_corr[gtx].fill_value
        FPB_mean_sigma[gtx].data[FPB_mean_sigma[gtx].mask] = FPB_mean_sigma[gtx].fill_value
        FPB_median_corr[gtx].data[FPB_median_corr[gtx].mask] = FPB_median_corr[gtx].fill_value
        FPB_median_sigma[gtx].data[FPB_median_sigma[gtx].mask] = FPB_median_sigma[gtx].fill_value
        FPB_n_corr[gtx].data[FPB_n_corr[gtx].mask] = FPB_n_corr[gtx].fill_value
        FPB_cal_corr[gtx].data[FPB_cal_corr[gtx].mask] = FPB_cal_corr[gtx].fill_value
        TPS_mean_corr[gtx].data[TPS_mean_corr[gtx].mask] = TPS_mean_corr[gtx].fill_value
        TPS_median_corr[gtx].data[TPS_median_corr[gtx].mask] = TPS_median_corr[gtx].fill_value

        #-- save tep and dead time information and statistics
        IS2_atl03_fit['ancillary_data'][gtx] = {}
        IS2_atl03_attrs['ancillary_data'][gtx] = {}
        #-- tep time of day
        IS2_atl03_fit['ancillary_data'][gtx]['tep_tod'] = np.array(tep[gtx]['tep_tod'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod']['units'] = "seconds since 2018-01-01"
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod']['long_name'] = "TEP Time Of Day"
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod']['standard_name'] = "time"
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tep_tod']['description'] = ("The time of day "
            "at of the start of the data within the TEP histogram, in seconds since the "
            "ATLAS SDP GPS Epoch. The ATLAS Standard Data Products (SDP) epoch offset is "
            "defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. "
            "By adding the offset contained within atlas_sdp_gps_epoch to delta time "
            "parameters, the time in gps_seconds relative to the GPS epoch can be computed.")
        #-- tep window start
        IS2_atl03_fit['ancillary_data'][gtx]['tx_start'] = np.array(tep[gtx]['tx_start'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['long_name'] = "Start of the TEP Window"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['description'] = ("Starting time for the "
            "window centered around the primary TEP arrival for calculating the transmit pulse shape.")
        #-- tep window end
        IS2_atl03_fit['ancillary_data'][gtx]['tx_end'] = np.array(tep[gtx]['tx_end'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['long_name'] = "End of the TEP Window"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['description'] = ("Ending time for the "
            "window centered around the primary TEP arrival for calculating the transmit pulse shape.")
        #-- tep robust dispersion estimator
        IS2_atl03_fit['ancillary_data'][gtx]['tx_robust_sprd'] = np.array(tep[gtx]['tx_robust_sprd'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['long_name'] = "Robust Spread"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['description'] = ("Temporal width of "
            "the transmit pulse (sec), calculated from the RDE of the primary TEP waveform")
        #-- tep full width at half maximum
        IS2_atl03_fit['ancillary_data'][gtx]['sigma_tx'] = np.array(tep[gtx]['sigma_tx'])
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['long_name'] = "Duration of Transmit Pulse"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['description'] = ("Temporal duration of "
            "the transmit pulse (sec), calculated from the FWHM of the TEP waveform")
        #-- mean dead time
        IS2_atl03_fit['ancillary_data'][gtx]['t_dead'] = np.array(mean_dead_time[gtx])
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['long_name'] = "Dead-time"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['source'] = "CAL42"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['description'] = ("Mean dead-time for "
            "channels in the detector (sec)")

        #-- copy beam variables
        IS2_atl03_fit[gtx] = dict(land_ice_segments={})
        IS2_atl03_fill[gtx] = dict(land_ice_segments={})
        IS2_atl03_attrs[gtx] = dict(land_ice_segments={})
        #-- group attributes for beam
        IS2_atl03_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl03_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl03_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl03_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl03_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl03_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl03_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        #-- group attributes for land_ice_segments
        IS2_atl03_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        #-- geolocation, time and segment ID
        #-- delta time
        IS2_atl03_fit[gtx]['land_ice_segments']['delta_time'] = Segment_delta_time[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['delta_time'] = Segment_delta_time[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['units'] = "seconds since 2018-01-01"
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['long_name'] = "Elapsed GPS seconds"
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['standard_name'] = "time"
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['calendar'] = "standard"
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['description'] = ("Number of GPS "
            "seconds since the ATLAS SDP epoch. The ATLAS Standard Data Products (SDP) epoch offset "
            "is defined within /ancillary_data/atlas_sdp_gps_epoch as the number of GPS seconds "
            "between the GPS epoch (1980-01-06T00:00:00.000000Z UTC) and the ATLAS SDP epoch. By "
            "adding the offset contained within atlas_sdp_gps_epoch to delta time parameters, the "
            "time in gps_seconds relative to the GPS epoch can be computed.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['delta_time']['coordinates'] = \
            "segment_id latitude longitude"
        #-- latitude
        IS2_atl03_fit[gtx]['land_ice_segments']['latitude'] = Segment_Latitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['latitude'] = Segment_Latitude[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['units'] = "degrees_north"
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['long_name'] = "Latitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['standard_name'] = "latitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['description'] = ("Latitude of "
            "segment center")
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['valid_min'] = -90.0
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['valid_max'] = 90.0
        IS2_atl03_attrs[gtx]['land_ice_segments']['latitude']['coordinates'] = \
            "segment_id delta_time longitude"
        #-- longitude
        IS2_atl03_fit[gtx]['land_ice_segments']['longitude'] = Segment_Longitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['longitude'] = Segment_Longitude[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['units'] = "degrees_east"
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['long_name'] = "Longitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['standard_name'] = "longitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['description'] = ("Longitude of "
            "segment center")
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['valid_min'] = -180.0
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['valid_max'] = 180.0
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['coordinates'] = \
            "segment_id delta_time latitude"
        #-- segment ID
        IS2_atl03_fit[gtx]['land_ice_segments']['segment_id'] = Segment_ID[gtx][1:]
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']['long_name'] = "Along-track segment ID number"
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']['description'] = ("A 7 digit number "
            "identifying the along-track geolocation segment number.  These are sequential, starting with "
            "1 for the first segment after an ascending equatorial crossing node. Equal to the segment_id for "
            "the second of the two 20m ATL03 segments included in the 40m ATL06 segment")
        IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']['coordinates'] = \
            "delta_time latitude longitude"
        #-- land ice height corrected for first photon bias and transmit-pulse shape
        IS2_atl03_fit[gtx]['land_ice_segments']['h_li'] = Segment_Land_Ice[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['h_li'] = Segment_Land_Ice[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li']['long_name'] = "Land Ice height"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li']['description'] = ("Standard land-ice segment "
            "height determined by land ice algorithm, corrected for first-photon bias, representing the "
            "median-based height of the selected PEs")
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li']['coordinates'] = \
            "segment_id delta_time latitude longitude"
        #-- land ice height errors (max of fit or first photon bias uncertainties)
        IS2_atl03_fit[gtx]['land_ice_segments']['h_li_sigma'] = Segment_Land_Ice_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['h_li_sigma'] = Segment_Land_Ice_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma']['long_name'] = "Expected RMS segment misfit"
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma']['description'] = ("Propagated error due to "
            "sampling error and FPB correction from the land ice algorithm")
        IS2_atl03_attrs[gtx]['land_ice_segments']['h_li_sigma']['coordinates'] = \
            "segment_id delta_time latitude longitude"
        #-- vertical geolocation error due to PPD and POD
        IS2_atl03_fit[gtx]['land_ice_segments']['sigma_geo_h'] = Segment_sigma_geo[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['sigma_geo_h'] = Segment_sigma_geo[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h']['long_name'] = "Vertical Geolocation Error"
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h']['description'] = ("Total vertical geolocation error "
            "due to PPD and POD, including the effects of horizontal geolocation error on the segment vertical error.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['sigma_geo_h']['coordinates'] = \
            "segment_id delta_time latitude longitude"
        #-- segment quality summary
        IS2_atl03_fit[gtx]['land_ice_segments']['atl06_quality_summary'] = Segment_Summary[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['atl06_quality_summary'] = Segment_Summary[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['long_name'] = "ATL06 Quality Summary"
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['description'] = ("The ATL06_quality_summary "
            "parameter indicates the best-quality subset of all ATL06 data. A zero in this parameter implies that no "
            "data-quality tests have found a problem with the segment, a one implies that some potential problem has "
            "been found. Users who select only segments with zero values for this flag can be relatively certain of "
            "obtaining high-quality data, but will likely miss a significant fraction of usable data, particularly in "
            "cloudy, rough, or low-surface-reflectance conditions.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['flag_meanings'] = \
            "best_quality potential_problem"
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['valid_min'] = 0
        IS2_atl03_attrs[gtx]['land_ice_segments']['longitude']['valid_max'] = 1
        IS2_atl03_attrs[gtx]['land_ice_segments']['atl06_quality_summary']['coordinates'] = \
            "segment_id delta_time latitude longitude"

        #-- dem variables
        IS2_atl03_fit[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['Description'] = ("The dem group "
            "contains the reference digital elevation model and geoid heights.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        #-- geoid height
        fv = fileID[gtx]['geophys_corr']['geoid'].attrs['_FillValue']
        geoid = np.ma.array(fileID[gtx]['geophys_corr']['geoid'][:], fill_value=fv)
        geoid.mask = geoid.data == geoid.fill_value
        geoid_h = (geoid[1:] + geoid[0:-1])/2.0
        geoid_h.data[geoid_h.mask] = geoid_h.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['dem']['geoid_h'] = geoid_h
        IS2_atl03_fill[gtx]['land_ice_segments']['dem']['geoid_h'] = geoid_h.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['long_name'] = "Geoid Height"
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['description'] = ("Geoid height above "
            "WGS-84 reference ellipsoid (range -107 to 86m)")
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['source'] = "EGM2008"
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['valid_min'] = -107
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['valid_max'] = 86
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['geoid_h']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        #-- geophysical variables
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['Description'] = ("The geophysical group "
            "contains parameters used to correct segment heights for geophysical effects, parameters "
            "related to solar background and parameters indicative of the presence or absence of clouds.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        #-- background rate
        bckgrd = (Segment_Background[gtx][1:] + Segment_Background[gtx][0:-1])/2.0
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['bckgrd'] = np.copy(bckgrd)
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['bckgrd'] = None
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd']['units'] = "counts / second"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd']['long_name'] = ("Background count "
            "rate based on the ATLAS 50-shot sum interpolated to the reference photon")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd']['description'] = ("The background "
            "count rate from the 50-shot altimetric histogram after removing the number of likely signal photons")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bckgrd']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- blowing snow PSC flag
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['bsnow_psc'] = \
            IS2_atl09_mds[pfl]['high_rate']['bsnow_psc'][1:]
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['bsnow_psc'] = None
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['long_name'] = "Blowing snow PSC flag"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['description'] = ("Indicates the "
            "potential for polar stratospheric clouds to affect the blowing snow retrieval, where 0=none and 3=maximum. "
            "This flag is a function of month and hemisphere and is only applied poleward of 60 north and south")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['flag_meanings'] = ("none slight "
            "moderate maximum_bsnow_PSC_affected")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['flag_values'] = [0,1,2,3]
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['valid_min'] = 0
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['valid_max'] = 3
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_psc']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- blowing snow confidence
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['bsnow_conf'] = \
            IS2_atl09_mds[pfl]['high_rate']['bsnow_con'][1:]
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['bsnow_conf'] = None
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['long_name'] = "Blowing snow confidence"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['description'] = ("Indicates the blowing snow "
            "confidence, where -3=surface not detected; 2=no surface wind;-1=no scattering layer found; 0=no top layer found; "
            "1=none-little; 2=weak; 3=moderate; 4=moderate-high; 5=high; 6=very high")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['flag_meanings'] = ("surface_not_detected "
            "no_surface_wind no_scattering_layer_found no_top_layer_found none_little weak moderate moderate_high high very_high")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['flag_values'] = [-3,-2,-1,0,1,2,3,4,5,6]
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['valid_min'] = -3
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['valid_max'] = 6
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_conf']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- blowing snow optical depth
        bsnow_od = np.ma.array(IS2_atl09_mds[pfl]['high_rate']['bsnow_od'][1:],
            fill_value=IS2_atl09_attrs[pfl]['high_rate']['bsnow_od']['_FillValue'])
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['bsnow_od'] = bsnow_od
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['bsnow_od'] = bsnow_od.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od']['long_name'] = "Blowing snow OD"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od']['description'] = "Blowing snow layer optical depth"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['bsnow_od']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- cloud flag ASR
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr'] = \
            IS2_atl09_mds[pfl]['high_rate']['cloud_flag_asr'][1:]
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr'] = None
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['long_name'] = "Cloud Flag ASR"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['description'] = ("Indicates the Cloud "
            "flag probability from apparent surface reflectance, where 0=clear with high confidence; 1=clear with medium "
            "confidence; 2=clear with low confidence; 3=cloudy with low confidence; 4=cloudy with medium confidence; 5=cloudy "
            "with high confidence; 6=unknown")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['flag_meanings'] = ("clear_with_high_confidence "
            "clear_with_medium_confidence clear_with_low_confidence cloudy_with_low_confidence cloudy_with_medium_confidence "
            "cloudy_with_high_confidence unknown")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['flag_values'] = [0,1,2,3,4,5,6]
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['valid_min'] = 0
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['valid_max'] = 6
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_asr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- cloud flag atm
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm'] = \
            IS2_atl09_mds[pfl]['high_rate']['cloud_flag_atm'][1:]
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm'] = None
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['long_name'] = "Cloud Flag Atm"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['description'] = ("Number of layers found "
            "from the backscatter profile using the DDA layer finder")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['flag_values'] = [0,1,2,3,4,5,6,7,8,9,10]
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['valid_min'] = 0
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['valid_max'] = 10
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['cloud_flag_atm']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- multiple scattering warning flag
        msw_flag = np.ma.array(IS2_atl09_mds[pfl]['high_rate']['msw_flag'][1:],
            fill_value=IS2_atl09_attrs[pfl]['high_rate']['msw_flag']['_FillValue'])
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['msw_flag'] = msw_flag
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['msw_flag'] = msw_flag.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['long_name'] = "Multiple Scattering Warning Flag"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['description'] = ("Combined flag indicating the "
            "risks of severe multiple scattering. The multiple scattering warning flag (ATL09 parameter msw_flag) has values from "
            "-1 to 5 where zero means no multiple scattering and 5 the greatest. If no layers were detected, then msw_flag = 0. "
            "If blowing snow is detected and its estimated optical depth is greater than or equal to 0.5, then msw_flag = 5. "
            "If the blowing snow optical depth is less than 0.5, then msw_flag = 4. If no blowing snow is detected but there are "
            "cloud or aerosol layers detected, the msw_flag assumes values of 1 to 3 based on the height of the bottom of the "
            "lowest layer: < 1 km, msw_flag = 3; 1-3 km, msw_flag = 2; > 3km, msw_flag = 1. A value of -1 indicates that the "
            "signal to noise of the data was too low to reliably ascertain the presence of cloud or blowing snow. We expect "
            "values of -1 to occur only during daylight.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['flag_meanings'] = ("cannot_determine "
            "no_layers layer_gt_3km layer_between_1_and_3_km layer_lt_1km blow_snow_od_lt_0.5 blow_snow_od_gt_0.5")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['flag_values'] = [-1,0,1,2,3,4,5]
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['valid_min'] = -1
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['valid_max'] = 5
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['msw_flag']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- range bias correction
        fv = fileID[gtx]['geolocation']['range_bias_corr'].attrs['_FillValue']
        range_bias_corr = np.ma.array(fileID[gtx]['geolocation']['range_bias_corr'][:], fill_value=fv)
        range_bias_corr.mask = range_bias_corr.data == range_bias_corr.fill_value
        segment_range_bias = (range_bias_corr[1:] + range_bias_corr[0:-1])/2.0
        segment_range_bias.data[segment_range_bias.mask] = segment_range_bias.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = segment_range_bias
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = segment_range_bias.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['long_name'] = "Range bias correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['description'] = "The range_bias estimated from geolocation analysis"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- total neutral atmosphere delay correction
        fv = fileID[gtx]['geolocation']['neutat_delay_total'].attrs['_FillValue']
        neutat_delay_total = np.ma.array(fileID[gtx]['geolocation']['neutat_delay_total'][:], fill_value=fv)
        neutat_delay_total.mask = neutat_delay_total.data == neutat_delay_total.fill_value
        segment_neutat_delay = (neutat_delay_total[1:] + neutat_delay_total[0:-1])/2.0
        segment_neutat_delay.data[segment_neutat_delay.mask] = segment_neutat_delay.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = segment_neutat_delay
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = segment_neutat_delay.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['long_name'] = "Total Neutral Atmospheric Delay"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['description'] = "Total neutral atmosphere delay correction (wet+dry)"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- solar elevation
        fv = fileID[gtx]['geolocation']['solar_elevation'].attrs['_FillValue']
        solar_elevation = np.ma.array(fileID[gtx]['geolocation']['solar_elevation'][:], fill_value=fv)
        solar_elevation.mask = solar_elevation.data == solar_elevation.fill_value
        segment_solar_elevation = (solar_elevation[1:] + solar_elevation[0:-1])/2.0
        segment_solar_elevation.data[segment_solar_elevation.mask] = segment_solar_elevation.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['solar_elevation'] = segment_solar_elevation
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['solar_elevation'] = segment_solar_elevation.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation']['units'] = "degrees"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation']['long_name'] = "Solar elevation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation']['description'] = ("Solar Angle above "
            "or below the plane tangent to the ellipsoid surface at the laser spot. Positive values mean the sun is above the "
            "horizon, while negative values mean it is below the horizon. The effect of atmospheric refraction is not included.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_elevation']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- solar azimuth
        fv = fileID[gtx]['geolocation']['solar_azimuth'].attrs['_FillValue']
        solar_azimuth = np.ma.array(fileID[gtx]['geolocation']['solar_azimuth'][:], fill_value=fv)
        solar_azimuth.mask = solar_azimuth.data == solar_azimuth.fill_value
        segment_solar_azimuth = (solar_azimuth[1:] + solar_azimuth[0:-1])/2.0
        segment_solar_azimuth.data[segment_solar_azimuth.mask] = segment_solar_azimuth.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['solar_azimuth'] = segment_solar_azimuth
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['solar_azimuth'] = segment_solar_azimuth.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth']['units'] = "degrees_east"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth']['long_name'] = "Solar azimuth"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth']['description'] = ("The direction, "
            "eastwards from north, of the sun vector as seen by an observer at the laser ground spot.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['solar_azimuth']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        #-- geophysical correction values at segment reference photons
        #-- dynamic atmospheric correction
        fv = fileID[gtx]['geophys_corr']['dac'].attrs['_FillValue']
        dac = np.ma.array(fileID[gtx]['geophys_corr']['dac'][:], fill_value=fv)
        dac.mask = dac.data == dac.fill_value
        segment_dac = (dac[1:] + dac[0:-1])/2.0
        segment_dac.data[segment_dac.mask] = segment_dac.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['dac'] = segment_dac
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['dac'] = segment_dac.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['long_name'] = "Dynamic Atmosphere Correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['description'] = ("Dynamic Atmospheric Correction "
            "(DAC) includes inverted barometer (IB) effect")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['source'] = 'Mog2D-G'
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['dac']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- solid earth tide
        fv = fileID[gtx]['geophys_corr']['tide_earth'].attrs['_FillValue']
        tide_earth = np.ma.array(fileID[gtx]['geophys_corr']['tide_earth'][:], fill_value=fv)
        tide_earth.mask = tide_earth.data == tide_earth.mask
        segment_earth_tide = (tide_earth[1:] + tide_earth[0:-1])/2.0
        segment_earth_tide.data[segment_earth_tide.mask] = segment_earth_tide.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_earth'] = segment_earth_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_earth'] = segment_earth_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['long_name'] = "Earth Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['description'] = "Solid Earth Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- load tide
        fv = fileID[gtx]['geophys_corr']['tide_load'].attrs['_FillValue']
        tide_load = np.ma.array(fileID[gtx]['geophys_corr']['tide_load'][:], fill_value=fv)
        tide_load.mask = tide_load.data == tide_load.fill_value
        segment_load_tide = (tide_load[1:] + tide_load[0:-1])/2.0
        segment_load_tide.data[segment_load_tide.mask] = segment_load_tide.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_load'] = segment_load_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_load'] = segment_load_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load']['long_name'] = "Load Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load']['description'] = ("Load Tide - Local "
            "displacement due to Ocean Loading (-6 to 0 cm)")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_load']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- ocean tide
        fv = fileID[gtx]['geophys_corr']['tide_ocean'].attrs['_FillValue']
        tide_ocean = np.ma.array(fileID[gtx]['geophys_corr']['tide_ocean'][:], fill_value=fv)
        tide_ocean.mask = tide_ocean.data == tide_ocean.fill_value
        segment_ocean_tide = (tide_ocean[1:] + tide_ocean[0:-1])/2.0
        segment_ocean_tide.data[segment_ocean_tide.mask] = segment_ocean_tide.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_ocean'] = segment_ocean_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_ocean'] = segment_ocean_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean']['long_name'] = "Ocean Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean']['description'] = ("Ocean Tides "
            "including diurnal and semi-diurnal (harmonic analysis), and longer period tides (dynamic and "
            "self-consistent equilibrium).")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_ocean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- ocean pole tide
        fv = fileID[gtx]['geophys_corr']['tide_oc_pole'].attrs['_FillValue']
        tide_oc_pole = np.ma.array(fileID[gtx]['geophys_corr']['tide_oc_pole'][:], fill_value=fv)
        tide_oc_pole.mask = tide_oc_pole.data == tide_oc_pole.fill_value
        segment_oc_pole_tide = (tide_oc_pole[1:] + tide_oc_pole[0:-1])/2.0
        segment_oc_pole_tide.data[segment_oc_pole_tide.mask] = segment_oc_pole_tide.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_oc_pole'] = segment_oc_pole_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_oc_pole'] = segment_oc_pole_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole']['long_name'] = "Ocean Pole Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole']['description'] = ("Oceanic surface "
            "rotational deformation due to polar motion (-2 to 2 mm).")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_oc_pole']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- pole tide
        fv = fileID[gtx]['geophys_corr']['tide_pole'].attrs['_FillValue']
        tide_pole = np.ma.array(fileID[gtx]['geophys_corr']['tide_pole'][:], fill_value=fv)
        tide_pole.mask = tide_pole.data == tide_pole.fill_value
        segment_pole_tide = (tide_pole[1:] + tide_pole[0:-1])/2.0
        segment_pole_tide.data[segment_pole_tide.mask] = segment_pole_tide.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_pole'] = segment_pole_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_pole'] = segment_pole_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole']['long_name'] = "Solid Earth Pole Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole']['description'] = ("Solid Earth Pole "
            "Tide - Rotational deformation due to polar motion  (-1.5 to 1.5 cm).")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_pole']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        #-- bias correction variables
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['Description'] = ("The bias_correction group "
            "contains information about the estimated first-photon bias, and the transmit-pulse-shape bias.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        #-- mean first photon bias
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr'] = FPB_mean_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr'] = FPB_mean_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr']['long_name'] = "first photon bias mean correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr']['description'] = ("Estimated first-photon-bias "
            "(fpb) correction to mean segment height")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- mean first photon bias uncertainty
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma'] = FPB_mean_sigma[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma'] = FPB_mean_sigma[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma']['long_name'] = ("first photon bias "
            "mean correction error")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma']['description'] = ("Estimated error in "
            "first-photon-bias (fpb) correction for mean segment heights")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_mean_corr_sigma']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- median first photon bias
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr'] = FPB_median_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr'] = FPB_median_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr']['long_name'] = "first photon bias median correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr']['description'] = ("Estimated first-photon-bias "
            "(fpb) correction giving the difference between the mean segment height and the corrected median height")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- median first photon bias correction
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma'] = FPB_median_sigma[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma'] = FPB_median_sigma[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma']['long_name'] = ("first photon bias median "
            "correction error")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma']['description'] = ("Estimated error in "
            "first-photon-bias (fpb) correction giving the difference between the mean segment height and the corrected median height")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_med_corr_sigma']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- first photon bias corrected number of photons
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr'] = FPB_n_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr'] = FPB_n_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr']['long_name'] = "corrected number of photons"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr']['description'] = ("Estimated photon count after "
            "first-photon-bias correction")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_n_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- CAL-19 first photon bias
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr'] = FPB_cal_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr'] = FPB_cal_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr']['long_name'] = "first photon bias calibrated correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr']['description'] = ("Estimated first-photon-bias "
            "(fpb) correction calculated using the ATL03 calibration products")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['fpb_cal_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- mean transmit pulse shape correction
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr'] = TPS_mean_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr'] = TPS_mean_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['long_name'] = "tx shape mean correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['description'] = ("Estimate of the difference "
            "between the mean of the full-waveform transmit-pulse and the mean of a broadened, truncated waveform consistent "
            "with the received pulse")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_mean_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- median transmit pulse shape correction
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['tx_med_corr'] = TPS_median_corr[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['tx_med_corr'] = TPS_median_corr[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['long_name'] = "tx shape median correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['description'] = ("Estimate of the difference "
            "between the median of the full-waveform transmit-pulse and the median of a broadened, truncated waveform consistent "
            "with the received pulse")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['tx_med_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- difference between the mean and median of the histogram fit residuals
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = Segment_Mean_Median[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = Segment_Mean_Median[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['long_name'] = "mean median residual"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['description'] = ("Difference between "
            "uncorrected mean and median of the histogram fit residuals")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        #-- fit statistics variables
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['Description'] = ("The fit_statistics group "
            "contains a variety of parameters that might indicate the quality of the fitted segment data. Data in "
            "this group are sparse, with dimensions matching the land_ice_segments group.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        #-- segment histogram fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = Segment_Height[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = Segment_Height[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['long_name'] = "Height Mean"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['description'] = ("Mean surface "
            "height from histogram decomposition, not corrected for first-photon bias or pulse truncation")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment minimum histogram fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_min'] = Segment_Minimum[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_min'] = Segment_Minimum[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min']['long_name'] = "Minimum Height"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min']['description'] = ("Minimum surface "
            "height from histogram decomposition, not corrected for first-photon bias or pulse truncation")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_min']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment maximum histogram fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_max'] = Segment_Maximum[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_max'] = Segment_Maximum[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max']['long_name'] = "Maximum Height"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max']['description'] = ("Maximum surface "
            "height from histogram decomposition, not corrected for first-photon bias or pulse truncation")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_max']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment fit along-track slopes
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx'] = Segment_dH_along[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx'] = Segment_dH_along[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx']['units'] = "meters/meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx']['long_name'] = "Along Track Slope"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx']['description'] = ("Along-track slope "
            "from along-track segment fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment fit across-track slopes
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy'] = Segment_dH_across[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy'] = Segment_dH_across[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy']['units'] = "meters/meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy']['long_name'] = "Across Track Slope"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy']['description'] = ("Across track slope "
            "from segment fits to weak and strong beam; the same slope is reported for both laser beams in each pair")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment fit height errors
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = Segment_Height_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = Segment_Height_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['long_name'] = "Height Error"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['description'] = ("Propagated height "
            "error due to PE-height sampling error for the height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment histogram fit errors (minimum)
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min'] = Segment_Minimum_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min'] = Segment_Minimum_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min']['long_name'] = "Minimum Height Error"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min']['description'] = ("Propagated height "
            "error due to PE-height sampling error for the minimum height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_min']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment histogram fit errors (maximum)
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max'] = Segment_Maximum_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max'] = Segment_Maximum_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max']['long_name'] = "Maximum Height Error"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max']['description'] = ("Propagated height "
            "error due to PE-height sampling error for the maximum height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_max']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment fit across-track slope errors
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma'] = Segment_dH_along_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma'] = Segment_dH_along_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma']['units'] = "meters/meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma']['long_name'] = "Sigma of Along Track Slope"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma']['description'] = ("Propagated error in "
            "the along-track segment slope")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dx_sigma']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment fit along-track slope errors
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma'] = Segment_dH_across_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma'] = Segment_dH_across_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma']['units'] = "meters/meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma']['long_name'] = "Sigma of Across Track Slope"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma']['description'] = ("Propagated error in "
            "the across-track segment slope calculated from segment fits to weak and strong beam")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['dh_fit_dy_sigma']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment histogram fit amplitudes
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['amp_mean'] = Segment_Amplitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['amp_mean'] = Segment_Amplitude[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean']['long_name'] = "Amplitude Mean"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean']['description'] = ("Amplitude of the "
            "mean surface height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_mean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment minimum histogram fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['amp_min'] = Segment_Minimum_Amplitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['amp_min'] = Segment_Minimum_Amplitude[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min']['long_name'] = "Minimum Amplitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min']['description'] = ("Amplitude of the "
            "minimum surface height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_min']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- segment maximum histogram fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['amp_max'] = Segment_Maximum_Amplitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['amp_max'] = Segment_Maximum_Amplitude[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max']['long_name'] = "Maximum Amplitude"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max']['description'] = ("Amplitude of the  "
            "maximum surface height from histogram decomposition")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['amp_max']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- number of photons in fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = Segment_N_Fit[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = Segment_N_Fit[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['long_name'] = "Number of Photons in Fit"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['description'] = ("Number of PEs used to "
            "determine mean surface height in the iterative histogram fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- number of peaks in histogram fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks'] = Segment_N_Peaks[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks'] = Segment_N_Peaks[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks']['long_name'] = "Number of Photons in Fit"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks']['description'] = ("Number of Peaks used to "
            "in the iterative histogram fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_peaks']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- size of the window used in the fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = Segment_Window[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = Segment_Window[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['long_name'] = "Surface Window Width"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['description'] = ("Width of the surface "
            "window from the histogram fit, top to bottom")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- signal-to-noise ratio
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['snr'] = Segment_SNR[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['snr'] = Segment_SNR[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr']['long_name'] = "SNR"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr']['description'] = ("Signal-to-noise "
            "ratio in the final refined window")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- robust dispersion estimator
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = Segment_RDE[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = Segment_RDE[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['long_name'] = "Robust Spread"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['description'] = ("RDE of misfit "
            "between PE heights and the histogram fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- number of iterations for fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations'] = Segment_Iterations[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations'] = Segment_Iterations[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations']['long_name'] = "Number of Iterations used in Fit"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations']['description'] = ("Number of Iterations when "
            "determining the mean surface height")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_iterations']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- signal source selection
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source'] = Segment_Source[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source'] = Segment_Source[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['long_name'] = "Signal Selection Source"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['description'] = ("Indicates the last "
            "algorithm attempted to select the signal for fitting. 1=Signal selection succeeded using ATL03 detected PE; 2=Signal "
            "selection failed using ATL03 detected PE but succeeded using all flagged ATL03 PE; 3=Signal selection failed using "
            "all flagged ATL03 PE, but succeeded using the backup algorithm; 4=All signal-finding strategies failed.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['flag_values'] = [1,2,3,4]
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['valid_min'] = 1
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['valid_max'] = 4
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['signal_selection_source']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- number potential segment pulses
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses'] = Segment_Pulses[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses'] = Segment_Pulses[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses']['long_name'] = "Number potential segment pulses"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses']['description'] = ("The number of pulses "
            "potentially included in the segment")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_seg_pulses']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        #-- ground track variables
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['Description'] = ("The ground_track group "
            "contains parameters describing the GT and RGT for each land ice segment, as well as angular "
            "information about the beams.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        #-- along-track X coordinates of segment fit
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track']['x_atc'] = Segment_X_atc[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track']['x_atc'] = Segment_X_atc[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc']['long_name'] = "X Along Track"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc']['description'] = ("The along-track "
            "x-coordinate of the segment, measured parallel to the RGT, measured from the ascending node of the equatorial "
            "crossing of a given RGT.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_atc']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- along-track Y coordinates of segment fit
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track']['y_atc'] = Segment_Y_atc[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track']['y_atc'] = Segment_Y_atc[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc']['long_name'] = "Y Along Track"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc']['description'] = ("The along-track "
            "y-coordinate of the segment, relative to the RGT, measured along the perpendicular to the RGT, "
            "positive to the right of the RGT.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['y_atc']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- along-track X coordinate spread of points used in segment fit
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track']['x_spread'] = Segment_X_spread[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track']['x_spread'] = Segment_X_spread[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread']['long_name'] = "X Along Track Spread"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread']['description'] = ("The spread of "
            "along-track x-coordinates for points used in the segment fit.  Coordinates measured parallel to the "
            "RGT, measured from the ascending node of the equatorial crossing of a given RGT.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['x_spread']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- elevation
        fv = fileID[gtx]['geolocation']['ref_elev'].attrs['_FillValue']
        ref_elev = np.ma.array(fileID[gtx]['geolocation']['ref_elev'][:], fill_value=fv)
        ref_elev.mask = ref_elev.data == ref_elev.fill_value
        segment_ref_elev = (ref_elev[1:] + ref_elev[0:-1])/2.0
        segment_ref_elev.data[segment_ref_elev.mask] = segment_ref_elev.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track']['ref_elev'] = segment_ref_elev
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track']['ref_elev'] = segment_ref_elev.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev']['units'] = "radians"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev']['long_name'] = "Elevation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev']['description'] = ("Elevation of the "
            "unit pointing vector for the reference photon in the local ENU frame in radians.  The angle is measured "
            "from East-North plane and positive towards Up")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_elev']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        #-- azimuth
        fv = fileID[gtx]['geolocation']['ref_azimuth'].attrs['_FillValue']
        ref_azimuth = np.ma.array(fileID[gtx]['geolocation']['ref_azimuth'][:], fill_value=fv)
        ref_azimuth.mask = ref_azimuth.data == ref_azimuth.fill_value
        segment_ref_azimuth = (ref_azimuth[1:] + ref_azimuth[0:-1])/2.0
        segment_ref_azimuth.data[segment_ref_azimuth.mask] = segment_ref_azimuth.fill_value
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track']['ref_azimuth'] = segment_ref_azimuth
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track']['ref_azimuth'] = segment_ref_azimuth.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth']['units'] = "radians"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth']['long_name'] = "Azimuth"
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth']['description'] = ("Azimuth of the "
            "unit pointing vector for the reference photon in the local ENU frame in radians.  The angle is measured "
            "from North and positive towards East")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['ref_azimuth']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

    #-- parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0):
        #-- use default output file name and path
        if args.output:
            output_file=os.path.expanduser(args.output)
        else:
            args=(SUB,'ATL86',YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX)
            file_format='{0}{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
            output_file=os.path.join(ATL03_dir,file_format.format(*args))
        #-- write to HDF5 file
        HDF5_ATL03_write(IS2_atl03_fit, IS2_atl03_attrs, COMM=comm,
            VERBOSE=args.verbose, INPUT=[args.ATL03,args.ATL09],
            FILL_VALUE=IS2_atl03_fill, CLOBBER=True, FILENAME=output_file)
        #-- change the permissions level to MODE
        os.chmod(output_file, args.mode)
    #-- close the input ATL03 file
    fileID.close()

#-- PURPOSE: read ICESat-2 ATL09 HDF5 data file for specific variables
def read_HDF5_ATL09(FILENAME, pfl, D, ATTRIBUTES=True, VERBOSE=False, COMM=None):
    #-- Open the HDF5 file for reading
    fileID = h5py.File(FILENAME, 'r', driver='mpio', comm=COMM)
    print(FILENAME) if VERBOSE and (COMM.rank == 0) else None

    #-- allocate python dictionaries for ICESat-2 ATL09 variables and attributes
    IS2_atl09_mds = {}
    IS2_atl09_attrs = {}

    #-- read profile reported for the ATLAS strong beams within the file
    IS2_atl09_mds[pfl] = dict(high_rate={})
    #-- extract delta_time for mapping ATL09 atmospheric parameters to ATL03
    delta_time = fileID[pfl]['high_rate']['delta_time'][:]
    #-- Calibrated Attenuated Backscatter at 25 hz
    high_rate_keys = ['aclr_true','bsnow_con','bsnow_dens','bsnow_h',
        'bsnow_h_dens','bsnow_od','bsnow_psc','cloud_flag_asr','cloud_flag_atm',
        'cloud_fold_flag','column_od_asr','column_od_asr_qf','msw_flag',
        'snow_ice','solar_azimuth','solar_elevation','surf_refl_true']
    #-- number of output ATL03 segments
    n_seg = len(D)
    #-- parallel indices for filling variables
    ii = np.arange(COMM.rank,n_seg,COMM.size)
    #-- extract variables of interest and map to ATL03 segments
    for key in high_rate_keys:
        val = np.copy(fileID[pfl]['high_rate'][key][:])
        fint = scipy.interpolate.interp1d(delta_time, val,
            kind='nearest', fill_value='extrapolate')
        IS2_atl09_mds[pfl]['high_rate'][key] = np.zeros((n_seg),dtype=val.dtype)
        IS2_atl09_mds[pfl]['high_rate'][key][ii] = fint(D[ii]).astype(val.dtype)

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- Getting attributes of IS2_atl09_mds profile variables
        IS2_atl09_attrs[pfl] = dict(high_rate={})
        #-- Global Group Attributes
        for att_name,att_val in fileID[pfl].attrs.items():
            IS2_atl09_attrs[pfl][att_name] = att_val
        #-- Variable Attributes
        for key in high_rate_keys:
            IS2_atl09_attrs[pfl]['high_rate'][key] = {}
            for att_name,att_val in fileID[pfl]['high_rate'][key].attrs.items():
                IS2_atl09_attrs[pfl]['high_rate'][key][att_name] = att_val

    #-- Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl09_attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl09_mds,IS2_atl09_attrs)

#-- PURPOSE: outputting the reduced and corrected ICESat-2 data to HDF5
def HDF5_ATL03_write(IS2_atl03_data, IS2_atl03_attrs, COMM=None, INPUT=None,
    FILENAME='', FILL_VALUE=None, CLOBBER=True, VERBOSE=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- open output HDF5 file
    fileID = h5py.File(FILENAME, clobber)#, driver='mpio', comm=COMM)
    print(FILENAME) if VERBOSE and (COMM.rank == 0) else None

    #-- create HDF5 records
    h5 = {}

    # #-- ICESat-2 spacecraft orientation at time
    # fileID.create_group('orbit_info')
    # h5['orbit_info'] = {}
    # for k,v in IS2_atl03_data['orbit_info'].items():
    #     #-- Defining the HDF5 dataset variables
    #     val = 'orbit_info/{0}'.format(k)
    #     h5['orbit_info'][k] = fileID.create_dataset(val, np.shape(v), data=v,
    #         dtype=v.dtype, compression='gzip')
    #     #-- add HDF5 variable attributes
    #     for att_name,att_val in IS2_atl03_attrs['orbit_info'][k].items():
    #         h5['orbit_info'][k].attrs[att_name] = att_val

    #-- information ancillary to the data product
    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k in ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc','end_cycle',
        'end_geoseg','end_gpssow','end_gpsweek','end_orbit','end_region',
        'end_rgt','granule_end_utc','granule_start_utc','release','start_cycle',
        'start_geoseg','start_gpssow','start_gpsweek','start_orbit','start_region',
        'start_rgt','version']:
        #-- Defining the HDF5 dataset variables
        v = IS2_atl03_data['ancillary_data'][k]
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype)
        #-- add HDF5 variable attributes
        for att_name,att_val in IS2_atl03_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    #-- land_ice_segments variable groups for each beam
    GROUPS=['fit_statistics','geophysical','ground_track','dem','bias_correction']
    #-- write each output beam
    for gtx in ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']:
        fileID.create_group(gtx)
        fileID['ancillary_data'].create_group(gtx)

        #-- add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl03_attrs[gtx][att_name]

        #-- add transmit pulse shape and dead time parameters
        h5['ancillary_data'][gtx] = {}
        for k,v in IS2_atl03_data['ancillary_data'][gtx].items():
            #-- attributes
            attrs = IS2_atl03_attrs['ancillary_data'][gtx][k]
            #-- Defining the HDF5 dataset variables
            val = 'ancillary_data/{0}/{1}'.format(gtx,k)
            h5['ancillary_data'][gtx][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype)
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5['ancillary_data'][gtx][k].attrs[att_name] = att_val

        #-- create land_ice_segments group
        fileID[gtx].create_group('land_ice_segments')
        h5[gtx] = dict(land_ice_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl03_attrs[gtx]['land_ice_segments'][att_name]
            fileID[gtx]['land_ice_segments'].attrs[att_name] = att_val

        #-- segment_id
        v = IS2_atl03_data[gtx]['land_ice_segments']['segment_id']
        attrs = IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']
        # #-- parallel indices for filling variables
        # pind = np.arange(COMM.rank,len(v),COMM.size)
        #-- Defining the HDF5 dataset variables
        val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments','segment_id')
        h5[gtx]['land_ice_segments']['segment_id'] = fileID.create_dataset(val,
            np.shape(v), data=v, dtype=v.dtype, compression='gzip')
        # with h5[gtx]['land_ice_segments']['segment_ID'].collective:
        #     h5[gtx]['land_ice_segments']['segment_id'][pind] = v[pind]
        #-- add HDF5 variable attributes
        for att_name,att_val in attrs.items():
            h5[gtx]['land_ice_segments']['segment_id'].attrs[att_name] = att_val

        #-- geolocation, time and height variables
        for k in ['latitude','longitude','delta_time','h_li','h_li_sigma',
            'sigma_geo_h','atl06_quality_summary']:
            #-- values and attributes
            v = IS2_atl03_data[gtx]['land_ice_segments'][k]
            attrs = IS2_atl03_attrs[gtx]['land_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][k]
            #-- Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments',k)
            h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                compression='gzip')
            # with h5[gtx]['land_ice_segments'][k].collective:
            #     h5[gtx]['land_ice_segments'][k][pind] = v[pind]
            #-- attach dimensions
            for dim in ['segment_id']:
                h5[gtx]['land_ice_segments'][k].dims.create_scale(
                    h5[gtx]['land_ice_segments'][dim], dim)
                h5[gtx]['land_ice_segments'][k].dims[0].attach_scale(
                    h5[gtx]['land_ice_segments'][dim])
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][k].attrs[att_name] = att_val

        #-- fit statistics, geophysical corrections, geolocation and dem
        for key in GROUPS:
            fileID[gtx]['land_ice_segments'].create_group(key)
            h5[gtx]['land_ice_segments'][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl03_attrs[gtx]['land_ice_segments'][key][att_name]
                fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
            for k,v in IS2_atl03_data[gtx]['land_ice_segments'][key].items():
                #-- attributes
                attrs = IS2_atl03_attrs[gtx]['land_ice_segments'][key][k]
                fillvalue = FILL_VALUE[gtx]['land_ice_segments'][key][k]
                #-- Defining the HDF5 dataset variables
                val = '{0}/{1}/{2}/{3}'.format(gtx,'land_ice_segments',key,k)
                if fillvalue:
                    h5[gtx]['land_ice_segments'][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                            dtype=v.dtype, fillvalue=fillvalue, compression='gzip')
                else:
                    h5[gtx]['land_ice_segments'][key][k] = \
                        fileID.create_dataset(val, np.shape(v), data=v,
                            dtype=v.dtype, compression='gzip')
                # with h5[gtx]['land_ice_segments'][key][k].collective:
                #     h5[gtx]['land_ice_segments'][key][k][pind] = v[pind]
                #-- attach dimensions
                for dim in ['segment_id']:
                    h5[gtx]['land_ice_segments'][key][k].dims.create_scale(
                        h5[gtx]['land_ice_segments'][dim], dim)
                    h5[gtx]['land_ice_segments'][key][k].dims[0].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
                #-- add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx]['land_ice_segments'][key][k].attrs[att_name] = att_val

    #-- HDF5 file title
    fileID.attrs['featureType'] = 'trajectory'
    fileID.attrs['title'] = 'ATLAS/ICESat-2 Land Ice Height'
    fileID.attrs['summary'] = ('Estimates of the ice-sheet mean surface height '
        'relative to the WGS-84 ellipsoid, and ancillary parameters needed to '
        'interpret and assess the quality of these height estimates.')
    fileID.attrs['description'] = ('Land ice surface heights for each beam, '
        'along and across-track slopes calculated for beam pairs.  All '
        'parameters are calculated for the same along-track increments for '
        'each beam and repeat.')
    date_created = datetime.datetime.today()
    fileID.attrs['date_created'] = date_created.isoformat()
    project = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = project
    platform = 'ICESat-2 > Ice, Cloud, and land Elevation Satellite-2'
    fileID.attrs['project'] = platform
    #-- add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'http://nsidc.org/data/icesat2/data.html'
    fileID.attrs['processing_level'] = '4'
    #-- add attributes for input ATL03 and ATL09 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    #-- find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']:
        lon = IS2_atl03_data[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl03_data[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl03_data[gtx]['land_ice_segments']['delta_time']
        #-- setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time.min() if (delta_time.min() < tmn) else tmn
        tmx = delta_time.max() if (delta_time.max() > tmx) else tmx
    #-- add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    #-- convert start and end time from ATLAS SDP seconds into UTC time
    time_utc = convert_delta_time(np.array([tmn,tmx]))
    #-- convert to calendar date with convert_julian.py
    YY,MM,DD,HH,MN,SS = convert_julian(time_utc['julian'],FORMAT='tuple')
    #-- add attributes with measurement date start, end and duration
    tcs = datetime.datetime(np.int(YY[0]), np.int(MM[0]), np.int(DD[0]),
        np.int(HH[0]), np.int(MN[0]), np.int(SS[0]), np.int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(np.int(YY[1]), np.int(MM[1]), np.int(DD[1]),
        np.int(HH[1]), np.int(MN[1]), np.int(SS[1]), np.int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = '{0:0.0f}'.format(tmx-tmn)
    #-- Closing the HDF5 file
    fileID.close()

#-- run main program
if __name__ == '__main__':
    main()
