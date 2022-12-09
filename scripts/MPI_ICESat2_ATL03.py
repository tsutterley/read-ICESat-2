#!/usr/bin/env python
u"""
MPI_ICESat2_ATL03.py (12/2022)
Read ICESat-2 ATL03 and ATL09 data files to calculate average segment surfaces
    ATL03 datasets: Global Geolocated Photons
    ATL09 datasets: Atmospheric Characteristics

CALLING SEQUENCE:
    mpiexec -np 6 python MPI_ICESat2_ATL03.py ATL03_file ATL09_file

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
    convert_delta_time.py: converts from delta time into Julian and year-decimal
    fit.py: Utilities for calculating fits from ATL03 Geolocated Photon Data
    time.py: Utilities for calculating time operations
    utilities.py: download and management utilities for syncing files
    classify_photons.py: Yet Another Photon Classifier for Geolocated Photon Data

UPDATE HISTORY:
    Updated 12/2022: single implicit import of altimetry tools
    Updated 06/2022: update classify photons to match current GSFC version
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 10/2021: using python logging for handling verbose output
        do not use possible TEP photons in photon classification calculation
        added parsing for converting file lines to arguments
        use reference photon geolocation as default in case of no valid segments
    Updated 08/2021: update classify photons to match current GSFC version
    Updated 05/2021: add photon classifier based on GSFC YAPC algorithms
        move surface fit operations into separate module
    Updated 02/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: time utilities for converting times from JD and to decimal
    Updated 12/2020: H5py deprecation warning change to use make_scale
    Updated 10/2020: using argparse to set parameters
    Updated 09/2020: using reference photon delta time to interpolate ATL09
    Updated 08/2020: using convert delta time function to convert to Julian days
    Updated 07/2020: "re-tiding" is no longer unnecessary
    Updated 06/2020: verify that complementary beam pair is in list of beams
        set masks of output arrays after reading from HDF5
        add additional beam check within heights groups
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: adding segment quality summary variable
    Updated 04/2019: updated backup algorithm for when the surface fit fails
        estimate both mean and median first photon bias corrections
        estimate both mean and median transmit pulse shape corrections
    Updated 03/2019: extract a set of ATL09 parameters for each ATL03 segment_ID
    Updated 02/2019: procedures following ATBD for first ATL03 release
    Written 05/2017
"""
from __future__ import print_function, division

import sys
import os
import re
import h5py
import logging
import argparse
import datetime
import numpy as np
import scipy.signal
import scipy.interpolate
import sklearn.neighbors
import sklearn.cluster
from mpi4py import MPI
import icesat2_toolkit as is2tk
import yapc.classify_photons

# PURPOSE: keep track of MPI threads
def info(rank, size):
    logging.info(f'Rank {rank+1:d} of {size:d}')
    logging.info(f'module name: {__name__}')
    if hasattr(os, 'getppid'):
        logging.info(f'parent process: {os.getppid():d}')
    logging.info(f'process id: {os.getpid():d}')

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Read ICESat-2 ATL03 and ATL09 data files to calculate
            average segment surfaces
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = is2tk.utilities.convert_arg_line_to_args
    # command line parameters
    # first file listed contains the ATL03 file
    # second file listed is the associated ATL09 file
    parser.add_argument('ATL03',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='ICESat-2 ATL03 file to run')
    parser.add_argument('ATL09',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='?',
        help='ICESat-2 ATL09 file to run')
    # use default output file name
    parser.add_argument('--output','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Name and path of output file')
    # verbosity settings
    # verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of output files')
    # return the parser
    return parser

# PURPOSE: read ICESat-2 geolocated photon height data (ATL03)
# and backscatter profiles/atmospheric layer characteristics (ATL09)
# Computes average heights over segments
def main():
    # start MPI communicator
    comm = MPI.COMM_WORLD

    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # create logger
    loglevel = logging.INFO if args.verbose else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # output module information for process
    info(comm.rank,comm.size)
    if (comm.rank == 0):
        logging.info(f'{args.ATL03} -->')
    # directory setup
    ATL03_dir = os.path.dirname(args.ATL03)

    # compile regular expression operator for extracting data from ATL03 files
    rx1 = re.compile(r'(processed_)?(ATL\d+)_(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})'
        r'(\d{2})_(\d{4})(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
    # universal variables
    # speed of light
    c = 299792458.0
    # associated beam pairs
    associated_beam_pair = dict(gt1l='gt1r',gt1r='gt1l',gt2l='gt2r',gt2r='gt2l',
        gt3l='gt3r',gt3r='gt3l')

    # read ICESat-2 ATL03 HDF5 files (extract base parameters)
    SUB,PRD,YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX=rx1.findall(args.ATL03).pop()

    # Open the HDF5 file for reading
    fileID = h5py.File(args.ATL03, 'r', driver='mpio', comm=comm)

    # read each input beam within the file
    IS2_atl03_beams = []
    for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
        # check if subsetted beam contains data
        # check in both the geolocation and heights groups
        try:
            fileID[gtx]['geolocation']['segment_id']
            fileID[gtx]['heights']['delta_time']
        except KeyError:
            pass
        else:
            IS2_atl03_beams.append(gtx)

    # number of GPS seconds between the GPS epoch
    # and ATLAS Standard Data Product (SDP) epoch
    atlas_sdp_gps_epoch = fileID['ancillary_data']['atlas_sdp_gps_epoch'][:]
    # which TEP to use for a given spot (convert to 0-based index)
    tep_valid_spot = fileID['ancillary_data']['tep']['tep_valid_spot'][:] - 1
    tep_pce = ['pce1_spot1','pce2_spot3']
    # valid range of times for each TEP histogram
    tep_range_prim = fileID['ancillary_data']['tep']['tep_range_prim'][:]
    # save tep parameters for a given beam
    tep = {}

    # variables of interest for generating corrected elevation estimates
    Segment_ID = {}
    Segment_Index_begin = {}
    Segment_PE_count = {}
    Segment_Distance = {}
    Segment_Length = {}
    Segment_Background = {}
    # fit parameters
    Segment_delta_time = {}
    Segment_Height = {}
    Segment_Land_Ice = {}
    Segment_dH_along = {}
    Segment_dH_across = {}
    Segment_Height_Error = {}
    Segment_Land_Ice_Error = {}
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
    Segment_Window = {}
    Segment_RDE = {}
    Segment_SNR = {}
    Segment_Photon_SNR = {}
    Segment_Summary = {}
    Segment_Iterations = {}
    Segment_Clusters = {}
    Segment_Source = {}
    Segment_Pulses = {}
    # correction parameters
    FPB_mean_corr = {}
    FPB_mean_sigma = {}
    FPB_median_corr = {}
    FPB_median_sigma = {}
    mean_dead_time = {}
    FPB_n_corr = {}
    FPB_cal_corr = {}
    TPS_mean_corr = {}
    TPS_median_corr = {}

    # for each input beam within the file
    for gtx in sorted(IS2_atl03_beams):
        logging.info(gtx) if (comm.rank == 0) else None
        # beam type (weak versus strong) for time
        atlas_beam_type = fileID[gtx].attrs['atlas_beam_type'].decode('utf-8')
        n_pixels = 16.0 if (atlas_beam_type == "strong") else 4.0
        # ATL03 Segment ID
        Segment_ID[gtx] = fileID[gtx]['geolocation']['segment_id'][:]
        # number of valid overlapping ATL03 segments
        n_seg = len(Segment_ID[gtx]) - 1
        # number of photon events
        n_pe, = fileID[gtx]['heights']['delta_time'].shape

        # first photon in the segment (convert to 0-based indexing)
        Segment_Index_begin[gtx] = fileID[gtx]['geolocation']['ph_index_beg'][:] - 1
        # number of photon events in the segment
        Segment_PE_count[gtx] = fileID[gtx]['geolocation']['segment_ph_cnt'][:]
        # along-track distance for each ATL03 segment
        Segment_Distance[gtx] = fileID[gtx]['geolocation']['segment_dist_x'][:]
        # along-track length for each ATL03 segment
        Segment_Length[gtx] = fileID[gtx]['geolocation']['segment_length'][:]
        # ocean tide
        fv = fileID[gtx]['geophys_corr']['tide_ocean'].attrs['_FillValue']
        tide_ocean = np.ma.array(fileID[gtx]['geophys_corr']['tide_ocean'][:],
            fill_value=fv)
        tide_ocean.mask = tide_ocean.data == tide_ocean.fill_value
        # interpolate background photon rate based on 50-shot summation
        background_delta_time = fileID[gtx]['bckgrd_atlas']['delta_time'][:]
        SPL = scipy.interpolate.UnivariateSpline(background_delta_time,
            fileID[gtx]['bckgrd_atlas']['bckgrd_rate'][:],k=3,s=0)
        Segment_Background[gtx] = SPL(fileID[gtx]['geolocation']['delta_time'][:])

        # ATLAS spot number for beam in current orientation
        spot = int(fileID[gtx].attrs['atlas_spot_number'])
        # get ATLAS impulse response variables for the transmitter echo path (TEP)
        tep1,tep2 = ('atlas_impulse_response','tep_histogram')
        # get appropriate transmitter-echo-path histogram for spot
        associated_pce = tep_valid_spot[spot-1]
        pce = tep_pce[associated_pce]
        # delta time of TEP histogram
        tep_tod, = fileID[tep1][pce][tep2]['tep_tod'][:]
        # truncate tep to primary histogram (reflection 43-50 ns)
        # and extract signal tep from noise tep.  calculate width of tep
        # ATL03 recommends subsetting between 15-30 ns to avoid secondary
        tep_hist_time = np.copy(fileID[tep1][pce][tep2]['tep_hist_time'][:])
        tep_hist = np.copy(fileID[tep1][pce][tep2]['tep_hist'][:])
        t_TX,p_TX,W_TX,FWHM,TXs,TXe = is2tk.fit.extract_tep_histogram(
            tep_hist_time, tep_hist, tep_range_prim)
        # save tep information and statistics
        tep[gtx] = {}
        tep[gtx]['pce'] = pce
        tep[gtx]['tep_tod'] = tep_tod
        tep[gtx]['tx_start'] = TXs
        tep[gtx]['tx_end'] = TXe
        tep[gtx]['tx_robust_sprd'] = W_TX
        tep[gtx]['sigma_tx'] = FWHM

        # channel dead time and first photon bias table for beam
        cal1,cal2 = ('ancillary_data','calibrations')
        channel_dead_time = fileID[cal1][cal2]['dead_time'][gtx]['dead_time'][:]
        mean_dead_time[gtx] = np.mean(channel_dead_time)
        fpb_dead_time = fileID[cal1][cal2]['first_photon_bias'][gtx]['dead_time'][:]
        fpb_strength = fileID[cal1][cal2]['first_photon_bias'][gtx]['strength'][:]
        fpb_width = fileID[cal1][cal2]['first_photon_bias'][gtx]['width'][:]
        fpb_corr = fileID[cal1][cal2]['first_photon_bias'][gtx]['ffb_corr'][:]
        # calculate first photon bias as a function of strength and width
        # for the calculated mean dead time of the beam
        ndt,ns,nw = np.shape(fpb_corr)
        fpb_corr_dead_time = np.zeros((ns,nw))
        for s in range(ns):
            for w in range(nw):
                SPL = scipy.interpolate.UnivariateSpline(fpb_dead_time/1e9,
                    fpb_corr[:,s,w],k=3,s=0)
                fpb_corr_dead_time[s,w] = SPL(mean_dead_time[gtx])
        # bivariate spline for estimating first-photon bias using CAL-19
        CAL19 = scipy.interpolate.RectBivariateSpline(fpb_strength[0,:],
            fpb_width[0,:]/1e9, fpb_corr_dead_time/1e12, kx=1, ky=1)

        # allocate for output segment fit data
        fill_value = fileID[gtx]['geolocation']['sigma_h'].attrs['_FillValue']
        # delta time of fit photons
        # use reference photon delta time as initial guess
        # will fill in with segment center delta times
        Distributed_delta_time = is2tk.fit.segment_mean(
            fileID[gtx]['geolocation']['delta_time'][:])
        # segment fit heights
        Distributed_Height = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Height.mask = np.ones((n_seg),dtype=bool)
        # land ice height corrected for first photon bias and transmit-pulse shape
        Distributed_Land_Ice = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Land_Ice.mask = np.ones((n_seg),dtype=bool)
        # segment fit along-track slopes
        Distributed_dH_along = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_along.mask = np.ones((n_seg),dtype=bool)
        # segment fit height errors
        Distributed_Height_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Height_Error.mask = np.ones((n_seg),dtype=bool)
        # land ice height errors (max of fit or first photon bias uncertainties)
        Distributed_Land_Ice_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Land_Ice_Error.mask = np.ones((n_seg),dtype=bool)
        # segment fit along-track slope errors
        Distributed_dH_along_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_along_Error.mask = np.ones((n_seg),dtype=bool)
        # difference between the mean and median of the residuals from fit height
        Distributed_Mean_Median = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Mean_Median.mask = np.ones((n_seg),dtype=bool)
        # along-track X coordinates of segment fit
        Distributed_X_atc = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_X_atc.mask = np.ones((n_seg),dtype=bool)
        # along-track X coordinate spread of points used in segment fit
        Distributed_X_spread = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_X_spread.mask = np.ones((n_seg),dtype=bool)
        # along-track Y coordinates of segment fit
        Distributed_Y_atc = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Y_atc.mask = np.ones((n_seg),dtype=bool)
        # longitude of fit photons
        Distributed_Longitude = is2tk.fit.segment_mean(
            fileID[gtx]['geolocation']['reference_photon_lon'][:])
        # latitude of fit photons
        Distributed_Latitude = is2tk.fit.segment_mean(
            fileID[gtx]['geolocation']['reference_photon_lat'][:])
        # number of photons in fit
        Distributed_N_Fit = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Distributed_N_Fit.mask = np.ones((n_seg),dtype=bool)
        # size of the window used in the fit
        Distributed_Window = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_Window.mask = np.ones((n_seg),dtype=bool)
        # robust dispersion estimator
        Distributed_RDE = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_RDE.mask = np.ones((n_seg),dtype=bool)
        # signal-to-noise ratio
        Distributed_SNR = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_SNR.mask = np.ones((n_seg),dtype=bool)
        # maximum signal-to-noise ratio from photon classifier
        Distributed_Photon_SNR = np.ma.zeros((n_seg),fill_value=0,dtype=int)
        Distributed_Photon_SNR.mask = np.ones((n_seg),dtype=bool)
        # segment quality summary
        Distributed_Summary = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Distributed_Summary.mask = np.ones((n_seg),dtype=bool)
        # number of iterations for fit
        Distributed_Iterations = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Distributed_Iterations.mask = np.ones((n_seg),dtype=bool)
        # number of estimated clusters of data
        Distributed_Clusters = np.ma.zeros((n_seg),fill_value=0,dtype=int)
        Distributed_Clusters.mask = np.ones((n_seg),dtype=bool)
        # signal source selection
        Distributed_Source = np.ma.zeros((n_seg),fill_value=4,dtype=int)
        Distributed_Source.mask = np.ones((n_seg),dtype=bool)
        # number of pulses in segment
        Distributed_Pulses = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Distributed_Pulses.mask = np.ones((n_seg),dtype=bool)
        # first photon bias estimates
        Distributed_FPB_mean_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_mean_corr.mask = np.ones((n_seg),dtype=bool)
        Distributed_FPB_mean_sigma = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_mean_sigma.mask = np.ones((n_seg),dtype=bool)
        Distributed_FPB_median_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_median_corr.mask = np.ones((n_seg),dtype=bool)
        Distributed_FPB_median_sigma = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_median_sigma.mask = np.ones((n_seg),dtype=bool)
        Distributed_FPB_n_corr = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Distributed_FPB_n_corr.mask = np.ones((n_seg),dtype=bool)
        Distributed_FPB_cal_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_FPB_cal_corr.mask = np.ones((n_seg),dtype=bool)
        # transmit pulse shape bias estimates
        Distributed_TPS_mean_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_TPS_mean_corr.mask = np.ones((n_seg),dtype=bool)
        Distributed_TPS_median_corr = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_TPS_median_corr.mask = np.ones((n_seg),dtype=bool)

        # along-track and across-track distance for photon events
        x_atc = fileID[gtx]['heights']['dist_ph_along'][:].copy()
        y_atc = fileID[gtx]['heights']['dist_ph_across'][:].copy()
        # photon event heights
        h_ph = fileID[gtx]['heights']['h_ph'][:].copy()
        # digital elevation model interpolated to photon events
        dem_h = np.zeros((n_pe))
        # for each 20m segment
        for j,_ in enumerate(Segment_ID[gtx]):
            # index for 20m segment j
            idx = Segment_Index_begin[gtx][j]
            # skip segments with no photon events
            if (idx < 0):
                continue
            # number of photons in 20m segment
            cnt = Segment_PE_count[gtx][j]
            # add segment distance to along-track coordinates
            x_atc[idx:idx+cnt] += Segment_Distance[gtx][j]
            # interpolate digital elevation model to photon events
            dem_h[idx:idx+cnt] = fileID[gtx]['geophys_corr']['dem_h'][j]

        # iterate over ATLAS major frames
        photon_mframes = fileID[gtx]['heights']['pce_mframe_cnt'][:].copy()
        # background ATLAS group variables are based upon 50-shot summations
        # PCE Major Frames are based upon 200-shot summations
        pce_mframe_cnt = fileID[gtx]['bckgrd_atlas']['pce_mframe_cnt'][:].copy()
        # find unique major frames and their indices within background ATLAS group
        # (there will 4 background ATLAS time steps for nearly every major frame)
        unique_major_frames,unique_index = np.unique(pce_mframe_cnt,return_index=True)
        # number of unique major frames in granule for beam
        major_frame_count = len(unique_major_frames)
        # height of each telemetry band for a major frame
        tlm_height = {}
        tlm_height['band1'] = fileID[gtx]['bckgrd_atlas']['tlm_height_band1'][:].copy()
        tlm_height['band2'] = fileID[gtx]['bckgrd_atlas']['tlm_height_band2'][:].copy()
        # elevation above ellipsoid of each telemetry band for a major frame
        tlm_top = {}
        tlm_top['band1'] = fileID[gtx]['bckgrd_atlas']['tlm_top_band1'][:].copy()
        tlm_top['band2'] = fileID[gtx]['bckgrd_atlas']['tlm_top_band2'][:].copy()
        # buffer to telemetry band to set as valid
        tlm_buffer = 100.0
        # flag denoting photon events as possible TEP
        if (int(RL) < 4):
            isTEP = np.any((fileID[gtx]['heights']['signal_conf_ph'][:]==-2),axis=1)
        else:
            isTEP = (fileID[gtx]['heights']['quality_ph'][:] == 3)
        # photon event weights
        Distributed_Weights = np.zeros((n_pe),dtype=np.float64)
        # run for each major frame (distributed over comm.size # of processes)
        for iteration in range(comm.rank, major_frame_count, comm.size):
            # background atlas index for iteration
            idx = unique_index[iteration]
            # photon indices for major frame (buffered by 1 frame on each side)
            # do not use possible TEP photons in photon classification
            i1, = np.nonzero((photon_mframes >= unique_major_frames[iteration]-1) &
                (photon_mframes <= unique_major_frames[iteration]+1) &
                np.logical_not(isTEP))
            # indices for the major frame within the buffered window
            i2, = np.nonzero(photon_mframes[i1] == unique_major_frames[iteration])
            # sum of telemetry band widths for major frame
            h_win_width = 0.0
            # check that each telemetry band is close to DEM
            for b in ['band1','band2']:
                # bottom of the telemetry band for major frame
                tlm_bot_band = tlm_top[b][idx] - tlm_height[b][idx]
                if np.any((dem_h[i1[i2]] >= (tlm_bot_band-tlm_buffer)) &
                    (dem_h[i1[i2]] <= (tlm_top[b][idx]+tlm_buffer))):
                    # add telemetry height to window width
                    h_win_width += tlm_height[b][idx]
            # calculate photon event weights
            Distributed_Weights[i1[i2]] = yapc.classify_photons(x_atc[i1],
                h_ph[i1], h_win_width, i2, K=0, min_knn=5, min_ph=3,
                min_xspread=1.0, min_hspread=0.01, win_x=15.0, win_h=6.0,
                method='linear')
        # photon event weights
        pe_weights = np.zeros((n_pe),dtype=np.float64)
        comm.Allreduce(sendbuf=[Distributed_Weights, MPI.DOUBLE], \
            recvbuf=[pe_weights, MPI.DOUBLE], op=MPI.SUM)
        Distributed_Weights = None
        # wait for all distributed processes to finish for beam
        comm.Barrier()

        # photon event weights scaled to a single byte
        weight_ph = np.array(255*pe_weights,dtype=np.uint8)
        # verify photon event weights
        np.clip(weight_ph, 0, 255, out=weight_ph)

        # iterate over valid ATL03 segments
        # in ATL03 1-based indexing: invalid == 0
        # here in 0-based indexing: invalid == -1
        segment_indices, = np.nonzero((Segment_Index_begin[gtx][:-1] >= 0) &
            (Segment_Index_begin[gtx][1:] >= 0))
        iteration_count = len(segment_indices)
        # run for each geoseg (distributed over comm.size # of processes)
        for iteration in range(comm.rank, iteration_count, comm.size):
            # indice for iteration (can run through a subset of segments)
            j = segment_indices[iteration]

            # iterate over valid ATL03 segments
            # in ATL03 1-based indexing: invalid == 0
            # here in 0-based indexing: invalid == -1
            if (Segment_Index_begin[gtx][j] >= 0):
                # index for segment j
                idx = Segment_Index_begin[gtx][j]
                # number of photons in segment (use 2 ATL03 segments)
                c1 = int(Segment_PE_count[gtx][j])
                c2 = int(Segment_PE_count[gtx][j+1])
                cnt = c1 + c2
                # time of each Photon event (PE)
                segment_times = np.copy(fileID[gtx]['heights']['delta_time'][idx:idx+cnt])
                # Photon event lat/lon and elevation (re-tided WGS84)
                segment_heights = np.copy(h_ph[idx:idx+cnt])
                # ATL03 pe heights no longer apply the ocean tide
                # and so "re-tiding" is no longer unnecessary
                # segment_heights[:c1] += tide_ocean[j]
                # segment_heights[c1:] += tide_ocean[j+1]
                segment_lats = np.copy(fileID[gtx]['heights']['lat_ph'][idx:idx+cnt])
                segment_lons = np.copy(fileID[gtx]['heights']['lon_ph'][idx:idx+cnt])
                # Photon event channel and identification
                ID_channel = np.copy(fileID[gtx]['heights']['ph_id_channel'][idx:idx+cnt])
                ID_pulse = np.copy(fileID[gtx]['heights']['ph_id_pulse'][idx:idx+cnt])
                n_pulses = np.unique(ID_pulse).__len__()
                frame_number = np.copy(fileID[gtx]['heights']['pce_mframe_cnt'][idx:idx+cnt])
                # vertical noise-photon density
                background_density = 2.0*n_pulses*Segment_Background[gtx][j]/c
                # along-track X and Y coordinates
                distance_along_X = np.copy(x_atc[idx:idx+cnt])
                distance_along_Y = np.copy(y_atc[idx:idx+cnt])
                # check the spread of photons along-track (must be > 20m)
                along_X_spread = distance_along_X.max() - distance_along_X.min()
                # check confidence level associated with each photon event
                # -2: TEP
                # -1: Events not associated with a specific surface type
                #  0: noise
                #  1: buffer but algorithm classifies as background
                #  2: low
                #  3: medium
                #  4: high
                # Surface types for signal classification confidence
                # 0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater
                ice_sig_conf = np.copy(fileID[gtx]['heights']['signal_conf_ph'][idx:idx+cnt,3])
                ice_sig_low_count = np.count_nonzero(ice_sig_conf > 1)
                # indices of TEP classified photons
                ice_sig_tep_pe, = np.nonzero(ice_sig_conf == -2)
                # photon event weights from photon classifier
                segment_weights = weight_ph[idx:idx+cnt]
                snr_norm = np.max(segment_weights)
                # photon event signal-to-noise ratio from photon classifier
                photon_snr = np.array(100.0*segment_weights/snr_norm,dtype=int)
                Distributed_Photon_SNR.data[j] = np.copy(snr_norm)
                Distributed_Photon_SNR.mask[j] = (snr_norm > 0)
                # photon confidence levels from classifier
                pe_sig_conf = np.zeros((cnt),dtype=int)
                # calculate confidence levels from photon classifier
                pe_sig_conf[photon_snr >= 25] = 2
                pe_sig_conf[photon_snr >= 60] = 3
                pe_sig_conf[photon_snr >= 80] = 4
                # copy classification for TEP photons
                pe_sig_conf[ice_sig_tep_pe] = -2
                pe_sig_low_count = np.count_nonzero(pe_sig_conf > 1)
                # check if segment has photon events classified for land ice
                # that are at or above low-confidence threshold
                # and that the spread of photons is greater than 20m
                if (pe_sig_low_count > 10) & (along_X_spread > 20):
                    # use density-based spatial clustering in segment
                    db = sklearn.cluster.DBSCAN(eps=0.5).fit(
                        np.c_[distance_along_X, segment_heights],
                        sample_weight=photon_snr)
                    labels = db.labels_
                    # number of noise photons
                    noise_photons = list(labels).count(-1)
                    noise_cluster = 1 if noise_photons else 0
                    # number of photon event clusters in segment
                    n_clusters = len(set(labels)) - noise_cluster
                    Distributed_Clusters.data[j] = n_clusters
                    Distributed_Clusters.mask[j] = (n_clusters > 0)
                    # perform a surface fit procedure
                    Segment_X = Segment_Distance[gtx][j] + Segment_Length[gtx][j]
                    valid,fit,centroid = is2tk.fit.try_surface_fit(
                        distance_along_X, distance_along_Y, segment_heights,
                        pe_sig_conf, Segment_X, SURF_TYPE='linear', ITERATE=20,
                        CONFIDENCE=[1,0])
                    # indices of points used in final iterated fit
                    ifit = fit['indices'] if valid else None
                    if bool(valid) & (np.abs(fit['error'][0]) < 20):
                        Distributed_Height.data[j] = fit['beta'][0]
                        Distributed_Height.mask[j] = False
                        Distributed_dH_along.data[j] = fit['beta'][1]
                        Distributed_dH_along.mask[j] = False
                        Distributed_Height_Error.data[j] = fit['error'][0]
                        Distributed_Height_Error.mask[j] = False
                        Distributed_dH_along_Error.data[j] = fit['error'][1]
                        Distributed_dH_along_Error.mask[j] = False
                        # along-track and cross-track coordinates
                        Distributed_X_atc.data[j] = np.copy(centroid['x'])
                        Distributed_X_atc.mask[j] = False
                        Distributed_X_spread.data[j] = np.copy(along_X_spread)
                        Distributed_X_spread.mask[j] = False
                        Distributed_Y_atc.data[j] = np.copy(centroid['y'])
                        Distributed_Y_atc.mask[j] = False
                        # fit geolocation to the along-track distance of segment
                        Distributed_delta_time[j] = \
                            is2tk.fit.fit_geolocation(segment_times[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Longitude[j] = \
                            is2tk.fit.fit_geolocation(segment_lons[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Latitude[j] = \
                            is2tk.fit.fit_geolocation(segment_lats[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        # number of photons used in fit
                        Distributed_N_Fit.data[j] = len(ifit)
                        Distributed_N_Fit.mask[j] = False
                        # size of the final window
                        Distributed_Window.data[j] = np.copy(fit['window'])
                        Distributed_Window.mask[j] = False
                        # robust dispersion estimator
                        Distributed_RDE.data[j] = np.copy(fit['RDE'])
                        Distributed_RDE.mask[j] = False
                        # signal to noise ratio
                        N_BG = background_density*Distributed_Window.data[j]
                        Distributed_SNR.data[j] = Distributed_N_Fit.data[j]/N_BG
                        Distributed_SNR.mask[j] = False
                        # number of iterations used in fit
                        Distributed_Iterations.data[j] = np.copy(fit['iterations'])
                        Distributed_Iterations.mask[j] = False
                        Distributed_Source.data[j] = np.copy(valid)
                        Distributed_Source.mask[j] = False
                        Distributed_Pulses.data[j] = np.copy(n_pulses)
                        Distributed_Pulses.mask[j] = False
                        # calculate residuals off of fit surface for all data
                        x_slope = Distributed_dH_along[j]*(distance_along_X-Distributed_X_atc[j])
                        height_residuals = segment_heights-Distributed_Height[j]-x_slope
                        temporal_residuals = -2.0*height_residuals/c
                        # calculate difference between the mean and the median from the fit
                        Distributed_Mean_Median.data[j] = np.mean(height_residuals[ifit]) - \
                            np.median(height_residuals[ifit])
                        Distributed_Mean_Median.mask[j] = False
                        # calculate flags for quality summary
                        VPD = Distributed_N_Fit.data[j]/Distributed_Window.data[j]
                        Distributed_Summary.data[j] = int(
                            (Distributed_RDE.data[j] >= 1) |
                            (Distributed_Height_Error.data[j] >= 1) |
                            (VPD <= (n_pixels/4.0)))
                        Distributed_Summary.mask[j] = False
                        # estimate first photon bias corrections
                        # step-size for histograms (50 ps ~ 7.5mm height)
                        ii, = np.nonzero((height_residuals >= -Distributed_Window.data[j]) &
                            (height_residuals <= Distributed_Window.data[j]))
                        try:
                            FPB = is2tk.fit.calc_first_photon_bias(
                                temporal_residuals[ii], n_pulses, n_pixels,
                                mean_dead_time[gtx], 5e-11, ITERATE=20)
                        except:
                            pass
                        else:
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
                            # first photon bias correction from CAL-19
                            FPB_calibrated = CAL19.ev(FPB['strength'],FPB['width'])
                            Distributed_FPB_cal_corr.data[j] = -0.5*FPB_calibrated*c
                            Distributed_FPB_cal_corr.mask[j] = False
                        # estimate transmit pulse shape correction
                        try:
                            W_RX = 2.0*Distributed_RDE.data[j]/c
                            dt_W = 2.0*Distributed_Window.data[j]/c
                            TPS = is2tk.fit.calc_transmit_pulse_shape(t_TX,
                                p_TX, W_TX, W_RX, dt_W, Distributed_SNR.data[j],
                                ITERATE=50)
                        except:
                            pass
                        else:
                            Distributed_TPS_mean_corr.data[j] = 0.5*TPS['mean']*c
                            Distributed_TPS_mean_corr.mask[j] = False
                            Distributed_TPS_median_corr.data[j] = 0.5*TPS['median']*c
                            Distributed_TPS_median_corr.mask[j] = False

            # some ATL03 segments will not result in a valid fit
            # backup algorithm uses 4 segments to find a valid surface
            if (j not in (0,n_seg-2,n_seg-1)) & Distributed_Height.mask[j] & \
                (Segment_Index_begin[gtx][j-1] > 0):
                # index for segment j
                idx = Segment_Index_begin[gtx][j-1]
                # number of photons in segment (use 4 ATL03 segments)
                c1 = Segment_PE_count[gtx][j-1].astype(int)
                c2 = Segment_PE_count[gtx][j].astype(int)
                c3 = Segment_PE_count[gtx][j+1].astype(int)
                c4 = Segment_PE_count[gtx][j+2].astype(int)
                cnt = c1 + c2 + c3 + c4
                # time of each Photon event (PE)
                segment_times = np.copy(fileID[gtx]['heights']['delta_time'][idx:idx+cnt])
                # Photon event lat/lon and elevation (re-tided WGS84)
                segment_heights = np.copy(h_ph[idx:idx+cnt])
                # ATL03 pe heights no longer apply the ocean tide
                # and so "re-tiding" is no longer unnecessary
                # segment_heights[:c1] += tide_ocean[j-1]
                # segment_heights[c1:c1+c2] += tide_ocean[j]
                # segment_heights[c1+c2:c1+c2+c3] += tide_ocean[j+1]
                # segment_heights[c1+c2+c3:] += tide_ocean[j+2]
                segment_lats = np.copy(fileID[gtx]['heights']['lat_ph'][idx:idx+cnt])
                segment_lons = np.copy(fileID[gtx]['heights']['lon_ph'][idx:idx+cnt])
                # Photon event channel and identification
                ID_channel = np.copy(fileID[gtx]['heights']['ph_id_channel'][idx:idx+cnt])
                ID_pulse = np.copy(fileID[gtx]['heights']['ph_id_pulse'][idx:idx+cnt])
                n_pulses = np.unique(ID_pulse).__len__()
                frame_number = np.copy(fileID[gtx]['heights']['pce_mframe_cnt'][idx:idx+cnt])
                # vertical noise-photon density
                background_density = 2.0*n_pulses*Segment_Background[gtx][j]/c
                # along-track X and Y coordinates
                distance_along_X = np.copy(x_atc[idx:idx+cnt])
                distance_along_Y = np.copy(y_atc[idx:idx+cnt])
                # check the spread of photons along-track (must be > 40m)
                along_X_spread = distance_along_X.max() - distance_along_X.min()
                # check confidence level associated with each photon event
                # -2: TEP
                # -1: Events not associated with a specific surface type
                #  0: noise
                #  1: buffer but algorithm classifies as background
                #  2: low
                #  3: medium
                #  4: high
                # Surface types for signal classification confidence
                # 0=Land; 1=Ocean; 2=SeaIce; 3=LandIce; 4=InlandWater
                ice_sig_conf = np.copy(fileID[gtx]['heights']['signal_conf_ph'][idx:idx+cnt,3])
                ice_sig_low_count = np.count_nonzero(ice_sig_conf > 1)
                # indices of TEP classified photons
                ice_sig_tep_pe, = np.nonzero(ice_sig_conf == -2)
                # photon event weights from photon classifier
                segment_weights = pe_weights[idx:idx+cnt]
                snr_norm = np.max(segment_weights)
                # photon event signal-to-noise ratio from photon classifier
                photon_snr = np.zeros((cnt),dtype=int)
                if (snr_norm > 0):
                    photon_snr[:] = 100.0*segment_weights/snr_norm
                # copy signal to noise ratio for segment
                Distributed_Photon_SNR.data[j] = np.copy(snr_norm)
                Distributed_Photon_SNR.mask[j] = (snr_norm > 0)
                # photon confidence levels from classifier
                pe_sig_conf = np.zeros((cnt),dtype=int)
                # calculate confidence levels from photon classifier
                pe_sig_conf[photon_snr >= 25] = 2
                pe_sig_conf[photon_snr >= 60] = 3
                pe_sig_conf[photon_snr >= 80] = 4
                # copy classification for TEP photons
                pe_sig_conf[ice_sig_tep_pe] = -2
                pe_sig_low_count = np.count_nonzero(pe_sig_conf > 1)
                # check if segment has photon events classified for land ice
                # that are at or above low-confidence threshold
                # and that the spread of photons is greater than 40m
                if (pe_sig_low_count > 10) & (along_X_spread > 40):
                    # use density-based spatial clustering in segment
                    db = sklearn.cluster.DBSCAN(eps=0.5).fit(
                        np.c_[distance_along_X, segment_heights],
                        sample_weight=photon_snr)
                    labels = db.labels_
                    # number of noise photons
                    noise_photons = list(labels).count(-1)
                    noise_cluster = 1 if noise_photons else 0
                    # number of photon event clusters in segment
                    n_clusters = len(set(labels)) - noise_cluster
                    Distributed_Clusters.data[j] = n_clusters
                    Distributed_Clusters.mask[j] = (n_clusters > 0)
                    # perform a surface fit procedure
                    Segment_X = Segment_Distance[gtx][j] + Segment_Length[gtx][j]
                    valid,fit,centroid = is2tk.fit.try_surface_fit(
                        distance_along_X, distance_along_Y, segment_heights,
                        pe_sig_conf, Segment_X, SURF_TYPE='quadratic',
                        ITERATE=20, CONFIDENCE=[0])
                    # indices of points used in final iterated fit
                    ifit = fit['indices'] if valid else None
                    if bool(valid) & (np.abs(fit['error'][0]) < 20):
                        Distributed_Height.data[j] = fit['beta'][0]
                        Distributed_Height.mask[j] = False
                        Distributed_dH_along.data[j] = fit['beta'][1]
                        Distributed_dH_along.mask[j] = False
                        Distributed_Height_Error.data[j] = fit['error'][0]
                        Distributed_Height_Error.mask[j] = False
                        Distributed_dH_along_Error.data[j] = fit['error'][1]
                        Distributed_dH_along_Error.mask[j] = False
                        # along-track and cross-track coordinates
                        Distributed_X_atc.data[j] = np.copy(centroid['x'])
                        Distributed_X_atc.mask[j] = False
                        Distributed_X_spread.data[j] = np.copy(along_X_spread)
                        Distributed_X_spread.mask[j] = False
                        Distributed_Y_atc.data[j] = np.copy(centroid['y'])
                        Distributed_Y_atc.mask[j] = False
                        # fit geolocation to the along-track distance of segment
                        Distributed_delta_time[j] = \
                            is2tk.fit.fit_geolocation(segment_times[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Longitude[j] = \
                            is2tk.fit.fit_geolocation(segment_lons[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        Distributed_Latitude[j] = \
                            is2tk.fit.fit_geolocation(segment_lats[ifit],
                            distance_along_X[ifit], Distributed_X_atc[j])
                        # number of photons used in fit
                        Distributed_N_Fit.data[j] = len(ifit)
                        Distributed_N_Fit.mask[j] = False
                        # size of the final window
                        Distributed_Window.data[j] = np.copy(fit['window'])
                        Distributed_Window.mask[j] = False
                        # robust dispersion estimator
                        Distributed_RDE.data[j] = np.copy(fit['RDE'])
                        Distributed_RDE.mask[j] = False
                        # signal to noise ratio
                        N_BG = background_density*Distributed_Window.data[j]
                        Distributed_SNR.data[j] = Distributed_N_Fit.data[j]/N_BG
                        Distributed_SNR.mask[j] = False
                        # number of iterations used in fit
                        Distributed_Iterations.data[j] = np.copy(fit['iterations'])
                        Distributed_Iterations.mask[j] = False
                        Distributed_Source.data[j] = 2 + np.copy(valid)
                        Distributed_Source.mask[j] = False
                        Distributed_Pulses.data[j] = np.copy(n_pulses)
                        Distributed_Pulses.mask[j] = False
                        # calculate residuals off of fit surface for all data
                        x_slope = Distributed_dH_along[j]*(distance_along_X-Distributed_X_atc[j])
                        height_residuals = segment_heights-Distributed_Height[j]-x_slope
                        temporal_residuals = -2.0*height_residuals/c
                        # calculate difference between the mean and the median from the fit
                        Distributed_Mean_Median.data[j] = np.mean(height_residuals[ifit]) - \
                            np.median(height_residuals[ifit])
                        Distributed_Mean_Median.mask[j] = False
                        # calculate flags for quality summary
                        VPD = Distributed_N_Fit.data[j]/Distributed_Window.data[j]
                        Distributed_Summary.data[j] = int(
                            (Distributed_RDE.data[j] >= 1) |
                            (Distributed_Height_Error.data[j] >= 1) |
                            (VPD <= (n_pixels/4.0)))
                        Distributed_Summary.mask[j] = False
                        # estimate first photon bias corrections
                        # step-size for histograms (50 ps ~ 7.5mm height)
                        try:
                            ii, = np.nonzero((height_residuals >= -Distributed_Window.data[j]) &
                                (height_residuals <= Distributed_Window.data[j]))
                            FPB = is2tk.fit.calc_first_photon_bias(
                                temporal_residuals[ii], n_pulses, n_pixels,
                                mean_dead_time[gtx], 5e-11, ITERATE=20)
                        except:
                            pass
                        else:
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
                            # first photon bias correction from CAL-19
                            FPB_calibrated = CAL19.ev(FPB['strength'],FPB['width'])
                            Distributed_FPB_cal_corr.data[j] = -0.5*FPB_calibrated*c
                            Distributed_FPB_cal_corr.mask[j] = False
                        # estimate transmit pulse shape correction
                        try:
                            W_RX = 2.0*Distributed_RDE.data[j]/c
                            dt_W = 2.0*Distributed_Window.data[j]/c
                            TPS = is2tk.fit.calc_transmit_pulse_shape(t_TX,
                                p_TX, W_TX, W_RX, dt_W, Distributed_SNR.data[j],
                                ITERATE=50)
                        except:
                            pass
                        else:
                            Distributed_TPS_mean_corr.data[j] = 0.5*TPS['mean']*c
                            Distributed_TPS_mean_corr.mask[j] = False
                            Distributed_TPS_median_corr.data[j] = 0.5*TPS['median']*c
                            Distributed_TPS_median_corr.mask[j] = False

            # if there is a valid land ice height
            if (~Distributed_Height.mask[j]):
                # land ice height corrected for first photon bias and transmit-pulse shape
                # segment heights have already been "re-tided"
                Distributed_Land_Ice.data[j] = Distributed_Height.data[j] + \
                    Distributed_FPB_median_corr.data[j] + Distributed_TPS_median_corr.data[j]
                Distributed_Land_Ice.mask[j] = False
                # land ice height errors (max of fit or first photon bias uncertainties)
                Distributed_Land_Ice_Error.data[j] = np.sqrt(np.max([
                    Distributed_Height_Error.data[j]**2,
                    Distributed_FPB_median_sigma.data[j]**2]))
                Distributed_Land_Ice_Error.mask[j] = False

        # communicate output MPI matrices between ranks
        # operations are element summations and logical "and" across elements

        # delta time of fit photons
        Segment_delta_time[gtx] = np.zeros((n_seg),fill_value=fill_value)
        comm.Allreduce(sendbuf=[Distributed_delta_time, MPI.DOUBLE], \
            recvbuf=[Segment_delta_time[gtx], MPI.DOUBLE], op=MPI.SUM)
        Distributed_delta_time = None
        # segment fit heights
        Segment_Height[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Height[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Height.data, MPI.DOUBLE], \
            recvbuf=[Segment_Height[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Height.mask, MPI.BOOL], \
            recvbuf=[Segment_Height[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Height = None
        # land ice height corrected for first photon bias and transmit-pulse shape
        Segment_Land_Ice[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Land_Ice[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice.data, MPI.DOUBLE], \
            recvbuf=[Segment_Land_Ice[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice.mask, MPI.BOOL], \
            recvbuf=[Segment_Land_Ice[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Land_Ice = None
        # segment fit along-track slopes
        Segment_dH_along[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_along[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_dH_along.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_along[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_along.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_along[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_along = None
        # segment fit height errors
        Segment_Height_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Height_Error[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Height_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Height_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Height_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Height_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Height_Error = None
        # land ice height errors (max of fit or first photon bias uncertainties)
        Segment_Land_Ice_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Land_Ice_Error[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_Land_Ice_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Land_Ice_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_Land_Ice_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Land_Ice_Error = None
        # segment fit along-track slope errors
        Segment_dH_along_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_along_Error[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_dH_along_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_along_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_along_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_along_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_along_Error = None
        # difference between the mean and median of the residuals from fit height
        Segment_Mean_Median[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Mean_Median[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Mean_Median.data, MPI.DOUBLE], \
            recvbuf=[Segment_Mean_Median[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Mean_Median.mask, MPI.BOOL], \
            recvbuf=[Segment_Mean_Median[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Mean_Median = None
        # along-track X coordinates of segment fit
        Segment_X_atc[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_X_atc[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_X_atc.data, MPI.DOUBLE], \
            recvbuf=[Segment_X_atc[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_X_atc.mask, MPI.BOOL], \
            recvbuf=[Segment_X_atc[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_X_atc = None
        # along-track X coordinate spread of points used in segment fit
        Segment_X_spread[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_X_spread[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_X_spread.data, MPI.DOUBLE], \
            recvbuf=[Segment_X_spread[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_X_spread.mask, MPI.BOOL], \
            recvbuf=[Segment_X_spread[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_X_spread = None
        # along-track Y coordinates of segment fit
        Segment_Y_atc[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Y_atc[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Y_atc.data, MPI.DOUBLE], \
            recvbuf=[Segment_Y_atc[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Y_atc.mask, MPI.BOOL], \
            recvbuf=[Segment_Y_atc[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Y_atc = None
        # longitude of fit photons
        Segment_Longitude[gtx] = np.zeros((n_seg),fill_value=fill_value)
        comm.Allreduce(sendbuf=[Distributed_Longitude, MPI.DOUBLE], \
            recvbuf=[Segment_Longitude[gtx], MPI.DOUBLE], op=MPI.SUM)
        # latitude of fit photons
        Segment_Latitude[gtx] = np.zeros((n_seg),fill_value=fill_value)
        comm.Allreduce(sendbuf=[Distributed_Latitude, MPI.DOUBLE], \
            recvbuf=[Segment_Latitude[gtx], MPI.DOUBLE], op=MPI.SUM)
        Distributed_Latitude = None
        # number of photons in fit
        Segment_N_Fit[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Segment_N_Fit[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_N_Fit.data, MPI.INT], \
            recvbuf=[Segment_N_Fit[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_N_Fit.mask, MPI.BOOL], \
            recvbuf=[Segment_N_Fit[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_N_Fit = None
        # size of the window used in the fit
        Segment_Window[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_Window[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Window.data, MPI.DOUBLE], \
            recvbuf=[Segment_Window[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Window.mask, MPI.BOOL], \
            recvbuf=[Segment_Window[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Window = None
        # robust dispersion estimator
        Segment_RDE[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_RDE[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_RDE.data, MPI.DOUBLE], \
            recvbuf=[Segment_RDE[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_RDE.mask, MPI.BOOL], \
            recvbuf=[Segment_RDE[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_RDE = None
        # signal-to-noise ratio
        Segment_SNR[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_SNR[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_SNR.data, MPI.DOUBLE], \
            recvbuf=[Segment_SNR[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_SNR.mask, MPI.BOOL], \
            recvbuf=[Segment_SNR[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_SNR = None
        # photon event signal-to-noise ratio from photon classifier
        Segment_Photon_SNR[gtx] = np.ma.zeros((n_seg),fill_value=0,dtype=int)
        Segment_Photon_SNR[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Photon_SNR.data, MPI.INT], \
            recvbuf=[Segment_Photon_SNR[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Photon_SNR.mask, MPI.BOOL], \
            recvbuf=[Segment_Photon_SNR[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Photon_SNR = None
        # segment quality summary
        Segment_Summary[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Segment_Summary[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Summary.data, MPI.INT], \
            recvbuf=[Segment_Summary[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Summary.mask, MPI.BOOL], \
            recvbuf=[Segment_Summary[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Summary = None
        # number of iterations for fit
        Segment_Iterations[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Segment_Iterations[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Iterations.data, MPI.INT], \
            recvbuf=[Segment_Iterations[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Iterations.mask, MPI.BOOL], \
            recvbuf=[Segment_Iterations[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Iterations = None
        # number of photon event clusters
        Segment_Clusters[gtx] = np.ma.zeros((n_seg),fill_value=0,dtype=int)
        Segment_Clusters[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Clusters.data, MPI.INT], \
            recvbuf=[Segment_Clusters[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Clusters.mask, MPI.BOOL], \
            recvbuf=[Segment_Clusters[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Clusters = None
        # signal source selection
        Segment_Source[gtx] = np.ma.zeros((n_seg),fill_value=4,dtype=int)
        Segment_Source[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Source.data, MPI.INT], \
            recvbuf=[Segment_Source[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Source.mask, MPI.BOOL], \
            recvbuf=[Segment_Source[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Source = None
        # number of pulses in segment
        Segment_Pulses[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        Segment_Pulses[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_Pulses.data, MPI.INT], \
            recvbuf=[Segment_Pulses[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_Pulses.mask, MPI.BOOL], \
            recvbuf=[Segment_Pulses[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_Pulses = None
        # first photon bias estimates
        FPB_mean_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_mean_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_mean_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_mean_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_mean_corr = None
        FPB_mean_sigma[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_mean_sigma[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_sigma.data, MPI.DOUBLE], \
            recvbuf=[FPB_mean_sigma[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_mean_sigma.mask, MPI.BOOL], \
            recvbuf=[FPB_mean_sigma[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_mean_sigma = None
        FPB_median_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_median_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_median_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_median_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_median_corr = None
        FPB_median_sigma[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_median_sigma[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_sigma.data, MPI.DOUBLE], \
            recvbuf=[FPB_median_sigma[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_median_sigma.mask, MPI.BOOL], \
            recvbuf=[FPB_median_sigma[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_median_sigma = None
        FPB_n_corr[gtx] = np.ma.zeros((n_seg),fill_value=-1,dtype=int)
        FPB_n_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_n_corr.data, MPI.INT], \
            recvbuf=[FPB_n_corr[gtx].data, MPI.INT], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_n_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_n_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_n_corr = None
        FPB_cal_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        FPB_cal_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_FPB_cal_corr.data, MPI.DOUBLE], \
            recvbuf=[FPB_cal_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_FPB_cal_corr.mask, MPI.BOOL], \
            recvbuf=[FPB_cal_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_FPB_cal_corr = None
        # transmit pulse shape bias estimates
        TPS_mean_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        TPS_mean_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_TPS_mean_corr.data, MPI.DOUBLE], \
            recvbuf=[TPS_mean_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_TPS_mean_corr.mask, MPI.BOOL], \
            recvbuf=[TPS_mean_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_TPS_mean_corr = None
        TPS_median_corr[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        TPS_median_corr[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_TPS_median_corr.data, MPI.DOUBLE], \
            recvbuf=[TPS_median_corr[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_TPS_median_corr.mask, MPI.BOOL], \
            recvbuf=[TPS_median_corr[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_TPS_median_corr = None
        # wait for all distributed processes to finish for beam
        comm.Barrier()

    # copy variables for outputting to HDF5 file
    IS2_atl03_fit = {}
    IS2_atl03_fill = {}
    IS2_atl03_attrs = {}

    # ICESat-2 spacecraft orientation at time
    IS2_atl03_fit['orbit_info'] = {}
    IS2_atl03_attrs['orbit_info'] = {}
    for key,val in fileID['orbit_info'].items():
        IS2_atl03_fit['orbit_info'][key] = val[:]
        # Getting attributes of group and included variables
        # Global Group Attributes
        for att_name,att_val in fileID['orbit_info'].attrs.items():
            IS2_atl03_attrs['orbit_info'][att_name] = att_val
        # Variable Attributes
        IS2_atl03_attrs['orbit_info'][key] = {}
        for att_name,att_val in val.attrs.items():
            IS2_atl03_attrs['orbit_info'][key][att_name] = att_val

    # information ancillary to the data product
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    # Add this value to delta time parameters to compute full gps_seconds
    # could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    # and add leap seconds since 2018-01-01T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl03_fit['ancillary_data'] = {}
    IS2_atl03_attrs['ancillary_data'] = {}
    for key in ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc','end_cycle',
        'end_geoseg','end_gpssow','end_gpsweek','end_orbit','end_region',
        'end_rgt','granule_end_utc','granule_start_utc','release','start_cycle',
        'start_geoseg','start_gpssow','start_gpsweek','start_orbit','start_region',
        'start_rgt','version']:
        # get each HDF5 variable
        IS2_atl03_fit['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        # Getting attributes of group and included variables
        IS2_atl03_attrs['ancillary_data'][key] = {}
        for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
            IS2_atl03_attrs['ancillary_data'][key][att_name] = att_val

    # for each output beam
    for gtx in sorted(IS2_atl03_beams):
        # atmospheric profile for beam gtx from ATL09 dataset
        pfl = fileID[gtx].attrs['atmosphere_profile']
        # complementary beam in pair
        cmp = associated_beam_pair[gtx]
        # extract and interpolate atmospheric parameters from ATL09
        dtime = fileID[gtx]['geolocation']['delta_time'][:]
        IS2_atl09_mds,IS2_atl09_attrs = read_HDF5_ATL09(args.ATL09, pfl,
            dtime, ATTRIBUTES=True, COMM=comm)

        # segment fit across-track slopes
        Distributed_dH_across = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_across.mask = np.ones((n_seg),dtype=bool)
        # segment fit across-track slope errors
        Distributed_dH_across_Error = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_dH_across_Error.mask = np.ones((n_seg),dtype=bool)
        # contribution of geolocation uncertainty to height error
        Distributed_sigma_geo = np.ma.zeros((n_seg),fill_value=fill_value)
        Distributed_sigma_geo.mask = np.ones((n_seg),dtype=bool)

        # iterate over valid ATL03 segments
        # in ATL03 1-based indexing: invalid == 0
        # here in 0-based indexing: invalid == -1
        segment_indices, = np.nonzero((Segment_Index_begin[gtx][:-1] >= 0) &
            (Segment_Index_begin[gtx][1:] >= 0))
        # verify that complementary beam pair is in list of beams
        iteration_count = len(segment_indices) if (cmp in IS2_atl03_beams) else 0
        # run for each geoseg (distributed over comm.size # of processes)
        for iteration in range(comm.rank, iteration_count, comm.size):
            # indice for iteration (can run through a subset of segments)
            j = segment_indices[iteration]
            # across track slopes for beam
            if ((~Segment_Height[gtx].mask[j]) & (~Segment_Height[cmp].mask[j])):
                # segment fit across-track slopes
                dY = (Segment_Y_atc[gtx].data[j] - Segment_Y_atc[cmp].data[j])
                Distributed_dH_across.data[j] = (Segment_Land_Ice[gtx].data[j] -
                    Segment_Land_Ice[cmp].data[j])/dY
                Distributed_dH_across.mask[j] = False
                # segment fit across-track slope errors
                Distributed_dH_across_Error.data[j] = np.sqrt(
                    Segment_Land_Ice_Error[gtx].data[j]**2 +
                    Segment_Land_Ice_Error[cmp].data[j]**2)/np.abs(dY)
                Distributed_dH_across_Error.mask[j] = False
                # geolocation uncertainty
                sigma_geo_across = fileID[gtx]['geolocation']['sigma_across'][j]
                sigma_geo_along = fileID[gtx]['geolocation']['sigma_along'][j]
                sigma_geo_h = fileID[gtx]['geolocation']['sigma_h'][j]
                # contribution of geolocation uncertainty to height errors
                Distributed_sigma_geo.data[j] = np.sqrt(sigma_geo_h**2 +
                    (sigma_geo_along*Segment_dH_along[gtx].data[j])**2 +
                    (sigma_geo_across*Distributed_dH_across.data[j])**2)
                Distributed_sigma_geo.mask[j] = False

        # segment fit across-track slopes
        Segment_dH_across[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_across[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_dH_across.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_across[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_across.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_across[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_across = None
        # segment fit across-track slope errors
        Segment_dH_across_Error[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_dH_across_Error[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_dH_across_Error.data, MPI.DOUBLE], \
            recvbuf=[Segment_dH_across_Error[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_dH_across_Error.mask, MPI.BOOL], \
            recvbuf=[Segment_dH_across_Error[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_dH_across_Error = None
        # contribution of geolocation uncertainty to height errors
        Segment_sigma_geo[gtx] = np.ma.zeros((n_seg),fill_value=fill_value)
        Segment_sigma_geo[gtx].mask = np.ones((n_seg),dtype=bool)
        comm.Allreduce(sendbuf=[Distributed_sigma_geo.data, MPI.DOUBLE], \
            recvbuf=[Segment_sigma_geo[gtx].data, MPI.DOUBLE], op=MPI.SUM)
        comm.Allreduce(sendbuf=[Distributed_sigma_geo.mask, MPI.BOOL], \
            recvbuf=[Segment_sigma_geo[gtx].mask, MPI.BOOL], op=MPI.LAND)
        Distributed_sigma_geo = None
        # wait for all distributed processes to finish for beam
        comm.Barrier()

        # set values for invalid segments to fill_value of each variable
        Segment_Height[gtx].data[Segment_Height[gtx].mask] = Segment_Height[gtx].fill_value
        Segment_Land_Ice[gtx].data[Segment_Land_Ice[gtx].mask] = Segment_Land_Ice[gtx].fill_value
        Segment_dH_along[gtx].data[Segment_dH_along[gtx].mask] = Segment_dH_along[gtx].fill_value
        Segment_dH_across[gtx].data[Segment_dH_across[gtx].mask] = Segment_dH_across[gtx].fill_value
        Segment_Height_Error[gtx].data[Segment_Height_Error[gtx].mask] = Segment_Height_Error[gtx].fill_value
        Segment_Land_Ice_Error[gtx].data[Segment_Land_Ice_Error[gtx].mask] = Segment_Land_Ice_Error[gtx].fill_value
        Segment_dH_along_Error[gtx].data[Segment_dH_along_Error[gtx].mask] = Segment_dH_along_Error[gtx].fill_value
        Segment_dH_across_Error[gtx].data[Segment_dH_across_Error[gtx].mask] = Segment_dH_across_Error[gtx].fill_value
        Segment_Mean_Median[gtx].data[Segment_Mean_Median[gtx].mask] = Segment_Mean_Median[gtx].fill_value
        Segment_X_atc[gtx].data[Segment_X_atc[gtx].mask] = Segment_X_atc[gtx].fill_value
        Segment_X_spread[gtx].data[Segment_X_spread[gtx].mask] = Segment_X_spread[gtx].fill_value
        Segment_Y_atc[gtx].data[Segment_Y_atc[gtx].mask] = Segment_Y_atc[gtx].fill_value
        Segment_sigma_geo[gtx].data[Segment_sigma_geo[gtx].mask] = Segment_sigma_geo[gtx].fill_value
        Segment_N_Fit[gtx].data[Segment_N_Fit[gtx].mask] = Segment_N_Fit[gtx].fill_value
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

        # save tep and dead time information and statistics
        IS2_atl03_fit['ancillary_data'][gtx] = {}
        IS2_atl03_attrs['ancillary_data'][gtx] = {}
        # tep time of day
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
        # tep window start
        IS2_atl03_fit['ancillary_data'][gtx]['tx_start'] = np.array(tep[gtx]['tx_start'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['long_name'] = "Start of the TEP Window"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_start']['description'] = ("Starting time for the "
            "window centered around the primary TEP arrival for calculating the transmit pulse shape.")
        # tep window end
        IS2_atl03_fit['ancillary_data'][gtx]['tx_end'] = np.array(tep[gtx]['tx_end'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['long_name'] = "End of the TEP Window"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_end']['description'] = ("Ending time for the "
            "window centered around the primary TEP arrival for calculating the transmit pulse shape.")
        # tep robust dispersion estimator
        IS2_atl03_fit['ancillary_data'][gtx]['tx_robust_sprd'] = np.array(tep[gtx]['tx_robust_sprd'])
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['long_name'] = "Robust Spread"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['tx_robust_sprd']['description'] = ("Temporal width of "
            "the transmit pulse (sec), calculated from the RDE of the primary TEP waveform")
        # tep full width at half maximum
        IS2_atl03_fit['ancillary_data'][gtx]['sigma_tx'] = np.array(tep[gtx]['sigma_tx'])
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['long_name'] = "Duration of Transmit Pulse"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['source'] = tep[gtx]['pce']
        IS2_atl03_attrs['ancillary_data'][gtx]['sigma_tx']['description'] = ("Temporal duration of "
            "the transmit pulse (sec), calculated from the FWHM of the TEP waveform")
        # mean dead time
        IS2_atl03_fit['ancillary_data'][gtx]['t_dead'] = np.array(mean_dead_time[gtx])
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead'] = {}
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['units'] = "seconds"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['long_name'] = "Dead-time"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['contentType'] = "auxiliaryInformation"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['source'] = "CAL42"
        IS2_atl03_attrs['ancillary_data'][gtx]['t_dead']['description'] = ("Mean dead-time for "
            "channels in the detector (sec)")

        # copy beam variables
        IS2_atl03_fit[gtx] = dict(land_ice_segments={})
        IS2_atl03_fill[gtx] = dict(land_ice_segments={})
        IS2_atl03_attrs[gtx] = dict(land_ice_segments={})
        # group attributes for beam
        IS2_atl03_attrs[gtx]['Description'] = fileID[gtx].attrs['Description']
        IS2_atl03_attrs[gtx]['atlas_pce'] = fileID[gtx].attrs['atlas_pce']
        IS2_atl03_attrs[gtx]['atlas_beam_type'] = fileID[gtx].attrs['atlas_beam_type']
        IS2_atl03_attrs[gtx]['groundtrack_id'] = fileID[gtx].attrs['groundtrack_id']
        IS2_atl03_attrs[gtx]['atmosphere_profile'] = fileID[gtx].attrs['atmosphere_profile']
        IS2_atl03_attrs[gtx]['atlas_spot_number'] = fileID[gtx].attrs['atlas_spot_number']
        IS2_atl03_attrs[gtx]['sc_orientation'] = fileID[gtx].attrs['sc_orientation']
        # group attributes for land_ice_segments
        IS2_atl03_attrs[gtx]['land_ice_segments']['Description'] = ("The land_ice_segments group "
            "contains the primary set of derived products. This includes geolocation, height, and "
            "standard error and quality measures for each segment. This group is sparse, meaning "
            "that parameters are provided only for pairs of segments for which at least one beam "
            "has a valid surface-height measurement.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['data_rate'] = ("Data within this group are "
            "sparse.  Data values are provided only for those ICESat-2 20m segments where at "
            "least one beam has a valid land ice height measurement.")

        # geolocation, time and segment ID
        # delta time
        IS2_atl03_fit[gtx]['land_ice_segments']['delta_time'] = Segment_delta_time[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['delta_time'] = None
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
        # latitude
        IS2_atl03_fit[gtx]['land_ice_segments']['latitude'] = Segment_Latitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['latitude'] = None
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
        # longitude
        IS2_atl03_fit[gtx]['land_ice_segments']['longitude'] = Segment_Longitude[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['longitude'] = None
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
        # segment ID
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
        # land ice height corrected for first photon bias and transmit-pulse shape
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
        # land ice height errors (max of fit or first photon bias uncertainties)
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
        # vertical geolocation error due to PPD and POD
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
        # segment quality summary
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

        # dem variables
        IS2_atl03_fit[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['Description'] = ("The dem group "
            "contains the reference digital elevation model and geoid heights.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['dem']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        # geoid height
        fv = fileID[gtx]['geophys_corr']['geoid'].attrs['_FillValue']
        geoid = np.ma.array(fileID[gtx]['geophys_corr']['geoid'][:], fill_value=fv)
        geoid_h = is2tk.fit.segment_mean(geoid)
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

        # geophysical variables
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['Description'] = ("The geophysical group "
            "contains parameters used to correct segment heights for geophysical effects, parameters "
            "related to solar background and parameters indicative of the presence or absence of clouds.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")

        # background rate
        bckgrd = is2tk.fit.segment_mean(Segment_Background[gtx])
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
        # blowing snow PSC flag
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
        # blowing snow confidence
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
        # blowing snow optical depth
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
        # cloud flag ASR
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
        # cloud flag atm
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
        # multiple scattering warning flag
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
        # range bias correction
        fv = fileID[gtx]['geolocation']['range_bias_corr'].attrs['_FillValue']
        range_bias_corr = np.ma.array(fileID[gtx]['geolocation']['range_bias_corr'][:], fill_value=fv)
        segment_range_bias = is2tk.fit.segment_mean(range_bias_corr)
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = segment_range_bias
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = segment_range_bias.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['long_name'] = "Range bias correction"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['description'] = "The range_bias estimated from geolocation analysis"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['range_bias_corr']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # total neutral atmosphere delay correction
        fv = fileID[gtx]['geolocation']['neutat_delay_total'].attrs['_FillValue']
        neutat_delay_total = np.ma.array(fileID[gtx]['geolocation']['neutat_delay_total'][:], fill_value=fv)
        segment_neutat_delay = is2tk.fit.segment_mean(neutat_delay_total)
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = segment_neutat_delay
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = segment_neutat_delay.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['long_name'] = "Total Neutral Atmospheric Delay"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['description'] = "Total neutral atmosphere delay correction (wet+dry)"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['neutat_delay_total']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # solar elevation
        fv = fileID[gtx]['geolocation']['solar_elevation'].attrs['_FillValue']
        solar_elevation = np.ma.array(fileID[gtx]['geolocation']['solar_elevation'][:], fill_value=fv)
        segment_solar_elevation = is2tk.fit.segment_mean(solar_elevation)
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
        # solar azimuth
        fv = fileID[gtx]['geolocation']['solar_azimuth'].attrs['_FillValue']
        solar_azimuth = np.ma.array(fileID[gtx]['geolocation']['solar_azimuth'][:], fill_value=fv)
        segment_solar_azimuth = is2tk.fit.segment_mean(solar_azimuth)
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

        # geophysical correction values at segment reference photons
        # dynamic atmospheric correction
        fv = fileID[gtx]['geophys_corr']['dac'].attrs['_FillValue']
        dac = np.ma.array(fileID[gtx]['geophys_corr']['dac'][:], fill_value=fv)
        segment_dac = is2tk.fit.segment_mean(dac)
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
        # solid earth tide
        fv = fileID[gtx]['geophys_corr']['tide_earth'].attrs['_FillValue']
        tide_earth = np.ma.array(fileID[gtx]['geophys_corr']['tide_earth'][:], fill_value=fv)
        segment_earth_tide = is2tk.fit.segment_mean(tide_earth)
        IS2_atl03_fit[gtx]['land_ice_segments']['geophysical']['tide_earth'] = segment_earth_tide
        IS2_atl03_fill[gtx]['land_ice_segments']['geophysical']['tide_earth'] = segment_earth_tide.fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['contentType'] = "referenceInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['long_name'] = "Earth Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['description'] = "Solid Earth Tide"
        IS2_atl03_attrs[gtx]['land_ice_segments']['geophysical']['tide_earth']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # load tide
        fv = fileID[gtx]['geophys_corr']['tide_load'].attrs['_FillValue']
        tide_load = np.ma.array(fileID[gtx]['geophys_corr']['tide_load'][:], fill_value=fv)
        segment_load_tide = is2tk.fit.segment_mean(tide_load)
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
        # ocean tide
        fv = fileID[gtx]['geophys_corr']['tide_ocean'].attrs['_FillValue']
        tide_ocean = np.ma.array(fileID[gtx]['geophys_corr']['tide_ocean'][:], fill_value=fv)
        segment_ocean_tide = is2tk.fit.segment_mean(tide_ocean)
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
        # ocean pole tide
        fv = fileID[gtx]['geophys_corr']['tide_oc_pole'].attrs['_FillValue']
        tide_oc_pole = np.ma.array(fileID[gtx]['geophys_corr']['tide_oc_pole'][:], fill_value=fv)
        segment_oc_pole_tide = is2tk.fit.segment_mean(tide_oc_pole)
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
        # pole tide
        fv = fileID[gtx]['geophys_corr']['tide_pole'].attrs['_FillValue']
        tide_pole = np.ma.array(fileID[gtx]['geophys_corr']['tide_pole'][:], fill_value=fv)
        segment_pole_tide = is2tk.fit.segment_mean(tide_pole)
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

        # bias correction variables
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['Description'] = ("The bias_correction group "
            "contains information about the estimated first-photon bias, and the transmit-pulse-shape bias.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        # mean first photon bias
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
        # mean first photon bias uncertainty
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
        # median first photon bias
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
        # median first photon bias correction
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
        # first photon bias corrected number of photons
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
        # CAL-19 first photon bias
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
        # mean transmit pulse shape correction
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
        # median transmit pulse shape correction
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
        # difference between the mean and median of fit residuals
        IS2_atl03_fit[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = Segment_Mean_Median[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = Segment_Mean_Median[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['long_name'] = "mean median residual"
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['description'] = ("Difference between "
            "uncorrected mean and median of linear fit residuals")
        IS2_atl03_attrs[gtx]['land_ice_segments']['bias_correction']['med_r_fit']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"

        # fit statistics variables
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['Description'] = ("The fit_statistics group "
            "contains a variety of parameters that might indicate the quality of the fitted segment data. Data in "
            "this group are sparse, with dimensions matching the land_ice_segments group.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        # segment fit heights
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = Segment_Height[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = Segment_Height[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['long_name'] = "Height Mean"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['description'] = ("Mean surface "
            "height, not corrected for first-photon bias or pulse truncation")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_mean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # segment fit along-track slopes
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
        # segment fit across-track slopes
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
        # segment fit height errors
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = Segment_Height_Error[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = Segment_Height_Error[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['long_name'] = "Height Error"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['description'] = ("Propagated height "
            "error due to PE-height sampling error for height from the along-track fit, not including geolocation-induced error")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['sigma_h_mean']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # segment fit across-track slope errors
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
        # segment fit along-track slope errors
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
        # number of photons in fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = Segment_N_Fit[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = Segment_N_Fit[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['long_name'] = "Number of Photons in Fit"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['description'] = ("Number of PEs used to "
            "determine mean surface height in the iterative surface fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_fit_photons']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # size of the window used in the fit
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = Segment_Window[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = Segment_Window[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['contentType'] = "physicalMeasurement"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['long_name'] = "Surface Window Width"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['description'] = ("Width of the surface "
            "window, top to bottom")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['w_surface_window_final']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # signal-to-noise ratio
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
        # segment photon signal-to-noise ratio from photon classifier
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph'] = Segment_Photon_SNR[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph'] = Segment_Photon_SNR[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph']['long_name'] = "Maximum SNR"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph']['description'] = ("Maximum "
            "signal-to-noise ratio from the photon event classifier used to normalize the photon weights")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['snr_norm_ph']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # robust dispersion estimator
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = Segment_RDE[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = Segment_RDE[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['units'] = "meters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['contentType'] = "qualityInformation"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['long_name'] = "Robust Spread"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['description'] = ("RDE of misfit "
            "between PE heights and the along-track segment fit")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['h_robust_sprd']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # number of iterations for fit
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
        # number of photon event clusters
        IS2_atl03_fit[gtx]['land_ice_segments']['fit_statistics']['n_clusters'] = Segment_Clusters[gtx]
        IS2_atl03_fill[gtx]['land_ice_segments']['fit_statistics']['n_clusters'] = Segment_Clusters[gtx].fill_value
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters']['units'] = "1"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters']['contentType'] = "modelResult"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters']['long_name'] = "Number of Estimated Clusters"
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters']['description'] = ("Number of clusters calculated "
            "using weighted density-based spatial clustering")
        IS2_atl03_attrs[gtx]['land_ice_segments']['fit_statistics']['n_clusters']['coordinates'] = \
            "../segment_id ../delta_time ../latitude ../longitude"
        # signal source selection
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
        # number potential segment pulses
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

        # ground track variables
        IS2_atl03_fit[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_fill[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track'] = {}
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['Description'] = ("The ground_track group "
            "contains parameters describing the GT and RGT for each land ice segment, as well as angular "
            "information about the beams.")
        IS2_atl03_attrs[gtx]['land_ice_segments']['ground_track']['data_rate'] = ("Data within this group "
            "are stored at the land_ice_segments segment rate.")
        # along-track X coordinates of segment fit
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
        # along-track Y coordinates of segment fit
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
        # along-track X coordinate spread of points used in segment fit
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
        # elevation
        fv = fileID[gtx]['geolocation']['ref_elev'].attrs['_FillValue']
        ref_elev = np.ma.array(fileID[gtx]['geolocation']['ref_elev'][:], fill_value=fv)
        segment_ref_elev = is2tk.fit.segment_mean(ref_elev)
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
        # azimuth
        fv = fileID[gtx]['geolocation']['ref_azimuth'].attrs['_FillValue']
        ref_azimuth = np.ma.array(fileID[gtx]['geolocation']['ref_azimuth'][:], fill_value=fv)
        segment_ref_azimuth = is2tk.fit.segment_mean(ref_azimuth)
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

    # parallel h5py I/O does not support compression filters at this time
    if (comm.rank == 0):
        # use default output file name and path
        if args.output:
            output_file=os.path.expanduser(args.output)
        else:
            fargs=(SUB,'ATL06',YY,MM,DD,HH,MN,SS,TRK,CYCL,GRAN,RL,VERS,AUX)
            file_format='{0}{1}_{2}{3}{4}{5}{6}{7}_{8}{9}{10}_{11}_{12}{13}.h5'
            output_file=os.path.join(ATL03_dir,file_format.format(*fargs))
        # write to HDF5 file
        HDF5_ATL03_write(IS2_atl03_fit, IS2_atl03_attrs, COMM=comm,
            INPUT=[args.ATL03,args.ATL09], FILL_VALUE=IS2_atl03_fill,
            CLOBBER=True, FILENAME=output_file)
        # change the permissions level to MODE
        os.chmod(output_file, args.mode)
    # close the input ATL03 file
    fileID.close()

# PURPOSE: read ICESat-2 ATL09 HDF5 data file for specific variables
def read_HDF5_ATL09(FILENAME, pfl, D, ATTRIBUTES=True, COMM=None):
    # Open the HDF5 file for reading
    fileID = h5py.File(FILENAME, 'r', driver='mpio', comm=COMM)
    logging.info(FILENAME) if (COMM.rank == 0) else None

    # allocate python dictionaries for ICESat-2 ATL09 variables and attributes
    IS2_atl09_mds = {}
    IS2_atl09_attrs = {}

    # read profile reported for the ATLAS strong beams within the file
    IS2_atl09_mds[pfl] = dict(high_rate={})
    # extract delta_time for mapping ATL09 atmospheric parameters to ATL03
    delta_time = fileID[pfl]['high_rate']['delta_time'][:]
    # Calibrated Attenuated Backscatter at 25 hz
    high_rate_keys = ['aclr_true','bsnow_con','bsnow_dens','bsnow_h',
        'bsnow_h_dens','bsnow_od','bsnow_psc','cloud_flag_asr','cloud_flag_atm',
        'cloud_fold_flag','column_od_asr','column_od_asr_qf','msw_flag',
        'snow_ice','solar_azimuth','solar_elevation','surf_refl_true']
    # number of output ATL03 segments
    n_seg = len(D)
    # parallel indices for filling variables
    ii = np.arange(COMM.rank,n_seg,COMM.size)
    # extract variables of interest and map to ATL03 segments
    for key in high_rate_keys:
        val = np.copy(fileID[pfl]['high_rate'][key][:])
        fint = scipy.interpolate.interp1d(delta_time, val,
            kind='nearest', fill_value='extrapolate')
        IS2_atl09_mds[pfl]['high_rate'][key] = np.zeros((n_seg),dtype=val.dtype)
        IS2_atl09_mds[pfl]['high_rate'][key][ii] = fint(D[ii]).astype(val.dtype)

    # Getting attributes of included variables
    if ATTRIBUTES:
        # Getting attributes of IS2_atl09_mds profile variables
        IS2_atl09_attrs[pfl] = dict(high_rate={})
        # Global Group Attributes
        for att_name,att_val in fileID[pfl].attrs.items():
            IS2_atl09_attrs[pfl][att_name] = att_val
        # Variable Attributes
        for key in high_rate_keys:
            IS2_atl09_attrs[pfl]['high_rate'][key] = {}
            for att_name,att_val in fileID[pfl]['high_rate'][key].attrs.items():
                IS2_atl09_attrs[pfl]['high_rate'][key][att_name] = att_val

    # Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl09_attrs[att_name] = att_val

    # Closing the HDF5 file
    fileID.close()
    # Return the datasets and variables
    return (IS2_atl09_mds,IS2_atl09_attrs)

# PURPOSE: outputting the reduced and corrected ICESat-2 data to HDF5
def HDF5_ATL03_write(IS2_atl03_data, IS2_atl03_attrs, COMM=None, INPUT=None,
    FILENAME='', FILL_VALUE=None, CLOBBER=True):
    # setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # open output HDF5 file
    fileID = h5py.File(FILENAME, clobber)#, driver='mpio', comm=COMM)
    logging.info(FILENAME) if (COMM.rank == 0) else None

    # create HDF5 records
    h5 = {}

    # # ICESat-2 spacecraft orientation at time
    # fileID.create_group('orbit_info')
    # h5['orbit_info'] = {}
    # for k,v in IS2_atl03_data['orbit_info'].items():
    #     # Defining the HDF5 dataset variables
    #     val = 'orbit_info/{0}'.format(k)
    #     h5['orbit_info'][k] = fileID.create_dataset(val, np.shape(v), data=v,
    #         dtype=v.dtype, compression='gzip')
    #     # add HDF5 variable attributes
    #     for att_name,att_val in IS2_atl03_attrs['orbit_info'][k].items():
    #         h5['orbit_info'][k].attrs[att_name] = att_val

    # information ancillary to the data product
    # number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    # and ATLAS Standard Data Product (SDP) epoch (2018-01-01T00:00:00Z UTC)
    h5['ancillary_data'] = {}
    for k in ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc','end_cycle',
        'end_geoseg','end_gpssow','end_gpsweek','end_orbit','end_region',
        'end_rgt','granule_end_utc','granule_start_utc','release','start_cycle',
        'start_geoseg','start_gpssow','start_gpsweek','start_orbit','start_region',
        'start_rgt','version']:
        # Defining the HDF5 dataset variables
        v = IS2_atl03_data['ancillary_data'][k]
        val = 'ancillary_data/{0}'.format(k)
        h5['ancillary_data'][k] = fileID.create_dataset(val, np.shape(v), data=v,
            dtype=v.dtype)
        # add HDF5 variable attributes
        for att_name,att_val in IS2_atl03_attrs['ancillary_data'][k].items():
            h5['ancillary_data'][k].attrs[att_name] = att_val

    # land_ice_segments variable groups for each beam
    GROUPS=['fit_statistics','geophysical','ground_track','dem','bias_correction']
    # write each output beam
    for gtx in ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']:
        fileID.create_group(gtx)
        fileID['ancillary_data'].create_group(gtx)

        # add HDF5 group attributes for beam
        for att_name in ['Description','atlas_pce','atlas_beam_type',
            'groundtrack_id','atmosphere_profile','atlas_spot_number',
            'sc_orientation']:
            fileID[gtx].attrs[att_name] = IS2_atl03_attrs[gtx][att_name]

        # add transmit pulse shape and dead time parameters
        h5['ancillary_data'][gtx] = {}
        for k,v in IS2_atl03_data['ancillary_data'][gtx].items():
            # attributes
            attrs = IS2_atl03_attrs['ancillary_data'][gtx][k]
            # Defining the HDF5 dataset variables
            val = 'ancillary_data/{0}/{1}'.format(gtx,k)
            h5['ancillary_data'][gtx][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype)
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5['ancillary_data'][gtx][k].attrs[att_name] = att_val

        # create land_ice_segments group
        fileID[gtx].create_group('land_ice_segments')
        h5[gtx] = dict(land_ice_segments={})
        for att_name in ['Description','data_rate']:
            att_val = IS2_atl03_attrs[gtx]['land_ice_segments'][att_name]
            fileID[gtx]['land_ice_segments'].attrs[att_name] = att_val

        # segment_id
        v = IS2_atl03_data[gtx]['land_ice_segments']['segment_id']
        attrs = IS2_atl03_attrs[gtx]['land_ice_segments']['segment_id']
        # # parallel indices for filling variables
        # pind = np.arange(COMM.rank,len(v),COMM.size)
        # Defining the HDF5 dataset variables
        val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments','segment_id')
        h5[gtx]['land_ice_segments']['segment_id'] = fileID.create_dataset(val,
            np.shape(v), data=v, dtype=v.dtype, compression='gzip')
        # with h5[gtx]['land_ice_segments']['segment_ID'].collective:
        #     h5[gtx]['land_ice_segments']['segment_id'][pind] = v[pind]
        # make dimension
        h5[gtx]['land_ice_segments']['segment_id'].make_scale('segment_id')
        # add HDF5 variable attributes
        for att_name,att_val in attrs.items():
            h5[gtx]['land_ice_segments']['segment_id'].attrs[att_name] = att_val

        # geolocation, time and height variables
        for k in ['latitude','longitude','delta_time','h_li','h_li_sigma',
            'sigma_geo_h','atl06_quality_summary']:
            # values and attributes
            v = IS2_atl03_data[gtx]['land_ice_segments'][k]
            attrs = IS2_atl03_attrs[gtx]['land_ice_segments'][k]
            fillvalue = FILL_VALUE[gtx]['land_ice_segments'][k]
            # Defining the HDF5 dataset variables
            val = '{0}/{1}/{2}'.format(gtx,'land_ice_segments',k)
            h5[gtx]['land_ice_segments'][k] = fileID.create_dataset(val,
                np.shape(v), data=v, dtype=v.dtype, fillvalue=fillvalue,
                compression='gzip')
            # with h5[gtx]['land_ice_segments'][k].collective:
            #     h5[gtx]['land_ice_segments'][k][pind] = v[pind]
            # attach dimensions
            for i,dim in enumerate(['segment_id']):
                h5[gtx]['land_ice_segments'][k].dims[i].attach_scale(
                    h5[gtx]['land_ice_segments'][dim])
            # add HDF5 variable attributes
            for att_name,att_val in attrs.items():
                h5[gtx]['land_ice_segments'][k].attrs[att_name] = att_val

        # fit statistics, geophysical corrections, geolocation and dem
        for key in GROUPS:
            fileID[gtx]['land_ice_segments'].create_group(key)
            h5[gtx]['land_ice_segments'][key] = {}
            for att_name in ['Description','data_rate']:
                att_val=IS2_atl03_attrs[gtx]['land_ice_segments'][key][att_name]
                fileID[gtx]['land_ice_segments'][key].attrs[att_name] = att_val
            for k,v in IS2_atl03_data[gtx]['land_ice_segments'][key].items():
                # attributes
                attrs = IS2_atl03_attrs[gtx]['land_ice_segments'][key][k]
                fillvalue = FILL_VALUE[gtx]['land_ice_segments'][key][k]
                # Defining the HDF5 dataset variables
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
                # attach dimensions
                for i,dim in enumerate(['segment_id']):
                    h5[gtx]['land_ice_segments'][key][k].dims[i].attach_scale(
                        h5[gtx]['land_ice_segments'][dim])
                # add HDF5 variable attributes
                for att_name,att_val in attrs.items():
                    h5[gtx]['land_ice_segments'][key][k].attrs[att_name] = att_val

    # HDF5 file title
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
    # add attribute for elevation instrument and designated processing level
    instrument = 'ATLAS > Advanced Topographic Laser Altimeter System'
    fileID.attrs['instrument'] = instrument
    fileID.attrs['source'] = 'Spacecraft'
    fileID.attrs['references'] = 'https://nsidc.org/data/icesat-2'
    fileID.attrs['processing_level'] = '4'
    # add attributes for input ATL03 and ATL09 files
    fileID.attrs['input_files'] = ','.join([os.path.basename(i) for i in INPUT])
    # find geospatial and temporal ranges
    lnmn,lnmx,ltmn,ltmx,tmn,tmx = (np.inf,-np.inf,np.inf,-np.inf,np.inf,-np.inf)
    for gtx in ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']:
        lon = IS2_atl03_data[gtx]['land_ice_segments']['longitude']
        lat = IS2_atl03_data[gtx]['land_ice_segments']['latitude']
        delta_time = IS2_atl03_data[gtx]['land_ice_segments']['delta_time']
        # setting the geospatial and temporal ranges
        lnmn = lon.min() if (lon.min() < lnmn) else lnmn
        lnmx = lon.max() if (lon.max() > lnmx) else lnmx
        ltmn = lat.min() if (lat.min() < ltmn) else ltmn
        ltmx = lat.max() if (lat.max() > ltmx) else ltmx
        tmn = delta_time.min() if (delta_time.min() < tmn) else tmn
        tmx = delta_time.max() if (delta_time.max() > tmx) else tmx
    # add geospatial and temporal attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    fileID.attrs['date_type'] = 'UTC'
    fileID.attrs['time_type'] = 'CCSDS UTC-A'
    # convert start and end time from ATLAS SDP seconds into UTC time
    time_utc = is2tk.convert_delta_time(np.array([tmn,tmx]))
    # convert to calendar date
    YY,MM,DD,HH,MN,SS = is2tk.time.convert_julian(time_utc['julian'],
        format='tuple')
    # add attributes with measurement date start, end and duration
    tcs = datetime.datetime(int(YY[0]), int(MM[0]), int(DD[0]),
        int(HH[0]), int(MN[0]), int(SS[0]), int(1e6*(SS[0] % 1)))
    fileID.attrs['time_coverage_start'] = tcs.isoformat()
    tce = datetime.datetime(int(YY[1]), int(MM[1]), int(DD[1]),
        int(HH[1]), int(MN[1]), int(SS[1]), int(1e6*(SS[1] % 1)))
    fileID.attrs['time_coverage_end'] = tce.isoformat()
    fileID.attrs['time_coverage_duration'] = f'{tmx-tmn:0.0f}'
    # add software information
    fileID.attrs['software_reference'] = is2tk.version.project_name
    fileID.attrs['software_version'] = is2tk.version.full_version
    fileID.attrs['software_revision'] = is2tk.utilities.get_git_revision_hash()
    # Closing the HDF5 file
    fileID.close()

# run main program
if __name__ == '__main__':
    main()
