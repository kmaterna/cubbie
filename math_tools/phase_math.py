import numpy as np
import math
import sys


def real_imag2phase_amp(real, imag):
    # Simple trig math functions to convert real/imag to phase and amp.
    # Operates on 1D arrays.
    phase = [];
    amp = [];
    count = 0;
    for x in range(len(imag)):
        phase.append(np.arctan2(imag[x], real[x]));
        amp.append(np.sqrt(real[x] * real[x] + imag[x] * imag[x]));
        if math.isnan(np.arctan2(imag[x], real[x])):
            count = count + 1;
    if count > 0:
        print("Warning: %d of %d are nans" % (count, len(real)));
    else:
        print("Converted real,imag to phase,amp with no nans.");

    return [phase, amp];


def phase_amp2real_imag(phase, amp):
    # Math functions that operate on 1D arrays.
    real = [];
    imag = [];
    count = 0;
    for i in range(len(phase)):
        real.append(amp[i] * np.cos(phase[i]));
        imag.append(amp[i] * np.sin(phase[i]));
        if math.isnan(amp[i] * np.sin(phase[i])):
            count = count + 1
    if count > 0:
        print("Warning: %d of %d are nans" % (count, len(real)));
    else:
        print("Converted phase,amp to real,imag with no nans.");
    return [real, imag];


def develop_mean_phase(phase_array):
    # This function takes a 1D array of phase values, and determines the mean phase value (sensitive to cycle slips)
    # It uses the "mean of circular quantities" technique.
    xarray = [np.cos(i) for i in phase_array];
    yarray = [np.sin(i) for i in phase_array];
    xmean = np.mean(xarray);
    ymean = np.mean(yarray);
    tolerance = 0.00001;
    if abs(xmean) < tolerance and abs(ymean) < tolerance:
        print("Error! The mean phase is undefined!");
        sys.exit(0);
    else:
        meanphase = np.arctan2(ymean, xmean);
    # Mean phase already wrapped into the -pi to pi range.
    return meanphase;


def develop_median_phase(phase_array):
    # This function takes a 1D array of phase values, and determines the median phase value (sensitive to cycle slips)
    # It uses the "median of circular quantities" technique.
    xarray = [np.cos(i) for i in phase_array];
    yarray = [np.sin(i) for i in phase_array];
    xmedian = np.median(xarray);
    ymedian = np.median(yarray);
    tolerance = 0.00001;
    if abs(xmedian) < tolerance and abs(ymedian) < tolerance:
        print("Error! The median phase is undefined!");
        sys.exit(0);
    else:
        medianphase = np.arctan2(ymedian, xmedian);
    return medianphase;
