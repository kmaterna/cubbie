import numpy as np
import math
import sys


def real_imag2phase_amp(real, imag):
    """
    Simple math function operating on 1-d arrays.

    :param real: 1-d array
    :param imag: 1-d array
    :returns: [1-d array phase, 1-d array amp]
    """
    phase, amp = [], []
    count = 0
    for x in range(len(imag)):
        phase.append(np.arctan2(imag[x], real[x]))
        amp.append(np.sqrt(real[x] * real[x] + imag[x] * imag[x]))
        if math.isnan(np.arctan2(imag[x], real[x])):
            count = count + 1
    if count > 0:
        print("Warning: %d of %d are nans" % (count, len(real)))
    else:
        print("Converted real,imag to phase,amp with no nans.")
    return [phase, amp]


def phase_amp2real_imag(phase, amp):
    """
    Simple math function operating on 1-d arrays.

    :param phase: 1-d array
    :param amp: 1-d array
    :returns: [1-d array real, 1-d array imag]
    """
    real, imag = [], []
    count = 0
    for i in range(len(phase)):
        real.append(amp[i] * np.cos(phase[i]))
        imag.append(amp[i] * np.sin(phase[i]))
        if math.isnan(amp[i] * np.sin(phase[i])):
            count = count + 1
    if count > 0:
        print("Warning: %d of %d are nans" % (count, len(real)))
    else:
        print("Converted phase,amp to real,imag with no nans.")
    return [real, imag]


def develop_mean_phase(phase_array):
    """
    Takes 1D array of phase values, and determines mean phase value (sensitive to cycle slips)
    Uses the "mean of circular quantities" technique.

    :param phase_array: 1-D array of phase values
    :returns: mean, float
    """
    xarray = [np.cos(i) for i in phase_array]
    yarray = [np.sin(i) for i in phase_array]
    xmean = np.mean(xarray)
    ymean = np.mean(yarray)
    tolerance = 0.00001
    if abs(xmean) < tolerance and abs(ymean) < tolerance:
        print("Error! The mean phase is undefined!")
        sys.exit(0)
    else:
        meanphase = np.arctan2(ymean, xmean)
    # Mean phase already wrapped into the -pi to pi range.
    return meanphase


def develop_median_phase(phase_array):
    """
    Takes 1D array of phase values, and determines median phase value (sensitive to cycle slips)
    Uses the "median of circular quantities" technique.

    :param phase_array: 1-D array of phase values
    :returns: median, float
    """
    xarray = [np.cos(i) for i in phase_array]
    yarray = [np.sin(i) for i in phase_array]
    xmedian = np.median(xarray)
    ymedian = np.median(yarray)
    tolerance = 0.00001
    if abs(xmedian) < tolerance and abs(ymedian) < tolerance:
        print("Error! The median phase is undefined!")
        sys.exit(0)
    else:
        medianphase = np.arctan2(ymedian, xmedian)
    return medianphase
