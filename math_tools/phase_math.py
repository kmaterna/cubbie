import numpy as np
import sys


def real_imag2phase_amp(real, imag):
    """
    Simple math function operating on arrays.

    :param real: np.array
    :param imag: np.array
    :returns: [np.array phase, np.array amp]
    """
    phase = np.arctan2(imag, real)
    amp = np.sqrt(np.square(real) + np.square(imag))
    return phase, amp


def phase_amp2real_imag(phase, amp):
    """
    Simple math function operating on numpy arrays.

    :param phase: np.array
    :param amp: np.array
    :returns: np.array real, np.array imag
    """
    real = amp * np.cos(phase)
    imag = amp * np.sin(phase)
    return real, imag


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
