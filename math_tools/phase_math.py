import numpy as np
import math


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
