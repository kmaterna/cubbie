#!/usr/bin/python
# Perform a DEM error correction from Fattahi and Amelung, 2013, IEEE Proceedings.
# This correction operates in the time domain after the SBAS inversion.
# It solves for a term proportional to baseline history.

import numpy as np
import datetime as dt
import sentinel_utilities

def driver(ts_vector, datestrs, baseline_file):
    # A function to implement Fattahi and Amelung's 2013 paper
    # Right now, this assumes a linear velocity, although more complicated time histories can be implemented.
    # stems format: 'S1_20190105_ALL_F2'
    # times format: float years
    # xbaselines format: meters (first one 0 by definition)
    # datestrs format: '2015134'

    [_, times, baselines, _] = sentinel_utilities.read_baseline_table(baseline_file);
    dtarray = [];
    for i in range(len(times)):
        dtarray.append(dt.datetime.strptime(str(int(times[i] + 1)), '%Y%j'));

    # Re-order times and baselines in chronological order
    baselines = [x for _, x in sorted(zip(dtarray, baselines))];
    dtarray = sorted(dtarray);

    # design matrix: phase(t) = v(t-t0) + other terms + .... (4pi/lamda B(ti)/rsin(theta) z_error)
    # Baseline history: Bdot(i) = B(t_i)-B(t_i-t_i-1) / (t_i-t_i-1), i=[1-N]
    # velocity history: v(i) = phi(t_i)-phi(t_i-1)  / (t_i-t_i-1), i=[1-N]
    # The model we're solving for is [velocity, (4pi/lamda z_error/rsin(theta))].
    # I call that constant K_z_error
    G = np.ones((len(datestrs)-1, 2))
    v = np.zeros((len(datestrs)-1));
    Bdot = np.zeros((len(datestrs)-1));
    for i in range(0,len(datestrs)-1):
        v[i] = (ts_vector[i+1] - ts_vector[i]) / (dtarray[i+1]-dtarray[i]).days / 365.24;
        Bdot[i] = (baselines[i+1] - baselines[i]) / (dtarray[i+1]-dtarray[i]).days / 365.24;
        G[i,1] = Bdot[i];

    model = np.linalg.lstsq(G, v, rcond=0.001);  # rcond helps the solution converge
    K_z_error = model[0][1]  # constant time z_error;

    topo_phase = [K_z_error * (x-baselines[0]) for x in baselines];

    corrected_ts_vector = np.subtract(ts_vector,topo_phase);
    return corrected_ts_vector;