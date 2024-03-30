"""
Functions for reading/writing testing pixels
"""

import numpy as np
import collections
import datetime as dt

Igrams = collections.namedtuple("Igrams", ["dt1", "dt2", "juldays", "datestrs", "x_axis_days", "phase", "corr"])


def read_testing_pixel(ifile, coherence=True):
    print("Reading file %s " % ifile)
    if coherence:
        [juldays, phase, corr] = np.loadtxt(ifile, usecols=(1, 2, 3), unpack=True,
                                            dtype={'names': ('juldays', 'phase', 'corr'),
                                                   'formats': ('U15', np.float, np.float)})
    else:
        [juldays, phase] = np.loadtxt(ifile, usecols=(1, 2), unpack=True,
                                      dtype={'names': ('juldays', 'phase'), 'formats': ('U15', np.float)})
        corr = None
    day1 = [dt.datetime.strptime(x[0:7], "%Y%j") for x in juldays]
    day2 = [dt.datetime.strptime(x[8:], "%Y%j") for x in juldays]
    datestr_first = [x[0:7] for x in juldays]
    datestr_second = [x[8:] for x in juldays]
    datestrs = sorted(set(datestr_first + datestr_second))
    x_axis_datetimes = [dt.datetime.strptime(x, "%Y%j") for x in datestrs]
    x_axis_days = [(x - x_axis_datetimes[0]).days for x in
                   x_axis_datetimes]  # number of days since first acquisition.
    Test_Igrams = Igrams(dt1=day1, dt2=day2, juldays=juldays, datestrs=datestrs, x_axis_days=x_axis_days,
                         phase=phase, corr=corr)
    return Test_Igrams


def take_coherent_igrams(full_Igrams, corr_limit):
    # A function to only take the most coherent interferograms
    # This is usually unnecessary, since the NSBAS routine filters out nans anyway.
    new_day1, new_day2, new_juldays, new_phase, new_corr = [], [], [], [], []
    for i in range(len(full_Igrams.dt1)):
        if np.isnan(full_Igrams.corr[i]) or full_Igrams.corr[i] <= corr_limit:
            continue
        else:
            new_day1.append(full_Igrams.dt1[i])
            new_day2.append(full_Igrams.dt2[i])
            new_juldays.append(full_Igrams.juldays[i])
            new_phase.append(full_Igrams.phase[i])
            new_corr.append(full_Igrams.corr[i])
    new_Igrams = Igrams(dt1=new_day1, dt2=new_day2, juldays=new_juldays, datestrs=full_Igrams.datestrs,
                        x_axis_days=full_Igrams.x_axis_days, phase=new_phase, corr=new_corr)
    return new_Igrams


def write_testing_pixel(intf_tuple, pixel_value, coh_value, filename):
    # Outputting a specific pixel for using its values later in testing
    print("Writing %s " % filename)
    ofile = open(filename, 'w')
    for i in range(len(intf_tuple.filepaths)):
        if coh_value is not None:
            ofile.write("%s %s %.4f %.4f\n" % (intf_tuple.filepaths[i], intf_tuple.date_pairs_julian[i],
                                               pixel_value[i], coh_value[i]))
        else:
            ofile.write("%s %s %.4f\n" % (intf_tuple.filepaths[i], intf_tuple.date_pairs_julian[i], pixel_value[i]))
    ofile.close()
    return
