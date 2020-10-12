# October 2020, testing pixels

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import collections
import nsbas
import sentinel_utilities


Igrams = collections.namedtuple("Igrams", ["dt1", "dt2", "juldays", "datestrs", "x_axis_days", "phase", "corr"]);


def read_test_pixel(ifile):
    [juldays, phase, corr] = np.loadtxt(ifile, usecols=(1, 2, 3), unpack=True,
                                        dtype={'names': ('juldays', 'phase', 'corr'),
                                               'formats': ('U15', np.float, np.float)});
    day1 = [dt.datetime.strptime(x[0:7], "%Y%j") for x in juldays];
    day2 = [dt.datetime.strptime(x[8:], "%Y%j") for x in juldays];
    datestr_first = [x[0:7] for x in juldays];
    datestr_second = [x[8:] for x in juldays];
    datestrs = sorted(set(datestr_first + datestr_second));
    x_axis_datetimes = [dt.datetime.strptime(x, "%Y%j") for x in datestrs];
    x_axis_days = [(x - x_axis_datetimes[0]).days for x in
                   x_axis_datetimes];  # number of days since first acquisition.
    Test_Igrams = Igrams(dt1=day1, dt2=day2, juldays=juldays, datestrs=datestrs, x_axis_days=x_axis_days,
                         phase=phase, corr=corr);
    return Test_Igrams;


def take_coherent_igrams(full_Igrams, corr_limit):
    # A function to only take the most coherent interferograms
    # This is usually unnecessary, since the NSBAS routine filters out nans anyway.
    new_day1, new_day2, new_juldays, new_phase, new_corr = [], [], [], [], [];
    for i in range(len(full_Igrams.dt1)):
        if np.isnan(full_Igrams.corr[i]) or full_Igrams.corr[i] <= corr_limit:
            continue;
        else:
            new_day1.append(full_Igrams.dt1[i])
            new_day2.append(full_Igrams.dt2[i])
            new_juldays.append(full_Igrams.juldays[i])
            new_phase.append(full_Igrams.phase[i])
            new_corr.append(full_Igrams.corr[i])
    new_Igrams = Igrams(dt1=new_day1, dt2=new_day2, juldays=new_juldays, datestrs=full_Igrams.datestrs,
                        x_axis_days=full_Igrams.x_axis_days, phase=new_phase, corr=new_corr);
    return new_Igrams;


def outputs(x_axis_days, ts):
    plt.figure()
    plt.plot(x_axis_days, ts, '.');
    plt.savefig('test.png');
    return;


if __name__ == "__main__":
    ifile = 'Testing_Data/testing_pixel_3.txt';
    baseline_table = 'Testing_Data/baseline_table.dat'
    [stems, times, baselines, missiondays] = sentinel_utilities.read_baseline_table(baseline_table);
    full_Igrams = read_test_pixel(ifile);
    full_Igrams = take_coherent_igrams(full_Igrams, 0.385)
    sentinel_utilities.make_network_plot(full_Igrams.juldays, stems, times, baselines, "pixel_baseline_plot.png");
    ts = nsbas.do_nsbas_pixel(full_Igrams.phase, full_Igrams.juldays, 0, 56, full_Igrams.datestrs, coh_value=full_Igrams.corr);
    outputs(full_Igrams.x_axis_days, ts);
