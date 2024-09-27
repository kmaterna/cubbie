"""
The purpose of this script is to run some diagnostics on the mean coherence.
How did the coherence vary with other things, like perpendicular baseline and season?
What do histograms of coherence look like?
"""

import matplotlib.pyplot as plt
import datetime as dt
import os
import glob as glob
from subprocess import call, check_output
from ..intf_generating import sentinel_utilities


# ------------- DRIVERS ------------------ #

def analyze_coherence_function(igram_dir, raw_dir):
    [corr_file, baseline_table, corr_dirlist, plotname] = correlation_config(igram_dir, raw_dir)

    # Correlation vs. Other Things.
    if os.path.isfile(corr_file):
        print("Not computing file again")
    else:
        calc_write_corr_results(corr_dirlist, corr_file)  # Only need to do at the beginning
    [stem1, stem2, mean_corr] = sentinel_utilities.read_corr_results(corr_file)
    bl_tuples = sentinel_utilities.read_baseline_table(baseline_table)
    blt_datetimes = [x[1] for x in bl_tuples]
    xbaseline = [x[0] for x in bl_tuples]
    make_coh_vs_others_plots(stem1, stem2, mean_corr, blt_datetimes, xbaseline, plotname)
    return


# ------------ CORR vs OTHER STUFF -------------- # 

def correlation_config(igram_dir, raw_dir):
    corr_file = 'corr_results.txt'
    baseline_table = raw_dir + 'raw/baseline_table.dat'
    corr_dirlist = glob.glob(igram_dir + '???????_???????')
    plotname = "coherence_stats.eps"
    print("Config for analysis of coherence vs. other quantities: %s, %s, %s" % (corr_file, baseline_table, plotname))
    return [corr_file, baseline_table, corr_dirlist, plotname]


def calc_write_corr_results(dir_list, filename):
    """Write average coherence for each interferogram, in format: 2018005_2018017 0.152738928795 """
    ifile = open(filename, 'w')
    for item in dir_list:
        directname = item.split('/')[-1]  # format: 2015153_2015177
        call("gmt grdmath " + item + "/corr.grd MEAN = " + item + "/out.grd", shell=True)
        corr = check_output("gmt grdinfo " + item + "/out.grd | grep z | awk \'{print $3}\'", shell=True)
        corr = float(corr.split()[0])
        print("Computing for %s with mean coherence %f " % (directname, corr))
        call("rm " + item + "/out.grd", shell=True)
        ifile.write('%s %s\n' % (directname, corr))
    ifile.close()
    return


def make_coh_vs_others_plots(stem1, stem2, mean_corr, blt_datetimes, xbaseline, plotname):
    b_perp_baseline, temporal_baseline, season = [], [], []
    for i in range(len(mean_corr)):  # define a bunch of things for each interferogram.
        # How long in time?
        date1 = dt.datetime.strptime(sentinel_utilities.yj2ymd(stem1[i]), "%Y%m%d")
        date2 = dt.datetime.strptime(sentinel_utilities.yj2ymd(stem2[i]), "%Y%m%d")
        temporal_baseline.append(abs((date1 - date2).days))

        # Perpendicular Baseline?
        myindex1 = blt_datetimes.index(date1)
        myindex2 = blt_datetimes.index(date2)
        bl1 = xbaseline[myindex1]
        bl2 = xbaseline[myindex2]
        b_perp_baseline.append(abs(bl1 - bl2))

        # Which season?
        season.append(int(date1.strftime("%j")))

    summer0 = 160
    summer1 = 300

    f, axarr = plt.subplots(3, figsize=(15, 15))
    axarr[0].set_title('Mean Coherence of Sentinel-1 Interferograms', fontsize=24)
    axarr[0].plot(temporal_baseline, mean_corr, '.', markersize=13)
    axarr[0].set_xlabel('Time (days)', fontsize=20)
    axarr[0].set_ylabel('Mean Coherence', fontsize=20)
    axarr[0].tick_params(axis='both', labelsize=20)

    for i in range(len(season)):
        if summer0 < season[i] < summer1:
            if temporal_baseline[i] > 300:
                line1a, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '*', color='red', markersize=13,
                                        label='summer, 1+yr')
            else:
                line1, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '.', color='red', markersize=13,
                                       label='summer')
        else:
            if temporal_baseline[i] > 300:
                line2a, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '*', color='blue', markersize=13,
                                        label='not summer, 1+yr')
            else:
                line2, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '.', color='blue', markersize=13,
                                       label='not summer')
    axarr[1].set_xlabel('Baseline (m)', fontsize=20)
    axarr[1].set_ylabel('Mean Coherence', fontsize=20)
    axarr[1].tick_params(axis='both', labelsize=20)
    axarr[1].legend(handles=[line1, line2, line1a, line2a], loc=1, fontsize=16)

    for i in range(len(season)):
        if temporal_baseline[i] == 12:
            line1, = axarr[2].plot(season[i], mean_corr[i], '.', color='magenta', markersize=13, label='12 days')
        else:
            line2, = axarr[2].plot(season[i], mean_corr[i], '.', color='green', markersize=13, label='>12 days')
    axarr[2].plot([summer0, summer0], [0.05, 0.35], '--k')
    axarr[2].plot([summer1, summer1], [0.05, 0.35], '--k')
    axarr[2].set_xlabel('Day of Year', fontsize=20)
    axarr[2].set_ylabel('Mean Coherence', fontsize=20)
    axarr[2].legend(handles=[line1, line2], loc=2, fontsize=18)
    axarr[2].tick_params(axis='both', labelsize=20)
    print("Saving figure as %s " % plotname)
    plt.savefig(plotname)

    return
