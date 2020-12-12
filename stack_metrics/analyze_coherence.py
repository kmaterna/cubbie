# The purpose of this script is to run some diagnostics on the mean coherence. 
# How did the coherence vary with other things, like perpendicular baseline and season? 
# What do histograms of coherence look like? 
# Run from the F1 directory, for example

import matplotlib.pyplot as plt
import datetime as dt
import sys
import glob as glob
from subprocess import call, check_output
from intf_generating import sentinel_utilities


# ------------- DRIVERS ------------------ #

def analyze_coherence_function():
    [corr_file, baseline_table, corr_dirlist, plotname] = correlation_config();

    # Correlation vs. Other Things.
    calc_write_corr_results(corr_dirlist, corr_file);  # ONLY NEED TO DO AT THE BEGINNING
    [stem1, stem2, mean_corr] = sentinel_utilities.read_corr_results(corr_file);
    [stems_blt, _, xbaseline, _] = sentinel_utilities.read_baseline_table(baseline_table);
    make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline, plotname);

    # # Histograms, in 12-panel figures, one small panel for each interferogram
    # [file_names]=configure_histograms();
    # [xdata,ydata,corr_all,date_pairs]=rwr.reader_simple_format(file_names);
    # stack_metrics_tools.all_gridded_histograms(corr_all,date_pairs);
    return;


# ------------ CORR vs OTHER STUFF -------------- # 

def correlation_config():
    corr_file = 'corr_results.txt';
    baseline_table = 'raw/baseline_table.dat';
    corr_dirlist = glob.glob('intf_all/2*');
    plotname = "coherence_stats.eps";
    print(
        "Config for analysis of coherence versus other quantities: %s, %s, %s" % (corr_file, baseline_table, plotname));
    return [corr_file, baseline_table, corr_dirlist, plotname];


def calc_write_corr_results(dir_list, filename):
    ifile = open(filename, 'w');
    for item in dir_list:
        directname = item.split('/')[-1];  # format: 2015153_2015177
        call("gmt grdmath " + item + "/corr.grd MEAN = " + item + "/out.grd", shell=True);
        corr = check_output("gmt grdinfo " + item + "/out.grd | grep z | awk \'{print $3}\'", shell=True);
        corr = float(corr.split()[0]);
        print("Computing for %s with mean coherence %f " % (directname, corr));
        SLCs = check_output("ls " + item + "/*.SLC", shell=True);
        SLCs = SLCs.decode('utf-8')
        slc1 = SLCs.split('\n')[0];
        slc1 = slc1.split('/')[-1];
        slc2 = SLCs.split()[1];
        slc2 = slc2.split('/')[-1];
        call("rm " + item + "/out.grd", shell=True);
        ifile.write('%s %s %s %s\n' % (directname, slc1, slc2, corr))
    # Format: 2018005_2018017 S1A20180106_ALL_F1.SLC S1A20180118_ALL_F1.SLC 0.152738928795
    ifile.close();
    return;


def make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline, plotname):
    b_perp_baseline = [];
    temporal_baseline = [];
    season = [];
    for i in range(len(mean_corr)):  # define a bunch of things for each interferogram.
        # How long in time?
        datestr1 = stem1[i].split('_')[1];  # stem1 has format "S1_20180615_ALL_F1"
        date1 = dt.datetime.strptime(datestr1, '%Y%m%d');
        datestr2 = stem2[i].split('_')[1];
        date2 = dt.datetime.strptime(datestr2, '%Y%m%d');
        temporal_baseline.append(abs((date1 - date2).days));

        # Perpendicular Baseline?
        myindex1 = stems_blt.index(stem1[i]);
        myindex2 = stems_blt.index(stem2[i]);
        bl1 = xbaseline[myindex1];
        bl2 = xbaseline[myindex2];
        b_perp_baseline.append(abs(bl1 - bl2));

        # Which season?
        season.append(int(date1.strftime("%j")));

    summer0 = 160;
    summer1 = 300;

    f, axarr = plt.subplots(3, figsize=(15, 15));
    axarr[0].set_title('Mean Coherence of Sentinel-1 Interferograms', fontsize=24)
    axarr[0].plot(temporal_baseline, mean_corr, '.', markersize=13);
    axarr[0].set_xlabel('Time (days)', fontsize=20)
    axarr[0].set_ylabel('Mean Coherence', fontsize=20)
    axarr[0].tick_params(axis='both', labelsize=20)

    for i in range(len(season)):
        if summer0 < season[i] < summer1:
            if temporal_baseline[i] > 300:
                line1a, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '*', color='red', markersize=13,
                                        label='summer, 1+yr');
            else:
                line1, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '.', color='red', markersize=13,
                                       label='summer');
        else:
            if temporal_baseline[i] > 300:
                line2a, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '*', color='blue', markersize=13,
                                        label='not summer, 1+yr');
            else:
                line2, = axarr[1].plot(b_perp_baseline[i], mean_corr[i], '.', color='blue', markersize=13,
                                       label='not summer');
    axarr[1].set_xlabel('Baseline (m)', fontsize=20)
    axarr[1].set_ylabel('Mean Coherence', fontsize=20)
    axarr[1].tick_params(axis='both', labelsize=20)
    axarr[1].legend(handles=[line1, line2, line1a, line2a], loc=1, fontsize=16);

    for i in range(len(season)):
        if temporal_baseline[i] == 12:
            line1, = axarr[2].plot(season[i], mean_corr[i], '.', color='magenta', markersize=13, label='12 days');
        else:
            line2, = axarr[2].plot(season[i], mean_corr[i], '.', color='green', markersize=13, label='>12 days');
    axarr[2].plot([summer0, summer0], [0.05, 0.35], '--k');
    axarr[2].plot([summer1, summer1], [0.05, 0.35], '--k');
    axarr[2].set_xlabel('Day of Year', fontsize=20)
    axarr[2].set_ylabel('Mean Coherence', fontsize=20)
    axarr[2].legend(handles=[line1, line2], loc=2, fontsize=18)
    axarr[2].tick_params(axis='both', labelsize=20)
    print("Saving figure as %s " % plotname);
    plt.savefig(plotname);

    return;


# ------------- HISTOGRAMS ------------ # 
# This makes a bunch of 12-panel figures with the coherence of each interferogram as a histogram
# Not particularly interesting if there's great coherence everywhere. 

def configure_histograms():
    file_dir = "intf_all";
    file_type = "corr.grd";

    file_names = glob.glob(file_dir + "/*/" + file_type);
    if len(file_names) == 0:
        print("Error! No files matching search pattern.");
        sys.exit(1);
    print("Reading " + str(len(file_names)) + " files.");
    return [file_names];


if __name__ == "__main__":
    analyze_coherence_function();
