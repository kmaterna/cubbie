# Implementing some "Exclude" and "Include" criteria for particular interferograms before doing TS analysis.
# The internal intf_tuple is: (d1, d2, intf_filename, corr_filename)

import sys
import datetime as dt
import matplotlib.pyplot as plt
from . import stacking_utilities


def exclude_intfs_manually(total_intf_tuple, skip_file):
    """
    :param total_intf_tuple: list of tuples (d1, d2, intf_filename, corr_filename)
    :param skip_file: string, filename. Use "" to skip.
    """
    print("Excluding intfs based on manual_exclude file %s." % skip_file);
    print(" Started with %d total interferograms. " % (len(total_intf_tuple)));
    select_intf_tuple, manual_removes = [], [];
    if skip_file == "":
        print(" No manual exclude file provided.\n Returning all %d interferograms. " % (len(total_intf_tuple)));
        select_intf_tuple = total_intf_tuple;
    else:
        print(" Excluding the following interferograms based on SkipFile %s: " % skip_file);
        ifile = open(skip_file, 'r');
        for line in ifile:
            manual_removes.append(line.split()[0]);
        ifile.close();
        print(manual_removes);

        if not manual_removes:
            select_intf_tuple = total_intf_tuple;
        else:
            # Checking to see if each interferogram should be included.
            for igram in total_intf_tuple:
                include_flag = 1;
                for scene in manual_removes:
                    if scene in igram[2]:
                        include_flag = 0;
                if include_flag == 1:
                    select_intf_tuple.append(igram);
        print(" Returning %d interferograms " % len(select_intf_tuple));
    return select_intf_tuple;


def include_only_coseismic_intfs(total_intf_tuple, coseismic):
    """
    Implements a filter for spanning a coseismic interval, if included.
    :param total_intf_tuple: list of tuples (d1, d2, intf_filename, corr_filename)
    :param coseismic: string representing coseismic epoch, YYYYMMDD. Use "" to skip.
    """
    select_intf_tuple = [];
    if coseismic == "":
        return total_intf_tuple;
    else:
        print("Returning only interferograms that cross coseismic event at %s " % (
            dt.datetime.strftime(coseismic, "%Y-%m-%d")))
        for mytuple in total_intf_tuple:
            if mytuple[0] < coseismic < mytuple[1]:
                select_intf_tuple.append(mytuple);  # in the case of a coseismic constraint
        print(" Returning %d interferograms " % len(select_intf_tuple));
        return select_intf_tuple;


def include_intfs_by_time_range(total_intf_tuple, start_time, end_time):
    """
    Look for each interferogram that falls totally within time range given in config file.
    :param total_intf_tuple: list of tuples (d1, d2, intf_filename, corr_filename)
    :param start_time: string, YYYY-MM-DD. Use "" to skip.
    :param end_time: string, YYYY-MM-DD. Use "" to skip.
    """
    if start_time == "" and end_time == "":
        return total_intf_tuple;
    print("Including only interferograms in time range %s to %s." % (dt.datetime.strftime(start_time, "%Y-%m-%d"),
                                                                     dt.datetime.strftime(end_time, "%Y-%m-%d")));
    print(" Starting with %d interferograms " % len(total_intf_tuple))
    select_intf_tuple = [];
    for mytuple in total_intf_tuple:
        if start_time <= mytuple[0] <= end_time:
            if start_time <= mytuple[1] <= end_time:
                select_intf_tuple.append(mytuple);  # in the case of no coseismic constraint
    print(" Returning %d interferograms " % len(select_intf_tuple));
    return select_intf_tuple;


def include_timeinterval_intfs(total_intf_tuple, intf_timespan):
    """
    Only include interferograms of a certain time interval (such as shorter than one year, or longer than one year);
    intf_timespan is a string with the format '300+' for longer than 300 days etc.
    Use "" to skip.
    """
    if intf_timespan == "":
        return total_intf_tuple;
    select_intf_tuple = [];
    if intf_timespan[-1] == '+':
        criterion = 'longer'
    else:
        criterion = 'shorter'
    days = int(intf_timespan[0:-1]);
    print("Only including interferograms %s than %d days " % (criterion, days));
    for mytuple in total_intf_tuple:
        datedelta = (mytuple[1] - mytuple[0]).days;
        if criterion == "longer":
            if datedelta > days:
                select_intf_tuple.append(mytuple);
        else:
            if datedelta < days:
                select_intf_tuple.append(mytuple);
    print(" Returning %d interferograms " % len(select_intf_tuple));
    return select_intf_tuple;


def write_intf_record(intf_tuple_list, record_file):
    """Write a record of interferogram files selected for analysis. """
    print("Writing out list of %d interferograms used in this run to %s" % (len(intf_tuple_list), record_file));
    ofile = open(record_file, 'w');
    ofile.write("List of %d interferograms used in this run:\n" % (len(intf_tuple_list)));
    for mytuple in intf_tuple_list:
        ofile.write("%s\n" % (mytuple[2]));
    ofile.close();
    return;


# Metrics and Plotting
def make_igram_stick_plot(intf_file_tuples, ts_output_dir):
    """
    :param intf_file_tuples: list of (d1, d2, filename, filename) tuples
    :param ts_output_dir: string
    """
    print("Making simple plot of interferograms used.")
    plt.figure(dpi=300, figsize=(8, 7));
    for i in range(len(intf_file_tuples)):
        plt.plot([intf_file_tuples[i][0], intf_file_tuples[i][1]], [i, i], '.', markersize=7, linestyle=None,
                 color='gray');
        plt.plot([intf_file_tuples[i][0], intf_file_tuples[i][1]], [i, i], markersize=5);
    plt.title(str(len(intf_file_tuples)) + ' Interferograms Used');
    plt.xlabel('Time');
    plt.ylabel('Interferogram Number');
    plt.savefig(ts_output_dir + "/intf_record.png");
    return;


# DRIVER FOR SELECTING INTERFEROGRAMS
def make_selection_of_intfs(config_params):
    """
    HERE IS WHERE YOU SELECT WHICH INTERFEROGRAMS YOU WILL BE USING.
    WE MIGHT APPLY A MANUAL EXCLUDE, OR A TIME CONSTRAINT, ETC IN CONFIG SETTINGS.
    The working internal intf_tuple is: (d1, d2, intf_filename, corr_filename)
    :returns: list of intf filenames, list of corr filenames
    """

    if config_params.ts_format == "velocities_from_timeseries":
        intf_files = stacking_utilities.get_list_of_ts_grids(config_params);
        return intf_files, [];

    intf_tuples = stacking_utilities.get_list_of_intf_all(config_params);

    # Use the config file to excluse certain time ranges and implement coseismic constraints
    select_intf_tuples = include_intfs_by_time_range(intf_tuples, config_params.start_time, config_params.end_time);
    select_intf_tuples = include_only_coseismic_intfs(select_intf_tuples, config_params.coseismic);

    # Do you only want to include long or short interferograms?
    select_intf_tuples = include_timeinterval_intfs(select_intf_tuples, config_params.intf_timespan);

    # Manual Excludes?
    select_intf_tuples = exclude_intfs_manually(select_intf_tuples, config_params.skip_file);

    # Writing the exact interferograms used in this run, and returning file names.
    write_intf_record(select_intf_tuples, config_params.ts_output_dir + "/intf_record.txt")
    make_igram_stick_plot(select_intf_tuples, config_params.ts_output_dir);  # always make stick plot
    select_intf_list = [mytuple[2] for mytuple in select_intf_tuples]
    select_corr_list = [mytuple[3] for mytuple in select_intf_tuples]
    if len(select_intf_tuples) == 0:
        print("Error! Not starting with any interferograms. Exiting.");
        sys.exit(0);
    return select_intf_list, select_corr_list;
