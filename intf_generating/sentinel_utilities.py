# Sentinel Utilities

import subprocess
import sys
import glob, os, re
import datetime as dt
from subprocess import check_output
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
matplotlib.use('Agg')

"""
Note: Baseline_tuple_list has format like: [(150.2, dt, 2015230, S1_20150514_ALL_F2), (etc)...]
Baseline_tuple_list: [(baseline, dt.dt, datedoystr, stem),]
Note: intf_list has format like: 'S1A20150310_ALL_F1:S1A20150403_ALL_F1'
"""


def get_all_xml_tiff_names(directory, polarization, swath, filetype='xml'):
    """Returns a matching list of xml or tiff filenames, and datestrs(yyyymmdd)"""
    pathname1 = directory + "/*-"+polarization+"-*-00" + swath + "." + filetype;
    pathname2 = directory + "/*-"+polarization+"-*-00" + str(int(swath) + 3) + "." + filetype;
    list_of_images_temp = glob.glob(pathname1) + glob.glob(pathname2);
    list_of_images = []
    list_of_datestrs = [];
    for item in list_of_images_temp:
        list_of_images.append(item[:])
        list_of_datestrs.append(get_datestr_from_xml(item[:]));
    return list_of_images, list_of_datestrs;

def get_all_safes_in_dir(directory):
    """A directory that contains a bunch of safe files.
    Give us the list of files and the datetimes they go with (matching lengths). """
    dirlist = glob.glob(directory + '/*.SAFE');
    datelist = [];
    for item in dirlist:
        datelist.append(safe_to_date(item));
    return dirlist, datelist;

def safe_to_date(filename):
    temp = filename.split('/')[-1];
    datestr = temp[17:25];
    dtobj = dt.datetime.strptime(datestr, "%Y%m%d");
    return dtobj;

def get_safes_of_date(dirname, one_date):
    """Return all safes that happened on a particular datetime """
    dirlist, datelist = get_all_safes_in_dir(dirname);
    retval = [];
    for item, date in zip(dirlist, datelist):
        if date == one_date:
            retval.append(item);
    retval = sorted(retval);
    return retval;

def get_SAFE_list_for_raw_orig(config_params):
    """Return the files we're going to put into raw_orig
    Location depends on whether we've made frames or not."""
    if config_params.frame1 != '':
        file_list, dt_list = get_all_safes_in_dir(config_params.FRAMES_dir);
        # if we're assembling frames, we use the FRAMES directory.
    else:
        file_list, dt_list = get_all_safes_in_dir(config_params.DATA_dir);
        # if we're not assembling frames, we use the DATA directory.
    return file_list, dt_list;

def get_previous_and_following_day(mydate):
    """ This is a function that takes a dt object and generates
    [dt1, dt2]: the day before and the day after the date in question. """
    tomorrow = mydate + dt.timedelta(days=1);
    yesterday = mydate - dt.timedelta(days=1);
    return [yesterday, tomorrow];

def get_datestr_from_xml(xml_name):
    """
    xml file has name like s1a-iw1-slc-vv-20150121t134413-20150121t134424-004270-005317-001.xml
    We want to return 20150121. 
    """
    xml_name = xml_name.split('/')[-1];
    mydate = xml_name[15:23];
    return mydate;

def get_sat_from_xml(xml_name):
    xml_name = xml_name.split('/')[-1];
    sat = xml_name[0:3];
    return sat;

def ymd2yj(ymd):
    # Turn something like "20150524" into "2015100", useful for file naming conventions
    tdate = dt.datetime.strptime(ymd, "%Y%m%d");
    tdate = tdate - dt.timedelta(days=1);
    return dt.datetime.strftime(tdate, "%Y%j");

def yj2ymd(yj):
    yj = int(yj) + 1;
    tdate = dt.datetime.strptime(str(yj), "%Y%j");
    return dt.datetime.strftime(tdate, "%Y%m%d");

def yj_to_prm_name(yj, swath):
    # takes string like "2019166" and converts to "S1_20190516_ALL_F1"
    ymd = yj2ymd(yj);
    prm_name = "S1_" + ymd + "_ALL_F"+swath;
    return prm_name;

def get_eof_from_date_sat(mydate, sat, eof_dir):
    """ This returns something like S1A_OPER_AUX_POEORB_OPOD_20160930T122957_V20160909T225943_20160911T005943.EOF.
        eof_dir can be relative or absolute
        It takes something like 20171204, s1a, eof_dir
        It can also take something like dt.datetime, s1a, eof_dir
    """
    if type(mydate) == str:   # if you pass a string, convert to datetime
        mydate = dt.datetime.strptime(mydate, "%Y%m%d");
    [previous_day, following_day] = get_previous_and_following_day(mydate);
    previous_day = dt.datetime.strftime(previous_day, "%Y%m%d");
    following_day = dt.datetime.strftime(following_day, "%Y%m%d");
    eof_name = glob.glob(eof_dir + "/" + sat.upper() + "*" + previous_day + "*" + following_day + "*.EOF");
    if not eof_name:
        print("ERROR: did not find any EOF files matching the pattern " + eof_dir + "/" + sat.upper() +
              "*" + previous_day + "*" + following_day + "*.EOF");
        print("Exiting...")
        sys.exit(1);
    else:
        eof_name = eof_name[0];
    return eof_name;


def glob_intf_computed(parent_dir):
    # Return format "2016332_2017308" (julian days)
    full_names = glob.glob(parent_dir+"/intf_all/*");
    intf_computed = [];
    for item in full_names:
        intf_computed.append(item.split('/')[-1]);
    return intf_computed;


def get_common_intfs(desired_swaths):
    # After each swath has made intfs, which interferograms do we have in all desired swaths?
    # Returns in format '2019082_2019106'
    print("Checking the progress of common intfs after intf formation. ")
    intfs_array, all_intfs = [], [];
    for item in desired_swaths:
        new_intfs = glob_intf_computed('F'+item);
        intfs_array.append(new_intfs);
        print("------- F"+item+" ------- %d intfs" % len(new_intfs))
        all_intfs = list(set(new_intfs+all_intfs))  # growing list of total interferograms from all swaths
    print("Total unique interferograms: %d " % len(all_intfs));
    common_intfs = [];
    for item in all_intfs:
        if item in intfs_array[0]:
            if len(intfs_array) == 2 and item in intfs_array[1]:
                common_intfs.append(item);
            else:
                if len(intfs_array) == 3 and item in intfs_array[1] and item in intfs_array[2]:
                    common_intfs.append(item);
    print("Total common interferograms: %d " % len(common_intfs))
    return common_intfs;


def compare_frames_with_safes(config_params):
    # checking for output consistency after frames are created.
    print("Checking the frames and acquisitions and see if all your data has been included. ");
    dirlist, datelist = get_all_safes_in_dir(config_params.DATA_dir);
    frames_dirlist, frames_datelist = get_all_safes_in_dir(config_params.FRAMES_dir);
    print("BEGIN: %d SAFES from %d unique dates." % (len(dirlist), len(set(datelist))));
    print("END: %d SAFES from %d unique dates." % (len(frames_dirlist), len(set(frames_datelist))));
    return;

def format_image_name_as_datestr(master_image):
    # Takes something like 'S1_20190616_ALL_F1' or 's1b-iw1-slc-vv-20190616t014910-20190616t015028-016715-01f75e-004'
    # Returns YYYYMMDD
    master_datestr = re.findall(r"\d\d\d\d\d\d\d\d", master_image);  # example: ['20150607']
    master_datestr = master_datestr[0];
    return master_datestr;


def write_data_in(polarization, swath, master_date="", target_dir="F1"):
    """
    data.in is a reference table that links the xml file with the correct orbit file.
    swath is str
    master_date has format S1_20170204_ALL or 20170204
    """
    list_of_images, list_of_datestrs = get_all_xml_tiff_names("F" + swath + "/raw_orig", polarization, swath);
    outfile = open(target_dir + "/data.in", 'w');
    if master_date == "":
        print("No master date selected by the user. Printing in random order.")
        for item, mydate in zip(list_of_images, list_of_datestrs):
            item = item.split("/")[-1];  # getting rid of the directory
            sat = get_sat_from_xml(item);
            eof_name = get_eof_from_date_sat(mydate, sat, "F" + swath + "/raw_orig");
            outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
    else:
        # write the master date first.
        master_date = format_image_name_as_datestr(master_date);
        for item, mydate in zip(list_of_images, list_of_datestrs):
            if mydate == master_date:
                print("Found master date %s. Putting it on first line." % master_date);
                item = item.split("/")[-1];  # getting rid of the directory
                sat = get_sat_from_xml(item);
                eof_name = get_eof_from_date_sat(mydate, sat, "F" + swath + "/raw_orig");
                outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
                break;
        # then write the other dates. 
        for item, mydate in zip(list_of_images, list_of_datestrs):
            if mydate != master_date:
                item = item.split("/")[-1];  # getting rid of the directory
                sat = get_sat_from_xml(item);
                eof_name = get_eof_from_date_sat(mydate, sat, "F" + swath + "/raw_orig");
                outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
    outfile.close();
    print(target_dir+"/data.in successfully printed.")
    return;


def read_baseline_table(baselinefilename):
    # Returns a chronologically sorted list of tuples of (baseline, datetime, datedoystr, stem) values
    # Example: [(150.2, dt, 2015230, S1_20150514_ALL_F2), (etc)...]
    if baselinefilename == '':
        print("Error! No baseline file provided. Exiting...");
        sys.exit(0);
    baselineFile = np.genfromtxt(baselinefilename, dtype=str);
    stems = baselineFile[:, 0].astype(str);
    if len(stems[0]) > 60:  # adapting for a different format of baseline_table.dat sometimes happens.
        new_stems = [];
        for item in stems:
            swath = int(item[-1]);
            if swath >= 4:
                swath = swath - 3;
            new_stems.append("S1_" + item[15:23] + "_ALL_F" + str(swath));
        stems = new_stems;
    times = baselineFile[:, 1].astype(float);
    missiondays = baselineFile[:, 2].astype(float);
    baselines = baselineFile[:, 4].astype(float);

    dtarray, datestrs = [], [];
    for i in range(len(times)):
        dtarray.append(dt.datetime.strptime(str(int(times[i] + 1)), '%Y%j'));
        datestrs.append(str(times[i] + 1)[0:7]);   # string with format "2014361"

    # Re-order times and baselines in chronological order, and re-package
    baselines = [x for _, x in sorted(zip(dtarray, baselines))];
    datestrs = [x for _, x in sorted(zip(dtarray, datestrs))];
    stems = [x for _, x in sorted(zip(dtarray, stems))];
    missiondays = [x for _, x in sorted(zip(dtarray, missiondays))];
    dtarray = sorted(dtarray);
    baseline_tuple_list = [];
    for i in range(len(baselines)):
        baseline_tuple_list.append((baselines[i], dtarray[i], datestrs[i], stems[i], missiondays[i]));

    return baseline_tuple_list;


def read_intf_table(tablefilename):
    tablefile = np.genfromtxt(tablefilename, dtype=str);
    intf_all = tablefile[:].astype(str);
    return intf_all;


def write_intf_table(intf_all, tablefilename):
    ofile = open(tablefilename, 'w');
    for i in intf_all:
        ofile.write("%s\n" % i);
    ofile.close();
    return;


# after running the baseline calculation from the first pre_proc_batch,
# choose a new master that is close to the median baseline and timespan.
def choose_master_image(baseline_table_file):
    print("Selecting master image; default_master is none");

    # load baseline table
    baseline_tuple_list = read_baseline_table(baseline_table_file);
    time = [x[4] for x in baseline_tuple_list];  # in missiondays (float)
    baseline = [x[0] for x in baseline_tuple_list];  # meters
    shortform_names = [x[3] for x in baseline_tuple_list];  # format S1_20170304_F1_ALL

    # calculate shortest distance from median to scenes
    consider_time = True
    if consider_time:
        time_baseline_scale = 1  # arbitrary scaling factor, units of (meters/day)
        sceneDistance = np.sqrt(
            ((time - np.median(time)) / time_baseline_scale) ** 2 + (baseline - np.median(baseline)) ** 2)
    else:
        sceneDistance = np.sqrt((baseline - np.median(baseline)) ** 2)

    minID = np.argmin(sceneDistance);   # sometimes this line is finicky, whether it needs [0] or not.
    master_shortform = shortform_names[minID];
    return master_shortform


def write_super_master_batch_config(masterid):
    ifile = open('batch.config', 'r');
    ofile = open('batch.config.new', 'w');
    for line in ifile:
        if 'master_image' in line:
            ofile.write('master_image = ' + masterid + '\n');
        else:
            ofile.write(line);
    ifile.close();
    ofile.close();
    subprocess.call(['mv', 'batch.config.new', 'batch.config'], shell=False);
    print("Writing master_image into batch.config");
    return;


def set_up_merge_unwrap(desired_swaths):
    print("Setting up merged unwrapping for swaths:");
    print(desired_swaths);
    outdir = 'merged'
    subprocess.call(["mkdir", "-p", outdir], shell=False);
    subprocess.call(["cp", "F"+desired_swaths[0]+"/topo/dem.grd", outdir], shell=False);  # need copy, not soft link.
    subprocess.call(["cp", "batch.config", outdir], shell=False);
    intf_all = get_common_intfs(desired_swaths);
    check_intf_all_sanity(desired_swaths, intf_all, 'phase.grd');  # defensive programming
    check_intf_all_sanity(desired_swaths, intf_all, 'corr.grd');  # defensive programming
    check_intf_all_sanity(desired_swaths, intf_all, 'mask.grd');  # defensive programming
    return intf_all, outdir;


def write_merge_batch_input(intf_all, master_image, desired_swaths=("1", "2", "3")):
    """
    # Necessary for merge swaths step.
    # The master of the first line should be the super master.
    """
    outfile = 'merged/inputfile.txt'
    print("\nWriting file " + outfile);
    master_image_in_format = ymd2yj(master_image[3:11])  # putting master image in the format taken by gmtsar dirs
    # Find the super master and stick it to the front.
    for item in intf_all:
        if master_image_in_format == item[0:7]:
            print("Found master image in intf %s" % item);
            print("Moving this image to the front");
            intf_all.insert(0, intf_all.pop(intf_all.index(item)));
            break;

    ofile = open(outfile, "w");
    for item in intf_all:
        first_part = yj_to_prm_name(item.split('_')[0], '1');  # first date into PRM format, swath 1
        second_part = yj_to_prm_name(item.split('_')[1], '1');  # second date into PRM format, swath 1
        if "1" in desired_swaths:
            ofile.write("../F1/intf_all/%s/:" % item);
            ofile.write(first_part + ".PRM:" + second_part + ".PRM,");
        if "2" in desired_swaths:
            ofile.write("../F2/intf_all/%s/:" % item);
            ofile.write(first_part.replace("_F1", "_F2") + ".PRM:" + second_part.replace("_F1", "_F2") + ".PRM");
        if "3" in desired_swaths:
            ofile.write(',');
            ofile.write("../F3/intf_all/%s/:" % item);
            ofile.write(first_part.replace("_F1", "_F3") + ".PRM:" + second_part.replace("_F1", "_F3") + ".PRM");
        ofile.write('\n');
    ofile.close();
    print("Done writing " + outfile);
    return;


def write_merge_unwrap(outfile):
    print("Writing unwrap instructions in file %s " % outfile)
    ofile = open(outfile, 'w');
    ofile.write('#!/bin/bash\n');
    ofile.write("cd merged\n");
    ofile.write("pwd\n");
    ofile.write("ls\n");
    ofile.write("merge_batch.csh inputfile.txt batch.config\n")
    ofile.write("cd ../\n");
    ofile.close();
    return;


def merge_wrapped(desired_swaths, master_image):
    """Merge desired swaths; don't unwrap or geocode."""
    if len(desired_swaths) == 1:  # for single swath
        igram_directory = "F"+str(desired_swaths[0])+"/intf_all/";
    else:
        common_intfs, igram_directory = set_up_merge_unwrap(desired_swaths);  # check directory + paths
        write_merge_batch_input(common_intfs, master_image, desired_swaths);  # put the master image in front
        subprocess.call(['merge_swaths_mod.csh', "inputfile.txt", "batch.config"], shell=False);  # from procdir
    return igram_directory;


def write_ordered_unwrapping(numproc, swath, sh_file, config_file):
    [stem1, stem2, mean_corr] = read_corr_results("corr_results.txt");

    stem1_ordered = [x for y, x in sorted(zip(mean_corr, stem1), reverse=True)];
    stem2_ordered = [x for y, x in sorted(zip(mean_corr, stem2), reverse=True)];

    outfile = open(sh_file, 'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch unwrap Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("cd F" + swath + "\n");
    outfile.write("rm intf?.in\n");
    for i, item in enumerate(stem1_ordered):
        outfile.write(
            'echo "' + stem1_ordered[i] + ":" + stem2_ordered[i] + '" >> intf' + str(np.mod(i, numproc)) + '.in\n');
    outfile.write("\n# Unwrap the interferograms.\n\n")
    outfile.write("ls intf?.in | parallel --eta 'unwrap_mod.csh {} " + config_file + "'\n\n\n");
    outfile.close();

    return;


def write_unordered_unwrapping(numproc, swath, sh_file, config_file):
    infile = 'F' + str(swath) + '/intf_record.in';
    intfs = [];
    for line in open(infile):
        intfs.append(line[0:-1]);
    outfile = open(sh_file, 'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch unwrap Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("cd F" + swath + "\n");
    outfile.write("rm intf?.in\n");
    for i, item in enumerate(intfs):
        outfile.write('echo "' + item + '" >> intf' + str(np.mod(i, numproc)) + '.in\n');
        # outfile.write("echo S1A20180106_ALL_F1:S1A20180118_ALL_F1 >> intf0.in\n"); break;   
    outfile.write("\n# Unwrap the interferograms.\n\n")
    outfile.write("ls intf?.in | parallel --eta 'unwrap_mod.csh {} " + config_file + "'\n\n\n");  # if you have parallel
    # outfile.write("unwrap_mod.csh intf_record.in "+config_file+"\n\n\n");  # if you don't have parallel 
    outfile.close();

    return;


def read_corr_results(corr_file):
    [stem1, stem2, mean_corr] = np.loadtxt(corr_file, unpack=True, usecols=(1, 2, 3),
                                           dtype={'names': ('stem1', 'stem2', 'mean_corr'),
                                                  'formats': ('U22', 'U22', float)})
    stem1 = [x.split('.')[0] for x in stem1];  # format: S1_20181025_ALL_F1
    stem2 = [x.split('.')[0] for x in stem2];  # format: S1_20181025_ALL_F1
    return [stem1, stem2, mean_corr];


def get_small_baseline_subsets(baseline_tuple_list, tbaseline_max, xbaseline_max):
    """ Grab all the pairs that are below the critical baselines in space and time. 
    Return format is a list of strings like 'S1A20150310_ALL_F1:S1A20150403_ALL_F1'. 
    """
    print("Getting small baseline subsets with bperp_max = %f m and t_max = %f days"
          % (xbaseline_max, tbaseline_max));
    nacq = len(baseline_tuple_list);
    intf_pairs = [];
    stems = [x[3] for x in baseline_tuple_list];
    datetimearray = [x[1] for x in baseline_tuple_list];
    xbaseline = [x[0] for x in baseline_tuple_list];
    for i in range(0, nacq - 1):
        for j in range(i + 1, nacq):
            dtdelta = datetimearray[i] - datetimearray[j];
            dtdeltadays = dtdelta.days;  # how many days exist between the two acquisitions?
            if abs(dtdeltadays) < tbaseline_max:
                if abs(xbaseline[i] - xbaseline[j]) < xbaseline_max:
                    img1_stem = stems[i];
                    img2_stem = stems[j];
                    img1_time = int(img1_stem[3:11]);
                    img2_time = int(img2_stem[3:11]);
                    if img1_time < img2_time:  # if the images are listed in chronological order
                        intf_pairs.append(stems[i] + ":" + stems[j]);
                    else:  # if the images are in reverse chronological order
                        intf_pairs.append(stems[j] + ":" + stems[i]);
                else:
                    print("WARNING: %s:%s rejected due to large perpendicular baseline of %f m."
                          % (stems[i], stems[j], abs(xbaseline[i] - xbaseline[j])));
    # The total number of pairs is (n*n-1)/2.  How many of them fit our small baseline criterion?
    total_possible_pairs = nacq * (nacq - 1) / 2;
    print("SBAS Pairs: Returning %d of %d possible interferograms to compute. "
          % (len(intf_pairs), total_possible_pairs));
    return intf_pairs;


def get_chain_subsets(baseline_tuple_list):
    """goal: order tbaselines ascending order. Then just take adjacent stems as the intf pairs.
    future idea: implement a bypass option, where we can ignore some acquisitions """
    print("CHAIN Pairs: Getting chain connections. ");
    intf_pairs = [];
    for i in range(len(baseline_tuple_list) - 1):
        intf_pairs.append(baseline_tuple_list[i][3] + ':' + baseline_tuple_list[i+1][3]);  # stems in ascending order
    print("CHAIN Pairs: Returning %d interferograms to compute from %d images. "
          % (len(intf_pairs), len(baseline_tuple_list)));
    return intf_pairs;


def filter_intf_start_end(intf_pairs, startdate, enddate):
    """ Take a list of interferograms and return the ones that fall within a given time window.
    Startdate, enddate : YYYYMMDD strings """
    if startdate == "" and enddate == "":
        return intf_pairs;
    startdate = dt.datetime.strptime(startdate, "%Y%m%d");
    enddate = dt.datetime.strptime(enddate, "%Y%m%d");
    intf_all = [];
    for item in intf_pairs:
        date1 = item[3:11];
        date2 = item[22:30];
        d1 = dt.datetime.strptime(date1, "%Y%m%d");
        d2 = dt.datetime.strptime(date2, "%Y%m%d");
        if startdate <= d1 <= enddate and startdate <= d2 <= enddate:
            intf_all.append(item);
    return intf_all;


def make_network_plot(intf_pairs, baseline_tuple_list, plotname):
    """ intf_pairs is a list of strings, has two possible formats.
    baseline_tuple_list has format given at top. """
    print("Printing network plot with %d intfs" % (len(intf_pairs)));
    if len(intf_pairs) == 0:
        print("Error! Cannot make network plot because there are no interferograms. ");
        sys.exit(1);
    stems = [x[3] for x in baseline_tuple_list];
    datetimearray = [x[1] for x in baseline_tuple_list];
    xbaseline = [x[0] for x in baseline_tuple_list];
    xstart, xend, tstart, tend = [], [], [], [];

    # If there's intf_pairs format like "S1A20160817_ALL_F2:S1A20160829_ALL_F2"
    if "S1" in intf_pairs[0]:
        for item in intf_pairs:
            scene1 = item[0:18];  # has some format like S1A20160817_ALL_F2
            scene2 = item[19:];  # has some format like S1A20160817_ALL_F2
            for x in range(len(stems)):
                if stems[x] == scene1:
                    xstart.append(xbaseline[x]);
                    tstart.append(datetimearray[x]);
                if stems[x] == scene2:
                    xend.append(xbaseline[x]);
                    tend.append(datetimearray[x]);

    # If there's intf_pairs format like "2017089:2018101"....
    if len(intf_pairs[0]) == 15:
        for i in range(len(intf_pairs)):
            scene1 = intf_pairs[i][0:7];
            scene2 = intf_pairs[i][8:15];
            im1_dt = dt.datetime.strptime(scene1, '%Y%j');
            im2_dt = dt.datetime.strptime(scene2, '%Y%j');
            for x in range(len(stems)):
                if datetimearray[x] == im1_dt:
                    xstart.append(xbaseline[x]);
                    tstart.append(datetimearray[x]);
                if datetimearray[x] == im2_dt:
                    xend.append(xbaseline[x]);
                    tend.append(datetimearray[x]);

    plt.figure();
    plt.plot_date(tstart, xstart, '.b');
    plt.plot_date(tend, xend, '.b');
    for i in range(len(tstart)):
        plt.plot_date([tstart[i], tend[i]], [xstart[i], xend[i]], 'b', linewidth=0.5);
    yrs_formatter = mdates.DateFormatter('%m-%y');
    plt.xlabel("Date");
    plt.gca().xaxis.set_major_formatter(yrs_formatter);
    plt.ylabel("Baseline (m)");
    plt.title("Network Geometry with %d Interferograms" % (len(intf_pairs)));
    plt.savefig(plotname);
    plt.close();
    print("Finished printing network plot");
    return;

#
# Reporting and defensive programming
#

def check_raw_orig_sanity(swath):
    number_of_tiffs = int(check_output('ls F' + swath + '/raw_orig/*.tiff | wc -l', shell=True));
    number_of_safes = int(check_output('ls F' + swath + '/raw_orig/*.safe | wc -l', shell=True));
    number_of_EOFs = int(check_output('ls F' + swath + '/raw_orig/*.EOF | wc -l', shell=True));
    print('number of tiffs is %d ' % number_of_tiffs);
    print('number of safes is %d ' % number_of_safes);
    print('number of EOFs is %d ' % number_of_EOFs);
    assert(number_of_tiffs == number_of_safes), DirectoryError("Error: raw_orig has non-matching tiff/safe files.")
    assert(number_of_tiffs == number_of_EOFs), DirectoryError("Error: raw_orig has non-matching tiff/EOF files.")
    return;


def check_intf_all_sanity(intended_swaths, common_intfs, filename='phase.grd'):
    """
    # Figure out whether all intended interferograms were made.
    # filename is like 'phase.grd'
    # common_intfs is in convenient format like '2019082_2019106'
    # Check that all common interferograms have a given file. """
    for swath in intended_swaths:
        for igram in common_intfs:
            if not os.path.isfile("F"+swath+'/intf_all/'+igram+'/'+filename):
                raise DirectoryError('error: file phase.grd not found in F%s/%s' % (swath, igram));
    print("Status: %d of %s have been found in all intended swaths. " % (len(common_intfs), filename));
    return;

#
# Exceptions and Exception handling
#

# A special exception for when a directory is poorly situated, and is going to fail. 
class DirectoryError(Exception):
    def __init__(self, value):
        self.value = value;

    def __str__(self):
        return repr(self.value);
