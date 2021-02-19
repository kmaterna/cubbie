# Sentinel Utilities

import subprocess
import os
import sys
import glob
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import collections
matplotlib.use('Agg')


def get_all_xml_names(directory, polarization, swath):
    # Returns a matching list of filenames and datestrs (yyyymmdd)
    pathname1 = directory + "/*-"+polarization+"-*-00" + swath + ".xml";
    pathname2 = directory + "/*-"+polarization+"-*-00" + str(int(swath) + 3) + ".xml";
    list_of_images_temp = glob.glob(pathname1) + glob.glob(pathname2);
    list_of_images = []
    list_of_datestrs = [];
    for item in list_of_images_temp:
        list_of_images.append(item[:])
        list_of_datestrs.append(get_date_from_xml(item[:]));
    return list_of_images, list_of_datestrs;

def get_all_safes_in_dir(directory):
    # A directory that contains a bunch of safe files. Give us the list of files and the datetimes they go with.
    # Matching lengths
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
    # Return all safes that happened on a particular datetime
    dirlist, datelist = get_all_safes_in_dir(dirname);
    retval = [];
    for item, date in zip(dirlist, datelist):
        if date == one_date:
            retval.append(item);
    retval = sorted(retval);
    return retval;

def get_SAFE_list_for_raw_orig(config_params):
    # Return the files we're going to put into raw_orig
    # Location depends on whether we've made frames or not.
    if config_params.frame1 != '':
        file_list, dt_list = get_all_safes_in_dir(config_params.FRAMES_dir);
        # if we're assembling frames, we use the FRAMES directory.
    else:
        file_list, dt_list = get_all_safes_in_dir(config_params.DATA_dir);
        # if we're not assembling frames, we use the DATA directory.
    return file_list, dt_list;

def get_all_tiff_names(directory, polarization, swath):
    pathname1 = directory + "/*-"+polarization+"-*-00" + swath + ".tiff";
    pathname2 = directory + "/*-"+polarization+"-*-00" + str(int(swath) + 3) + ".tiff";
    list_of_images_temp = glob.glob(pathname1) + glob.glob(pathname2);
    list_of_images = []
    for item in list_of_images_temp:
        list_of_images.append(item[:])
    return list_of_images;


def get_previous_and_following_day(datestring):
    """ This is a function that takes a date like 20160827 and generates 
    [20160826, 20160828]: the day before and the day after the date in question. """
    year = int(datestring[0:4]);
    month = int(datestring[4:6]);
    day = int(datestring[6:8]);
    mydate = dt.date(year, month, day);
    tomorrow = mydate + dt.timedelta(days=1);
    yesterday = mydate - dt.timedelta(days=1);
    previous_day = pad_string_zeros(yesterday.year) + pad_string_zeros(yesterday.month) + pad_string_zeros(
        yesterday.day);
    following_day = pad_string_zeros(tomorrow.year) + pad_string_zeros(tomorrow.month) + pad_string_zeros(tomorrow.day);
    return [previous_day, following_day];


def get_date_from_xml(xml_name):
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


def pad_string_zeros(num):
    if num < 10:
        numstring = "0" + str(num);
    else:
        numstring = str(num);
    return numstring;


def ymd2yj(ymd):
    # Turn something like "20150524" into "2015100", useful for file naming conventions
    tdate = dt.datetime.strptime(ymd, "%Y%m%d");
    tdate = tdate - dt.timedelta(days=1);
    return dt.datetime.strftime(tdate, "%Y%j");


def yj2ymd(yj):
    yj = int(yj) + 1;
    tdate = dt.datetime.strptime(str(yj), "%Y%j");
    return dt.datetime.strftime(tdate, "%Y%m%d");


def get_eof_from_date_sat(mydate, sat, eof_dir):
    """ This returns something like S1A_OPER_AUX_POEORB_OPOD_20160930T122957_V20160909T225943_20160911T005943.EOF.
        eof_dir can be relative or absolute
        It takes something like 20171204, s1a, eof_dir
        It can also take something like dt.datetime, s1a, eof_dir
    """
    if type(mydate) == dt.datetime:   # if you pass a datetime object
        mydate = dt.datetime.strftime(mydate, "%Y%m%d");
    [previous_day, following_day] = get_previous_and_following_day(mydate);
    eof_name = glob.glob(eof_dir + "/" + sat.upper() + "*" + previous_day + "*" + following_day + "*.EOF");
    if not eof_name:
        print("ERROR: did not find any EOF files matching the pattern " + eof_dir + "/" + sat.upper() +
              "*" + previous_day + "*" + following_day + "*.EOF");
        print("Exiting...")
        sys.exit(1);
    else:
        eof_name = eof_name[0];
    return eof_name;


def glob_intf_computed():
    full_names = glob.glob("intf_all/*");
    intf_computed = [];
    for item in full_names:
        intf_computed.append(item[9:]);
    return intf_computed;


def compare_frames_with_safes(config_params):
    # checking for output consistency after frames are created.
    print("Checking the frames and acquisitions and see if all your data has been included. ");
    dirlist, datelist = get_all_safes_in_dir(config_params.DATA_dir);
    frames_dirlist, frames_datelist = get_all_safes_in_dir(config_params.FRAMES_dir);
    print("BEGIN: %d SAFES from %d unique dates." % (len(dirlist), len(set(datelist))));
    print("END: %d SAFES from %d unique dates." % (len(frames_dirlist), len(set(frames_datelist))));
    return;


def make_data_in(polarization, swath, master_date=""):
    """
    data.in is a reference table that links the xml file with the correct orbit file.
    """
    list_of_images, list_of_datestrs = get_all_xml_names("F" + str(swath) + "/raw_orig", polarization, swath);
    outfile = open("F" + str(swath) + "/data.in", 'w');
    if master_date == "":
        print("No master date selected by the user. Printing in random order.")
        for item, mydate in zip(list_of_images, list_of_datestrs):
            item = item.split("/")[-1];  # getting rid of the directory
            sat = get_sat_from_xml(item);
            eof_name = get_eof_from_date_sat(mydate, sat, "F" + str(swath) + "/raw_orig");
            outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
    else:
        # write the master date first. 
        for item, mydate in zip(list_of_images, list_of_datestrs):
            if mydate == master_date:
                print("Found master date %s in config file. Putting it on first line." % master_date);
                item = item.split("/")[-1];  # getting rid of the directory
                sat = get_sat_from_xml(item);
                eof_name = get_eof_from_date_sat(mydate, sat, "F" + str(swath) + "/raw_orig");
                outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
                break;
        # then write the other dates. 
        for item, mydate in zip(list_of_images, list_of_datestrs):
            if mydate != master_date:
                item = item.split("/")[-1];  # getting rid of the directory
                sat = get_sat_from_xml(item);
                eof_name = get_eof_from_date_sat(mydate, sat, "F" + str(swath) + "/raw_orig");
                outfile.write(item[:-4] + ":" + eof_name.split("/")[-1] + "\n");
    outfile.close();
    print("data.in successfully printed.")
    return;


def read_baseline_table(baselinefilename):
    # Returns a chronologically sorted list of tuples of (baseline, datetime, datestr, stem) values
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
    baselines = baselineFile[:, 4].astype(float);

    dtarray, datestrs = [], [];
    for i in range(len(times)):
        dtarray.append(dt.datetime.strptime(str(int(times[i] + 1)), '%Y%j'));
        datestrs.append(str(times[i] + 1)[0:7]);   # string with format "2014361"

    # Re-order times and baselines in chronological order, and re-package
    baselines = [x for _, x in sorted(zip(dtarray, baselines))];
    datestrs = [x for _, x in sorted(zip(dtarray, datestrs))];
    stems = [x for _, x in sorted(zip(stems, datestrs))];
    dtarray = sorted(dtarray);
    baseline_tuple_list = [];
    for i in range(len(baselines)):
        baseline_tuple_list.append((baselines[i], dtarray[i], datestrs[i], stems[i]));

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
def choose_master_image(default_master, swath):
    print("Selecting master image; default_master is %s" % default_master);

    # load baseline table
    baselineFile = np.genfromtxt('F' + str(swath) + '/raw/baseline_table.dat', dtype=str)
    time = baselineFile[:, 1].astype(float)
    baseline = baselineFile[:, 4].astype(float)
    shortform_names = baselineFile[:, 0].astype(str);

    # GMTSAR (currently) guarantees that this file has the same order of lines as baseline_table.dat.
    dataDotIn = np.genfromtxt('F' + str(swath) + '/raw/data.in', dtype='str').tolist()
    print(dataDotIn);

    if default_master == "":
        # calculate shortest distance from median to scenes
        consider_time = True
        if consider_time:
            time_baseline_scale = 1  # arbitrary scaling factor, units of (meters/day)
            sceneDistance = np.sqrt(
                ((time - np.median(time)) / time_baseline_scale) ** 2 + (baseline - np.median(baseline)) ** 2)
        else:
            sceneDistance = np.sqrt((baseline - np.median(baseline)) ** 2)

        minID = np.argmin(sceneDistance)
        masterID = dataDotIn[minID]

    else:  # if the user has selected a master on purpose
        for x in range(len(dataDotIn)):
            if default_master in dataDotIn[x]:
                masterID = dataDotIn[x];
                minID = x;

    # put masterId in the first line of data.in
    dataDotIn.pop(dataDotIn.index(masterID))
    dataDotIn.insert(0, masterID)
    master_shortform = shortform_names[minID];  # because GMTSAR initially puts baseline_table / data.in in same order.

    os.rename('F' + str(swath) + '/raw/data.in', 'F' + str(swath) + '/raw/data.in.old')
    np.savetxt('F' + str(swath) + '/raw/data.in', dataDotIn, fmt='%s')
    np.savetxt('F' + str(swath) + '/data.in', dataDotIn, fmt='%s')
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


def set_up_merge_unwrap(config_params):
    subprocess.call(["mkdir", "-p", "merged"], shell=False);
    subprocess.call(["cp", "topo/dem.grd", "merged"],
                    shell=False);  # I needed to actually copy it, not just soft link.
    subprocess.call(["cp", "batch.config", "merged"], shell=False);
    return;


def write_merge_batch_input(intf_all, master_image):
    # Useful for the merge swaths step. 
    # The master of the first line should be the super master. 
    # In the future, would be good to defensively check whether each directory contains phasefilt, mask, and corr. 
    ofile = open("merged/inputfile.txt", "w");

    # Find the super master and stick it to the front.
    for item in intf_all:
        first_part = item.split(":")[0];
        if master_image.split()[0][0:-3] in first_part:
            print("Found master image in intf %s" % item);
            print("Moving this image to the front");
            intf_all.insert(0, intf_all.pop(intf_all.index(item)));
            break;

    for item in intf_all:
        first_part = item.split(":")[0];
        second_part = item.split(":")[1];
        jd1 = ymd2yj(item[3:11]);  # a string
        jd2 = ymd2yj(item[22:30]);  # a string
        ofile.write("../F1/intf_all/%s_%s/:" % (jd1, jd2));
        ofile.write(first_part + ".PRM:" + second_part + ".PRM,");
        ofile.write("../F2/intf_all/%s_%s/:" % (jd1, jd2));
        ofile.write(first_part.replace("_F1", "_F2") + ".PRM:" + second_part.replace("_F1", "_F2") + ".PRM,");
        ofile.write("../F3/intf_all/%s_%s/:" % (jd1, jd2));
        ofile.write(first_part.replace("_F1", "_F3") + ".PRM:" + second_part.replace("_F1", "_F3") + ".PRM\n");
    ofile.close();
    return;


def write_merge_unwrap(outfile):
    print("Writing unwrap instructions in file %s " % outfile)
    ofile = open(outfile, 'w');
    ofile.write("cd merged\n");
    ofile.write("pwd\n");
    ofile.write("ls\n");
    ofile.write("merge_batch.csh inputfile.txt batch.config\n")
    ofile.write("cd ../\n");
    ofile.close();
    return;


def write_ordered_unwrapping(numproc, swath, sh_file, config_file):
    [stem1, stem2, mean_corr] = read_corr_results("corr_results.txt");

    stem1_ordered = [x for y, x in sorted(zip(mean_corr, stem1), reverse=True)];
    stem2_ordered = [x for y, x in sorted(zip(mean_corr, stem2), reverse=True)];
    mean_corr_ordered = sorted(mean_corr, reverse=True);

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


def write_select_unwrapping(numproc, swath, sh_file, config_file):
    intfs = ["S1_20190411_ALL_F" + swath + ":S1_20190423_ALL_F" + swath,
             "S1_20190330_ALL_F" + swath + ":S1_20190411_ALL_F" + swath,
             "S1_20190330_ALL_F" + swath + ":S1_20190423_ALL_F" + swath];

    outfile = open(sh_file, 'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch unwrap Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("cd F" + swath + "\n");
    outfile.write("rm intf?.in\n");
    for i, item in enumerate(intfs):
        outfile.write('echo "' + intfs[i] + '" >> intf' + str(np.mod(i, numproc)) + '.in\n');
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


def remove_nans_array(myarray):
    numarray = [];
    for i in range(len(myarray)):
        if ~np.isnan(myarray[i]):
            numarray.append(myarray[i][0]);
    return numarray;


def get_small_baseline_subsets(stems, tbaseline, xbaseline, tbaseline_max, xbaseline_max):
    """ Grab all the pairs that are below the critical baselines in space and time. 
    Return format is a list of strings like 'S1A20150310_ALL_F1:S1A20150403_ALL_F1'. 
    You can adjust this if you have specific processing needs. 
    """
    print("SBAS: Getting small baseline subsets with bperp_max = %f m and t_max = %f days"
          % (xbaseline_max, tbaseline_max));
    nacq = len(stems);
    intf_pairs = [];
    datetimearray = [];
    for k in tbaseline:
        datetimearray.append(dt.datetime.strptime(str(int(k) + 1), "%Y%j"));  # convert to datetime arrays.
    # print(datetimearray);
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


def get_chain_subsets(stems, tbaseline):
    # goal: order tbaselines ascending order. Then just take adjacent stems as the intf pairs.
    # future idea: implement a bypass option, where we can ignore some acquisitions
    print("CHAIN Pairs: Getting chain connections. ");
    intf_pairs = [];
    sorted_stems = [x for _, x in sorted(zip(tbaseline, stems))];  # sort by increasing t value
    for i in range(len(sorted_stems) - 1):
        intf_pairs.append(sorted_stems[i] + ':' + sorted_stems[i + 1]);
    print("CHAIN Pairs: Returning %d interferograms to compute from %d images. " % (len(intf_pairs), len(stems)));
    return intf_pairs;


def reduce_by_start_end_time(intf_pairs, startdate, enddate):
    # Take a list of interferograms and return the ones that fall within a given time window. 
    intf_all = [];
    for item in intf_pairs:
        date1 = item[3:11];
        date2 = item[22:30];
        d1 = dt.datetime.strptime(date1, "%Y%m%d");
        d2 = dt.datetime.strptime(date2, "%Y%m%d");
        if startdate <= d1 <= enddate and startdate <= d2 <= enddate:
            intf_all.append(item);
    return intf_all;


def make_network_plot(intf_pairs, stems, tbaseline, xbaseline, plotname):
    # intf_pairs is a list of strings, has two possible formats
    # stems is a list of strings, has format "S1_20150514_ALL_F2"
    # tbaseline is list of floats, as format 2015133.5774740
    # xbaseline is list of floats, in meters
    print("Printing network plot with %d intfs" % (len(intf_pairs)));
    if len(intf_pairs) == 0:
        print("Error! Cannot make network plot because there are no interferograms. ");
        sys.exit(1);

    xstart, xend, tstart, tend = [], [], [], [];

    # For intf_pairs, if there's a format like "S1A20160817_ALL_F2:S1A20160829_ALL_F2"
    if "S1" in intf_pairs[0]:
        for item in intf_pairs:
            scene1 = item[0:18];  # has some format like S1A20160817_ALL_F2
            scene2 = item[19:];  # has some format like S1A20160817_ALL_F2
            for x in range(len(stems)):
                if stems[x] == scene1:
                    xstart.append(xbaseline[x]);
                    tstart.append(dt.datetime.strptime(str(int(tbaseline[x]) + 1), '%Y%j'));
                if stems[x] == scene2:
                    xend.append(xbaseline[x]);
                    tend.append(dt.datetime.strptime(str(int(tbaseline[x]) + 1), '%Y%j'));

    # If there's a format like "2017089:2018101"....
    if len(intf_pairs[0]) == 15:
        dtarray = [];
        im1_dt = [];
        im2_dt = [];
        # Making a list of acquisition dates
        for i in range(len(tbaseline)):
            dtarray.append(dt.datetime.strptime(str(int(tbaseline[i]+1)), '%Y%j'));

        # Make the list of datetimes for the images. 
        for i in range(len(intf_pairs)):
            scene1 = intf_pairs[i][0:7];
            scene2 = intf_pairs[i][8:15];
            im1_dt.append(dt.datetime.strptime(scene1, '%Y%j'));
            im2_dt.append(dt.datetime.strptime(scene2, '%Y%j'));

        # Find the appropriate image pairs and baseline pairs
        for i in range(len(intf_pairs)):
            for x in range(len(dtarray)):
                if dtarray[x] == im1_dt[i]:
                    xstart.append(xbaseline[x]);
                    tstart.append(dtarray[x]);
                if dtarray[x] == im2_dt[i]:
                    xend.append(xbaseline[x]);
                    tend.append(dtarray[x]);

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

def compare_intended_list_with_directory(intended_array, actual_array, errormsg):
    # This takes two lists of dates formatted as strings, such as ['2015321_2015345']
    # It prints out any members that exist in the intended_array but not the actual_array
    for item in intended_array:
        if item in actual_array:
            continue;
        else:
            print("ERROR! %s expected, but not found in actual array." % item);
            print(errormsg);
    return;


def check_intf_all_sanity():
    # Figure out whether all intended interferograms were made and unwrapped. 
    # This makes sense for a single swath, not a merged situation. 

    print("Checking the progress of intf and unwrap steps. ")

    # The hardest part: Fix the differences between datetime formats in intf_record.in
    intended_intfs = np.genfromtxt('intf_record.in', dtype='str');
    intended_intfs = [i[3:11] + '_' + i[22:30] for i in
                      intended_intfs];  # these come in formatted S1A20161206_ALL_F1:S1A20161230_ALL_F1
    date1 = [dt.datetime.strptime(i[0:8], "%Y%m%d") - dt.timedelta(days=1) for i in intended_intfs];
    date2 = [dt.datetime.strptime(i[9:17], "%Y%m%d") - dt.timedelta(days=1) for i in intended_intfs];
    date1 = [dt.datetime.strftime(i, "%Y%j") for i in date1];
    date2 = [dt.datetime.strftime(i, "%Y%j") for i in date2];
    intended_intfs = [];
    for i in range(len(date1)):
        intended_intfs.append(date1[i] + '_' + date2[i]);
    num_intended = len(set(intended_intfs));
    print("  intended interferograms: %d from intf_record.in" % len(intended_intfs));
    print("  unique intended interferograms: %d " % num_intended);

    # Check for duplicated items in intf_record.in (may exist);
    duplicates = [item for item, count in collections.Counter(intended_intfs).items() if count > 1];
    print("  duplicated elements in intf_record.in: ");
    print(duplicates);

    # Collect the actual intf_all directories
    actual_intfs = subprocess.check_output('ls -d intf_all/201*_201* ', shell=True);
    actual_intfs = actual_intfs.split('\n');
    actual_intfs = [value for value in actual_intfs if value != ''];
    actual_intfs = [i.split('/')[-1] for i in actual_intfs];
    print("  actual interferograms: %d from intf_all directory " % len(actual_intfs));

    # Collect the unwrap.grd files
    actual_unwraps = subprocess.check_output('ls intf_all/unwrap.grd/*_unwrap.grd', shell=True);
    actual_unwraps = actual_unwraps.split('\n');
    actual_unwraps = [value for value in actual_unwraps if value != ''];
    actual_unwraps = [i.split('/')[-1] for i in actual_unwraps];
    actual_unwraps = [i.split('_unwrap.grd')[0] for i in actual_unwraps];
    print('  unwrapped interferograms: %d from intf_all/unwrap.grd directory' % len(actual_unwraps))

    if num_intended == len(actual_intfs):
        print("Congratulations! All of your interferograms have been made. ");
    else:
        compare_intended_list_with_directory(intended_intfs, actual_intfs, 'is not made.');
    if num_intended == len(actual_unwraps):
        print("Congratulations! All of your interferograms have been unwrapped. ");
    else:
        compare_intended_list_with_directory(intended_intfs, actual_unwraps, 'is not unwrapped.');

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
