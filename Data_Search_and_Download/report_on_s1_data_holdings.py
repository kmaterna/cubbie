#!/usr/bin/env python3

"""
The purpose of this script is to compare a directory of SAFE files on one track with the region you might want.
Did you get all the data you wanted? How many bursts?
Did every day cover the entire domain of interest?
Are there extraneous files in your directory?
Output Table:
# Options: Polygon, DataDir, DesiredDates, TodaysDate
# DesiredDates PresentDates NumFiles NumBursts CoverWhole? Continuous?
# TABLE...
# TABLE...
"""
import sys, glob
import numpy as np
import datetime as dt
import xml.etree.ElementTree as ET

def parse_args(argv):
    if len(argv) == 1 or argv[1] == "--help" or argv[1] == '--h':
        print("\nWelcome to a reporting script on your S1 data holdings.");
        print("First argument (required): directory where data lives (ex: DATA/)");
        print("Second argument (required): polygon file");
        print("Third argument (optional): desired dates file\n");
        print("Usage: report_on_s1_data_holdings.py datadir polygon.txt list_of_dates.txt\n");
        sys.exit(0);
    else:
        if len(argv) != 3 and len(argv) != 4:
            print("Error! Wrong number of input arguments. Must have 2 or 3 inputs. See help. ");
            sys.exit(0);
        datadir = argv[1];
        polygon_file = argv[2];
        if len(argv) == 4:
            list_of_dates = argv[3];
        else:
            list_of_dates = ();
    return datadir, polygon_file, list_of_dates;

def make_summary(datadir, polygon_file, list_of_dates=()):
    list_of_desired_dates = read_list_of_desired_dates(list_of_dates);
    polygon = read_polygon(polygon_file);
    safe_files = glob_safe_files(datadir);
    unique_dates = safe_list_to_unique_dates(safe_files);
    information_tuples = get_information_tuples(unique_dates, safe_files, polygon);
    write_table(list_of_desired_dates, information_tuples, "data_holdings_table.txt");
    return;

def glob_safe_files(datadir):
    filelist = glob.glob(datadir+"/*.SAFE");
    return filelist;

def safe_list_to_unique_dates(safe_files):
    # Get unique dates from list of files on the file system
    dates = [];
    for item in safe_files:
        safe_name = item.split('/')[-1];
        full_date = (safe_name.split('_')[5]);
        dates.append(dt.datetime.strptime(full_date[0:8], "%Y%m%d"));
    unique_dates = list(sorted(set(dates)));
    return unique_dates;

def read_list_of_desired_dates(list_of_dates):
    # Will read list of dates desired by user, in text file with standard format.
    if list_of_dates == ():
        return [];
    date_list = [];
    ifile = open(list_of_dates, 'r');
    for line in ifile:
        newdate = dt.datetime.strptime(line.split()[0], "%Y%m%d");
        date_list.append(newdate);
    ifile.close();
    return date_list;

def read_polygon(polygon_file):
    # Read polygon file
    polygon = [];
    ifile = open(polygon_file, 'r');
    for line in ifile:
        polygon.append([float(line.split()[0]), float(line.split()[1])]);
    ifile.close();
    return polygon;

def get_files_by_date(safe_files, date):
    # filter by date
    selected_files = [];
    for file in safe_files:
        if dt.datetime.strftime(date, "%Y%m%d") in file:
            selected_files.append(file);
    return selected_files;

def determine_continuous_safes(safe_files):
    # determine if a few safe directories involve continuous observation (no missing bursts)
    total_burst_times = [];
    for item in safe_files:
        f1_annotation_file = glob.glob(item + '/annotation/*slc*-001.xml');
        tree = ET.parse(f1_annotation_file[0]);
        root = tree.getroot();

        for child in root.iter('burst'):  # an Element 'burst'
            datestr = child.find('sensingTime').text;   # an Element 'sensingTime'
            total_burst_times.append(dt.datetime.strptime(datestr[0:19], "%Y-%m-%dT%H:%M:%S"));

    continuous = 'continuous';
    total_burst_times = sorted(set(total_burst_times));
    for i in range(len(total_burst_times)-1):
        time_diff = total_burst_times[i+1] - total_burst_times[i];
        if time_diff.seconds > 5:  # If there are more than 5 seconds between bursts, we have a discontinuity.
            print(total_burst_times[i+1], total_burst_times[i], ": DISCONTINUOUS OBSERVATIONS");
            continuous = 'gaps';

    return len(total_burst_times), continuous;

def determine_total_coverage(safe_files, polygon):
    # Given a list of several SAFE files and a polygon, is the SAFE footprint larger?
    # This function is run once for each date.
    lons_all = [];
    lats_all = [];
    for i in range(len(safe_files)):
        lons, lats = read_kml(safe_files[i]+'/preview/map-overlay.kml');
        lons_all.extend(lons);
        lats_all.extend(lats);

    # Determine the min/max of the domain of the safe files
    lon_max = np.max(lons_all);
    lat_max = np.max(lats_all);
    lon_min = np.min(lons_all);
    lat_min = np.min(lats_all);

    # Determine whether the polygon is completely covered by the domain.
    polygon_lons = [x[0] for x in polygon];
    polygon_lats = [x[1] for x in polygon];
    print("data: ", lon_min, lon_max, lat_min, lat_max);
    print("polygon: ", np.min(polygon_lons), np.max(polygon_lons), np.min(polygon_lats), np.max(polygon_lats));

    if lon_min < np.min(polygon_lons) and lat_min < np.min(polygon_lats) and lon_max > np.max(polygon_lons) \
            and lat_max > np.max(polygon_lats):
        total_coverage = 'covered';
    else:
        total_coverage = 'missing';
        print("Lacking total coverage...")

    return total_coverage;

def read_kml(infile):
    lats, lons = [], [];
    ifile = open(infile, 'r');
    for line in ifile:
        if "coordinates" in line:
            temp = line.split('>')[-2];
            temp = temp.split('<')[0];
            temp = temp.split();
            for item in temp:
                lons.append(float(item.split(',')[0]))
                lats.append(float(item.split(',')[1]))
    return lons, lats;

def get_information_tuples(unique_dates, safe_files, polygon):
    # Construct a tuple: date, num_files, num_bursts, cover_whole, continuous
    list_of_all_tuples = [];
    for date in unique_dates:
        print("\nFor date %s: " % dt.datetime.strftime(date, "%Y-%m-%d"));
        selected_files = get_files_by_date(safe_files, date);
        cover_whole = determine_total_coverage(selected_files, polygon);  # inspect bursts for total range
        num_bursts, continuous = determine_continuous_safes(selected_files);  # inspect bursts for continuity
        information_tuple = (date, len(selected_files), num_bursts, cover_whole, continuous);
        list_of_all_tuples.append(information_tuple);
    return list_of_all_tuples;

def write_table(list_of_desired_dates, information_tuples, outfile):
    ofile = open(outfile, 'w');
    ofile.write("# Options: Polygon, DataDir, DesiredDates, TodaysDate\n");
    ofile.write("# Columns: DesiredDates | PresentDates NumFiles NumBursts CoverWhole? Continuous?\n")
    list_of_present_dates = [x[0] for x in information_tuples];
    total_dates = list(sorted(set(list_of_desired_dates + list_of_present_dates)));
    ofile.write("# Dates: desired=%d, downloaded=%d, total=%d\n" % (len(list_of_desired_dates),
                                                                    len(information_tuples), len(total_dates)));
    for one_date in total_dates:
        date_as_string = dt.datetime.strftime(one_date, "%Y%m%d");
        if one_date in list_of_desired_dates:
            ofile.write(date_as_string+" | ");
        else:
            ofile.write("         | ");
        if one_date in list_of_present_dates:
            idx = list_of_present_dates.index(one_date);
            ofile.write("%s %d %d %s %s\n" % (date_as_string, information_tuples[idx][1], information_tuples[idx][2],
                                              information_tuples[idx][3], information_tuples[idx][4]));
        else:
            ofile.write("\n");
    ofile.close();
    return;


if __name__ == "__main__":
    datadir, polygon, list_of_dates = parse_args(sys.argv);
    make_summary(datadir, polygon, list_of_dates);
