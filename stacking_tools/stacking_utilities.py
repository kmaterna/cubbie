# Stacking Utilities

import os
import sys
import glob
import re
import collections
import datetime as dt
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import haversine
from s1_batches.intf_generating import get_ra_rc_from_ll
from s1_batches.read_write_insar_utilities import isce_read_write


def get_list_of_intf_all(config_params):
    """
    This is mechanical: just takes the list of interferograms in intf_all.
    It is designed to work with both ISCE and GMTSAR files.
    The more advanced selection takes place in make_selection_of_intfs.
    :returns: list of tuples like : (dt1, dt2, intf_file, corr_file).
    """
    if config_params.SAT == "S1":
        total_intf_list = glob.glob(config_params.intf_dir + "/???????_???????/"+config_params.intf_filename)
        total_corr_list = glob.glob(config_params.intf_dir + "/???????_???????/"+config_params.corr_filename)
        intf_file_tuples = get_intf_datetuple_gmtsar(total_intf_list, total_corr_list)
    elif config_params.SAT == "UAVSAR":
        # Specific to the case of UAVSAR stacks with alt-unwrapped taking place
        total_intf_list = glob.glob("../Igrams/*/alt_unwrapped/filt*_fully_processed.uwrappedphase")
        total_corr_list = glob.glob("../Igrams/*/alt_unwrapped/filt*_fully_processed.cor")
        intf_file_tuples = get_intf_datetuple_isce(total_intf_list, total_corr_list)
    else:
        total_intf_list, total_corr_list, intf_file_tuples = [], [], []
    print("Identifying all unwrapped intfs in %s: " % config_params.intf_dir)
    print("  Found %d interferograms for stacking. " % (len(total_intf_list)))
    print("  Found %d coherence files for stacking. " % (len(total_corr_list)))
    if len(total_intf_list) != len(total_corr_list):
        print("ERROR!  Length of intfs does match length of coherence files!  Please fix this before proceeding. \n")
        sys.exit(1)
    return intf_file_tuples


def get_intf_datetuple_gmtsar(total_intf_list, total_corr_list):
    """
    Turn GMTSAR interferograms into date-date-filename tuples
    :param total_intf_list: list of strings, filenames
    :param total_corr_list: list of strings, filenames
    :returns: list of tuples (dt, dt, filename, filename)
    """
    intf_tuple_list = []
    for i in range(len(total_intf_list)):
        datesplit = re.findall(r"\d\d\d\d\d\d\d_\d\d\d\d\d\d\d", total_intf_list[i])[0]  # example: 2010040_2014064
        date1 = dt.datetime.strptime(str(int(datesplit[0:7]) + 1), "%Y%j")
        date2 = dt.datetime.strptime(str(int(datesplit[8:15]) + 1), "%Y%j")  # adding 1 because 000 = January 1
        intf_tuple_list.append((date1, date2, total_intf_list[i], total_corr_list[i]))
    return intf_tuple_list


def get_intf_datetuple_isce(total_intf_list, total_corr_list):
    """
    Turn ISCE interferograms into date-date-filename tuples
    :param total_intf_list: list of strings, filenames
    :param total_corr_list: list of strings, filenames
    :returns: list of tuples (dt, dt, filename, filename)
    """
    intf_tuple_list = []
    for i in range(len(total_intf_list)):
        datesplit = re.findall(r"\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\d\d", total_intf_list[i])[0]  # ex: 20100402_20140304
        date1 = dt.datetime.strptime(datesplit[0:8], "%Y%m%d")
        date2 = dt.datetime.strptime(datesplit[9:17], "%Y%m%d")
        intf_tuple_list.append((date1, date2, total_intf_list[i], total_corr_list[i]))
    return intf_tuple_list


def get_unique_dts_from_intf_dates(date_pairs_dt):
    """
    Construct a list of unique datetimes covered by a network of intf pairs
    :param date_pairs_dt: list of objects of form [dt1, dt2]
    :returns: list of datetime objects
    """
    total_dates = []
    for item in date_pairs_dt:
        total_dates.append(item[0])
        total_dates.append(item[1])
    xdates = sorted(set(total_dates))
    return xdates


def get_TS_dates(date_julstrings):
    """"
    Get the x axis associated with a certain set of interferograms
    :param date_julstrings: list of N date_julstrings in format [YYYYJJJ_YYYJJJ,...]
    :returns: lists, length M, associated with TS. Strings in format YYYYJJJ, dts, and number of days since first image.
    """
    dates_total = []
    for i in range(len(date_julstrings)):
        dates_total.append(date_julstrings[i][0:7])
        dates_total.append(date_julstrings[i][8:15])
    datestrs = sorted(set(dates_total))
    x_axis_datetimes = [dt.datetime.strptime(x, "%Y%j") for x in datestrs]
    x_axis_days = [(x - x_axis_datetimes[0]).days for x in
                   x_axis_datetimes]  # number of days since first acquisition.
    return datestrs, x_axis_datetimes, x_axis_days


def get_list_of_ts_grids(config_params):
    """
    Glob function. Used instead of regular reader for making velocities out of pre-existing time series grids.
    """
    ts_slice_files = glob.glob(config_params.intf_dir + "/????????.grd")
    if len(ts_slice_files) == 0:
        print("Error! Not starting with any time series slices. Exiting.")
        sys.exit(0)
    return ts_slice_files


#  ----------- Reference Pixel Math ------------  #
def get_ref_index(ref_loc, ref_idx, geocoded_flag, intf_files, signalspread_filename):
    """
    Get the index of the reference pixel (generally using merged-subswath files)
    If you don't have the reference pixel in the config file,
    the program will stop execution so you can write it there.
    :param ref_loc: string, lon/lat
    :param ref_idx: string, int/int
    :param geocoded_flag: bool
    :param intf_files: list of tuples of (dt, dt, filename, filename)
    :param signalspread_filename: string, filename of grid file
    """
    print("Identifying reference pixel:")
    if ref_idx == "" and ref_loc == "":
        get_100p_pixels_manually_choose(signalspread_filename)
        sys.exit(0)
    elif ref_idx != "":  # First preference: Get the reference pixel index from the config file
        rowref = int(ref_idx.split('/')[0])
        colref = int(ref_idx.split('/')[1])
        print("  Found rowref, colref = %d, %d from config file.\n  Moving on..." % (rowref, colref))
    else:  # here we have a reference pixel from lat/lon
        lon = float(ref_loc.split('/')[0])
        lat = float(ref_loc.split('/')[1])
        if geocoded_flag:   # Second preference: Extract the lat/lon from already-geocoded intf_files
            rowref, colref = get_reference_pixel_from_geocoded_grd(lon, lat, intf_files[0][2])
            # Would use uavsar_from_lonlat_get_rowcol if doing UAVSAR here.
        else:  # Last preference: Extract the lat/lon from radar coordinates using trans.dat of merged-subswath files
            trans_dat = "topo/trans.dat"   # could be topo/trans.dat or merged/trans.data
            rowref, colref = get_referece_pixel_from_radarcoord_grd(lon, lat, trans_dat, intf_files[0][2])
        print("\nSTOP! Please write the reference row/col %d/%d into your config file. \n" % (rowref, colref))
        sys.exit(1)
    return rowref, colref


def get_100p_pixels_manually_choose(signalspread_filename):
    """
    This iterative function helps you manually choose a reference pixel based on signalspread.nc.
    You might have to run through this function a number of times
    To select your boxes and your eventual reference pixel (row, col).
    I pick one in a stable area, outside of the deformation, ideally in a desert.
    """
    print("Finding the pixels that are 100 percent coherent from signalspread")
    [x, y, ss_data] = netcdf_read_write.read_any_grd(signalspread_filename)

    count = 0
    candidate_options, great_options = [], []
    for i in range(len(y)):
        for j in range(len(x)):
            if ss_data[i][j] == 100:
                count = count + 1
                candidate_options.append((i, j))
    total_pixels = np.multiply(np.shape(ss_data)[0], np.shape(ss_data[1]))
    print("%d of %d (%f percent) are totally coherent. " % (count, total_pixels, 100 * (count / total_pixels)))
    print("%d pixels are good options for the reference pixel. " % (len(candidate_options)))

    # Manual select: great options
    xrange_great = (50, 250)     # manual select here
    yrange_great = (900, 1000)  # manual select here
    for item in candidate_options:
        if yrange_great[0] < item[0] < yrange_great[1]:
            if xrange_great[0] < item[1] < xrange_great[1]:
                great_options.append(item)
    print("%d pixels are great options for the reference pixel. " % (len(great_options)))

    # Manual select: Get a single reference pixel from the great options
    refpix = great_options[40]   # manual select here

    # Make a plot that shows where those pixels are
    plt.figure()
    plt.imshow(ss_data, aspect=1, cmap='rainbow')
    plt.gca().invert_yaxis()
    plt.plot([item[1] for item in candidate_options], [item[0] for item in candidate_options], '.', color='k')
    plt.plot([item[1] for item in great_options], [item[0] for item in great_options], '.', color='g')
    plt.plot(refpix[1], refpix[0], '.', color='r')
    plt.savefig('best_pixels.png')
    plt.close()
    print("Based on 100p pixels, selecting reference pixel at row/col %d/%d " % (refpix[0], refpix[1]))
    print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.")
    return refpix[0], refpix[1]


def uavsar_from_lonlat_get_rowcol(config_params):
    """
    An alternative way of choosing a reference pixel with UAVSAR files
    Given a ref loc, get the geocoded grids and find the nearest point.
    Return its row and column.
    Step 1: Geocode properly based on a sample interferogram grid (done)
    Step 2: extract nearest pixel (code in Brawley repo)
    """

    # isce_geocode_tools.geocode_UAVSAR_stack(config_params, 'geocoded')

    reflon = float(config_params.ref_loc.split(',')[0])
    reflat = float(config_params.ref_loc.split(',')[1])
    # Next we get the nearest pixel from the rasters
    _, _, raster_lon = isce_read_write.read_scalar_data(config_params.ts_output_dir + "/cut_lon.gdal")
    _, _, raster_lat = isce_read_write.read_scalar_data(config_params.ts_output_dir + "/cut_lat.gdal")
    i_found, j_found = get_nearest_pixel_in_raster(raster_lon, raster_lat, reflon, reflat)
    print("From lon/lat, found Row and Column at %d, %d " % (i_found, j_found))
    print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.")
    return i_found, j_found


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat):
    """
    A very general function
    Take a 2D raster of lons and lats and find the grid location closest to the target location
    """
    dist = np.zeros(np.shape(raster_lon))
    lon_shape = np.shape(raster_lon)
    for i in range(lon_shape[0]):
        for j in range(lon_shape[1]):
            mypt = [raster_lat[i][j], raster_lon[i][j]]
            dist[i][j] = haversine.distance((target_lat, target_lon), mypt)
    minimum_distance = np.nanmin(dist)
    if minimum_distance < 0.25:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist))
        i_found = idx[0][0]
        j_found = idx[1][0]
        print(raster_lon[i_found][j_found], raster_lat[i_found][j_found])
    else:
        i_found, j_found = -1, -1  # error codes
    return i_found, j_found


def get_referece_pixel_from_radarcoord_grd(lon, lat, trans_dat, example_grd):
    print("  Finding coordinate %.4f, %.4f in radarcoord grd file %s" % (lon, lat, trans_dat))
    [ra, az] = get_ra_rc_from_ll.get_ra_from_ll(trans_dat, example_grd, lon, lat)
    if np.isnan(ra) or np.isnan(az):
        print("WARNING: Cannot Find %f %f in file." % (lon, lat))
        rowref = np.nan
        colref = np.nan
    else:
        [rowref, colref] = get_ra_rc_from_ll.get_nearest_row_col(example_grd, ra, az)
        print("  Found Coordinate at row, col %d, %d " % (rowref, colref))
    return rowref, colref


def get_reference_pixel_from_geocoded_grd(ref_lon, ref_lat, ifile):
    """ Find the nearest pixel to a reference point in a geocoded grid"""
    print("  Finding coordinate %.4f, %.4f in geocoded interferograms %s" % (ref_lon, ref_lat, ifile))
    [xdata, ydata, _] = netcdf_read_write.read_netcdf4(ifile)
    if xdata[0] > 180:   # If numbers are in the range above 180, turn them into -180 to 180
        xdata = [i-360 for i in xdata]
    row_idx = np.argmin(np.abs(np.array(ydata) - ref_lat))
    col_idx = np.argmin(np.abs(np.array(xdata) - ref_lon))
    if row_idx == 0 or row_idx == len(ydata) or col_idx == 0 or col_idx == len(xdata):
        print("WARNING: Coordinate %f %f may be near edge of domain." % (ref_lon, ref_lat))
        row_idx = np.nan
        col_idx = np.nan
    print("  Found Coordinates at row/col: %d/%d " % (row_idx, col_idx))
    return row_idx, col_idx


def check_clean_computation(rowref, colref, mytuple, signal_spread_data):
    """
    Check quality of reference pixel.
    :param rowref: int
    :param colref: int
    :param mytuple: intf_tuple**  of the form that may be re-factored soon.
    :param signal_spread_data: 2D array, size matches GRD array data
    """
    ref_pixel_values = mytuple.zvalues[:, rowref, colref]
    num_nans = np.sum(np.isnan(ref_pixel_values))
    num_intfs = len(ref_pixel_values)
    assert (num_nans / num_intfs < 0.5), ValueError("DataCube refpixel has >50% nans.")
    reference_ss = signal_spread_data[rowref, colref]
    print("Intf Stack has %f percent nan igrams for ref pixel %d, %d" % (100*num_nans/num_intfs, rowref, colref))
    print("Signal Spread has %f percent coherent igrams for ref pixel %d, %d" % (reference_ss, rowref, colref))
    assert (signal_spread_data[rowref, colref] > 50), ValueError("RefPix has <50% coherent interferograms.")
    return


def report_on_refpixel(rowref, colref, signal_spread_data, outdir):
    ofile = open(outdir+"/metrics_report.txt", 'w')
    ofile.write("Refpixel is Row/col %d %d \n" % (rowref, colref))
    ofile.write("Percentage of good interferograms at that pixel: %f \n" % signal_spread_data[rowref, colref])
    ofile.close()
    return


#  ----------- Connected Components ------------  #
def find_largest_connected_component(labels):
    """ Return the number of the largest component in a connected component graph. """
    counting_elements = collections.defaultdict(int)
    for item in labels:
        counting_elements[str(int(item))] += 1
    max_key = max(counting_elements, key=counting_elements.get)
    return max_key, counting_elements[max_key]


def find_connected_dates(date_pairs, sample_date):
    """
    :param date_pairs: list of strings with date pairs used, format '2015157_2018177' (real julian day)
    :param sample_date: desired date, format '2015157'
    """
    connected_dates = []
    for i in range(len(date_pairs)):
        if date_pairs[i][0:7] == sample_date:
            connected_dates.append(date_pairs[i][8:15])
        if date_pairs[i][8:15] == sample_date:
            connected_dates.append(date_pairs[i][0:7])
    return connected_dates


def connected_components_search(date_pairs, datestrs):
    """
    Are we inverting a complete network?
    This function will catch both 'disconnected networks' and 'bad day' cases.
    We want only one connected component with len==len(datestrs).
    Otherwise the network should fail.
    1. initialize the queue with the first date
    2. find anything connected to that date, add to the queue.
    3. When all of the dates connected to the first date have been added to queue, pop the first date from queue
    4. Repeat until the queue is gone.
    5. See if all dates have been added to the component.
    Returns the number of the biggest connected component, the number of elements of that component, and
    the total labels, for each date.
    Reduces to a single connected component with length len(datestrs) if we have one cc that touches every date.
    :param date_pairs: list of strings with date pairs used, format '2015157_2018177' (real julian day)
    :param datestrs: list of strings with dates desired for inversion, format '2015157'
    """

    # Initializing first time through
    label = np.zeros(np.shape(datestrs))
    unenqueued = set(datestrs)
    current_cc = 1

    while unenqueued:   # while we still have unlabeled nodes:
        # Put the first datestr with cc==0 into the queue, for traversal. Ugly but workable.
        queue = [unenqueued.pop()]

        # loop into the next connected component
        while len(queue) > 0:
            sample_date = queue.pop()  # go exploring the queue, starting at the back.
            idx_sample = datestrs.index(sample_date)
            label[idx_sample] = current_cc
            connected_dates = find_connected_dates(date_pairs, sample_date)
            for i in range(len(connected_dates)):
                if connected_dates[i] in unenqueued:
                    queue.append(connected_dates[i])
                    unenqueued.remove(connected_dates[i])

        current_cc = current_cc + 1
    # print(label)  # testing code
    cc_num, num_elements = find_largest_connected_component(label)
    return cc_num, num_elements, label


def reduce_graph_to_largest_cc(date_pairs, datestrs):
    """
    Shrink a network to only largest connected component of the graph
    :param date_pairs: list of strings with date pairs used, format '2015157_2018177' (real julian day)
    :param datestrs: list of strings with dates desired for inversion, format '2015157'
    """
    cc_num, num_elements, labels = connected_components_search(date_pairs, datestrs)
    if num_elements == len(datestrs):
        return date_pairs, datestrs
    else:
        new_date_pairs, new_datestrs = [], []
        for i in range(len(datestrs)):
            if labels[i] == float(cc_num):
                new_datestrs.append(datestrs[i])
        for item in date_pairs:
            if item[0:7] in new_datestrs:
                new_date_pairs.append(item)
        return new_date_pairs, new_datestrs


# Functions to get TS points in row/col coordinates
def drive_cache_ts_points(ts_points_file, intf_file_example, geocoded_flag):
    """ If you want to re-compute things, you need to delete the cache. """
    if os.path.isfile(ts_points_file + ".cache"):
        lons, lats, names, rows, cols = read_ts_points_file(ts_points_file + ".cache")
    elif os.path.isfile(ts_points_file):  # if there's no cache, we will make one.
        lons, lats, names, rows, cols = read_ts_points_file(ts_points_file)
        lons, lats, names, rows, cols = match_ts_points_row_col(lons, lats, names, rows, cols,
                                                                intf_file_example, geocoded_flag)
        write_ts_points_file(lons, lats, names, rows, cols, ts_points_file + ".cache")
    else:
        print(
            "Error! You ask for points but there's no ts_points_file %s . No points computed. " % ts_points_file)
        return None, None, None, None, None
    return lons, lats, names, rows, cols


def match_ts_points_row_col(lons, lats, names, rows, cols, example_grd, geocoded_flag):
    """Find each row and col that hasn't been found before, either in geocoded or radarcoords. """
    trans_dat = "merged/trans.dat"
    for i in range(len(lons)):
        if rows[i] == '':
            if geocoded_flag:
                irow, icol = get_reference_pixel_from_geocoded_grd(lons[i], lats[i], example_grd)
            else:
                irow, icol = get_referece_pixel_from_radarcoord_grd(lons[i], lats[i], trans_dat, example_grd)
            rows[i] = irow
            cols[i] = icol
    return lons, lats, names, rows, cols


def read_ts_points_file(ts_points_file):
    """
    Here we can use several formats simultaneously. Point name is required.
    Format 1:  -117.76 35.88 313 654 coso1
    Format 2:  -117.76 35.90 coso2
    """
    print("Reading file %s" % ts_points_file)
    lons, lats, names, rows, cols = [], [], [], [], []
    ifile = open(ts_points_file, 'r')
    for line in ifile:
        temp = line.split()
        if len(temp) == 3:  # we have provided the lat/lon/name
            lons.append(float(temp[0]))
            lats.append(float(temp[1]))
            names.append(temp[2])
            rows.append('')
            cols.append('')
        elif len(temp) == 5:  # we have provided the lat/lon/row/col/name
            if np.isnan(float(temp[2])):  # if the cache has nan for those pixels, we skip. 
                continue
            else:
                lons.append(float(temp[0]))
                lats.append(float(temp[1]))
                rows.append(int(temp[2]))
                cols.append(int(temp[3]))
                names.append(temp[4])
        else:
            print("Format of input ts_points file not recognized. ")
            sys.exit(0)
    print("  Computing time series at %d geographic points " % (len(lons)))
    ifile.close()
    return lons, lats, names, rows, cols


# Output Functions
def write_ts_points_file(lons, lats, names, rows, cols, ts_points_file):
    print("Writing %s " % ts_points_file)
    ofile = open(ts_points_file, 'w')
    for i in range(len(lons)):
        ofile.write("%.5f %.5f %s %s %s\n" % (lons[i], lats[i], str(rows[i]), str(cols[i]), names[i]))
    ofile.close()
    return


def get_axarr_numbers(cols, idx):
    """Given an incrementally counting idx number and a subplot dimension, where is our plot?
    total_plots = rows * cols """
    col_num = np.mod(idx, cols)
    row_num = int(np.floor(idx / cols))
    return row_num, col_num


def plot_full_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1):
    """ Make a nice time series plot. """
    tdata, xdata, ydata, TS_array = netcdf_read_write.read_3D_netcdf(TS_NC_file)
    num_rows_plots = 3
    num_cols_plots = 4

    f, axarr = plt.subplots(num_rows_plots, num_cols_plots, figsize=(16, 10), dpi=300)
    for i in range(len(xdates)):
        rownum, colnum = get_axarr_numbers(num_cols_plots, i)
        axarr[rownum][colnum].imshow(TS_array[i, :, :], aspect=aspect, cmap='rainbow', vmin=vmin, vmax=vmax)
        titlestr = dt.datetime.strftime(xdates[i], "%Y-%m-%d")
        axarr[rownum][colnum].get_xaxis().set_visible(False)
        axarr[rownum][colnum].set_title(titlestr, fontsize=20)

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical')
    cb.set_label('Displacement (mm)', fontsize=18)
    cb.ax.tick_params(labelsize=12)

    plt.savefig(TS_image_file)
    return


def plot_incremental_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1):
    """Make a nice incremental displacement time series plot. """
    tdata, xdata, ydata, TS_array = netcdf_read_write.read_3D_netcdf(TS_NC_file)
    num_rows_plots = 3
    num_cols_plots = 4

    # Combining the two shortest intervals into one. 
    print(np.shape(TS_array))
    selected = [0, 1, 3, 4, 5, 6, 7, 8, 9, 10]
    TS_array = TS_array[selected, :, :]
    xdates = [xdates[i] for i in range(11) if i in selected]
    print(np.shape(TS_array))
    print(len(xdates))

    f, axarr = plt.subplots(num_rows_plots, num_cols_plots, figsize=(16, 10), dpi=300)
    for i in range(1, len(xdates)):
        rownum, colnum = get_axarr_numbers(num_cols_plots, i)
        data = np.subtract(TS_array[i, :, :], TS_array[i - 1, :, :])
        axarr[rownum][colnum].imshow(data, aspect=aspect, cmap='rainbow', vmin=vmin, vmax=vmax)
        titlestr = dt.datetime.strftime(xdates[i], "%Y-%m-%d")
        axarr[rownum][colnum].get_xaxis().set_visible(False)
        axarr[rownum][colnum].set_title(titlestr, fontsize=20)

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical')
    cb.set_label('Displacement (mm)', fontsize=18)
    cb.ax.tick_params(labelsize=12)

    plt.savefig(TS_image_file)
    return
