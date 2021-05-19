# This is in Python

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import datetime as dt
import stacking_utilities
import dem_error_correction


# ------------ UTILITY FUNCTIONS ------------ #
def get_TS_dates(date_julstrings):
    """" Get the x axis associated with a certain set of interferograms
    Takes a list of N date_julstrings in format [YYYYJJJ_YYYJJJ,...]
    Returns lists of N long associated with each acquisition: a string in format YYYYJJJ,
    a dt object, and the number of days since first image. """
    dates_total = [];
    for i in range(len(date_julstrings)):
        dates_total.append(date_julstrings[i][0:7])
        dates_total.append(date_julstrings[i][8:15])
    datestrs = sorted(set(dates_total));
    x_axis_datetimes = [dt.datetime.strptime(x, "%Y%j") for x in datestrs];
    x_axis_days = [(x - x_axis_datetimes[0]).days for x in
                   x_axis_datetimes];  # number of days since first acquisition.
    return datestrs, x_axis_datetimes, x_axis_days;


def select_datestrs_for_pixel(i, j, param_dict, intf_tuple, signal_spread_tuple, coh_tuple, datestrs):
    """
    reduce the datestrs and x_axis_days for a certain pixel based on where its interferograms have nans.
    This reduces the likelihood of encountering a singular matrix during time series SBAS inversion.
    Update the datestrs and x_axis_days
    It guarantees that the pixel's network will fall within a single connected component for SBAS inversion.
    """
    valid_date_julstrings, select_datestrs, select_x_axis_days = [], [], [];
    ss, pixel_value, coh_value = pixel_extractor(i, j, param_dict, intf_tuple, signal_spread_tuple, coh_tuple);

    # Filter the interferograms for ones that contain real data.
    for i in range(len(pixel_value)):
        if coh_value:   # if we are using coherence
            if coh_value[i] > param_dict["signal_coh_cutoff"]:
                if not math.isnan(pixel_value[i]):
                    valid_date_julstrings.append(intf_tuple.date_pairs_julian[i]);
        else:
            if not math.isnan(pixel_value[i]):
                valid_date_julstrings.append(intf_tuple.date_pairs_julian[i]);

    # Here we filter interferograms again based on the largest connected component of the graph:
    valid_date_julstrings, _ = stacking_utilities.reduce_graph_to_largest_cc(valid_date_julstrings, datestrs);

    select_datestrs, _, select_x_axis_days = get_TS_dates(valid_date_julstrings);
    return select_datestrs, select_x_axis_days;


def initial_defensive_programming(intf_tuple, signal_spread_tuple, coh_tuple, param_dict):
    assert(np.shape(intf_tuple.zvalues[0]) == np.shape(signal_spread_tuple)), ValueError("SS and intf diff shapes");
    if coh_tuple is not None and np.shape(intf_tuple.zvalues[0]) != np.shape(coh_tuple.zvalues[0]):
        print("ERROR: coherence data does not match input data. Stopping immediately. ");
        print("Shape of coherence data:", np.shape(coh_tuple.zvalues[0]));
        print("Shape of data array:", np.shape(intf_tuple.zvalues[0]));
        sys.exit(1);
    if param_dict["ts_type"] == 'WNSBAS' and param_dict["dem_error"] == 1:
        print("Error! You cannot have weighted least squares and dem_error at the same time.");
        print("Please un-set one of these options. ");
        sys.exit(1);
    stacking_utilities.check_clean_computation(param_dict["rowref"], param_dict["colref"], intf_tuple,
                                               signal_spread_tuple);
    return;


def compute_velocity_math(TS, x_axis_days):
    vel = np.polyfit(x_axis_days, TS, 1);
    vel = vel[0] * 365.24;  # conversion from mm/day to mm/yr
    return vel;


# ------------ COMPUTE ------------ #
# The point here is to loop through each pixel, determine if there's enough data to use, and then
# make an NSBAS matrix describing each image that's a real number (not nan).

def Velocities(param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple):
    """
    Solve velocities directly for each pixel, since different pixels may have different datestrs in
    regions with bad coherence.
    """
    initial_defensive_programming(intf_tuple, signal_spread_tuple, coh_tuple, param_dict);
    retval_main = np.zeros([len(intf_tuple.yvalues), len(intf_tuple.xvalues)]);
    retval_metrics = [[{} for _i in range(len(intf_tuple.xvalues))] for _j in range(len(intf_tuple.yvalues))];
    datestrs, x_dts, x_axis_days = get_TS_dates(intf_tuple.date_pairs_julian);

    def packager_function(i, j, intf_tuple):
        # Giving access to all these variables
        return compute_vel(i, j, param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple,
                           datestrs);

    retval_main, retval_metrics = iterator_func(intf_tuple, packager_function, retval_main, retval_metrics,
                                                0, 200000);   # numbers are a test.
    return retval_main, retval_metrics;


def Full_TS(param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple):
    """ This is how you access Time Series solutions from NSBAS"""
    initial_defensive_programming(intf_tuple, signal_spread_tuple, coh_tuple, param_dict);
    datestrs, x_dts, _ = get_TS_dates(intf_tuple.date_pairs_julian);
    # Establishing the return array
    empty_vector = [np.empty(np.shape(datestrs))];
    retval_main = [[empty_vector for _i in range(len(intf_tuple.xvalues))] for _j in range(len(intf_tuple.yvalues))];
    retval_metrics = [[{} for _i in range(len(intf_tuple.xvalues))] for _j in range(len(intf_tuple.yvalues))];

    def packager_function(i, j, intf_tuple):
        # Giving access to all these variables.
        return compute_TS(i, j, param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple, datestrs);

    retval_main, retval_metrics = iterator_func(intf_tuple, packager_function, retval_main, retval_metrics,
                                                param_dict["start_index"], param_dict["end_index"]);
    return retval_main, retval_metrics;


def Velocities_from_TS(ts_tuple):
    """
    The easy function to take a timeseries saved on disk and construct velocities
    This one doesn't have a memory leak.
    """
    retval_main = np.zeros([len(ts_tuple.yvalues), len(ts_tuple.xvalues)]);
    retval_metrics = [[{} for _i in range(len(ts_tuple.xvalues))] for _j in range(len(ts_tuple.yvalues))];
    x_axis_days = [(i - ts_tuple.ts_dates[0]).days for i in ts_tuple.ts_dates];

    def packager_function(i, j, ts_tuple):
        return compute_velocity_from_ts(i, j, ts_tuple, x_axis_days);

    retval_main, retval_metrics = iterator_func(ts_tuple, packager_function, retval_main, retval_metrics);
    return retval_main;


def iterator_func(intf_tuple, func, retval, retval_metrics, start_index=0, end_index=None):
    """ This iterator performs a for loop. It assumes the return value can be stored in an array of ixj"""
    print("Performing iteration on %d files" % (len(intf_tuple.zvalues)));
    print("Started at: ");
    print(dt.datetime.now());
    previous_time = dt.datetime.now();
    c = 0;
    true_count = 1;
    if end_index is None:
        end_index = len(intf_tuple.yvalues) * len(intf_tuple.xvalues);
    # iterate through the 3D array of data
    it = np.nditer(intf_tuple.zvalues[0, :, :], flags=['multi_index'], order='F');
    while not it.finished:
        i = it.multi_index[0];
        j = it.multi_index[1];
        if c >= start_index:
            if c == end_index:
                break;
            retval[i][j], nanflag, retval_metrics[i][j] = func(i, j, intf_tuple);
            if np.mod(c, 10000) == 0:
                print('Done with ' + str(c) + ' out of ' + str(
                    len(intf_tuple.xvalues) * len(intf_tuple.yvalues)) + ' pixels')
                print("  working on pixel %d %d " % (i, j));
            if not nanflag:
                true_count = true_count + 1;  # how many pixels were actually inverted?
            if np.mod(true_count, 10000) == 0:
                current_time = dt.datetime.now();
                delta = current_time - previous_time;
                print("--> 10K inversions took: %.2f s" % delta.total_seconds());
                previous_time = current_time;
        # if c==60000:
        # 	break;
        c = c + 1;
        it.iternext();
    print("Finished at: ");
    print(dt.datetime.now());
    return retval, retval_metrics;


# ---------- LOWER LEVEL COMPUTE FUNCTIONS ---------- #
# Functions that go into the iterator
def compute_vel(i, j, param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple, datestrs):
    """
    For a given pixel, what is the velocity?  We will compute SBAS time series (in mm)
    But we do not guarantee that each pixel will have the same dates in the time series that we use to make velocity.
    Velocity is just computed over the dates that are available.
    A useful caveat in places with poor coherence that are likely to have disconnected SBAS networks.
    Right now only operates on the largest connected component of the graph.
    """
    # for each pixel, update datestrs based on pixel_value and intf_tuple.
    updated_datestrs, updated_x_axis_days = select_datestrs_for_pixel(i, j, param_dict, intf_tuple, signal_spread_tuple,
                                                                      coh_tuple, datestrs);

    TS, nanflag, output_metrics_dict = compute_TS(i, j, param_dict, intf_tuple, signal_spread_tuple, baseline_tuple,
                                                  coh_tuple, updated_datestrs);
    if not nanflag:
        vel = compute_velocity_math(TS[0], updated_x_axis_days);   # the velocity step
    else:
        vel = np.nan;
    return vel, nanflag, output_metrics_dict;


def compute_velocity_from_ts(i, j, ts_tuple, x_axis_days):
    """ What is the velocity of bunch of dates and displacements? ts_tuple in mm; velocity in mm/yr """
    ts_values = ts_tuple.zvalues[:, i, j];
    if sum(np.isnan(ts_values)) > 30:   # under certain degenerate conditions, nanflag=1
        return np.nan, 1, {};
    else:
        vel = compute_velocity_math(ts_values, x_axis_days);
    return vel, 0, {};


def compute_TS(i, j, param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple, datestrs):
    """
    For a given pixel, what is the SBAS time series (in mm)?
    Apply various corrections
    """
    empty_vector = np.empty(np.shape(datestrs));  # Length of the TS model
    empty_vector[:] = np.nan;
    output_metrics_dict = {};
    ss, pixel_value, coh_value = pixel_extractor(i, j, param_dict, intf_tuple, signal_spread_tuple, coh_tuple);
    if ss > param_dict["nsbas_good_perc"] and sum(np.isnan(pixel_value)) < len(pixel_value) * 0.5:
        # nan condition is for degenerate cases (happened on coastlines where the water was just coherent enough)

        # Produce an uncorrected time series
        ts_vector = do_nsbas_pixel(pixel_value, intf_tuple.date_pairs_julian, param_dict["wavelength"], datestrs,
                                   coh_value);

        # Applying corrections
        if param_dict["dem_error"]:  # If we are implementing a DEM error correction.
            ts_vector, Kz_error = dem_error_correction.driver(ts_vector, datestrs, baseline_tuple);
            output_metrics_dict["Kz_error"] = Kz_error;
        if param_dict["sbas_smoothing"] > 0:  # Smoothing after the time series has been created
            ts_vector = temporal_smoothing_ts(ts_vector, param_dict["sbas_smoothing"]);

        TS = [ts_vector];
        if np.sum(np.isnan(TS[0])) == len(TS[0]):
            nanflag = True;
        else:
            nanflag = False;
    else:
        TS = [empty_vector];
        nanflag = True;
        if param_dict["dem_error"]:
            output_metrics_dict["Kz_error"] = np.nan;
    # Here, I have determined that every single return has a Kz_error key.
    return TS, nanflag, output_metrics_dict;


# Functions at or near the lowest level
def pixel_extractor(i, j, param_dict, intf_tuple, signal_spread_tuple, coh_tuple):
    """ Extract a pixel from several 2D arrays, referencing it to the reference pixel """
    ss = signal_spread_tuple[i, j];
    pixel_value = intf_tuple.zvalues[:, i, j];
    reference_pixel_value = intf_tuple.zvalues[:, param_dict["rowref"], param_dict["colref"]];
    if coh_tuple is None:
        coh_value = None;
    else:
        coh_value = coh_tuple.zvalues[:, i, j];
    pixel_value = np.subtract(pixel_value, reference_pixel_value);  # with respect to the reference pixel.
    return [ss, pixel_value, coh_value];


def do_nsbas_pixel(pixel_value, date_pairs, wavelength, datestrs, coh_value=None):
    """"
    pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram
    date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image,
      format 2015157_2018177 (real julian day)
    datestrs: a list of the dates we want to invert on, in format 2015157
    This solves Gm = d for the movement of the pixel with smoothing.
    If coh_value is an array, we do weighted least squares
    This function expects the values in the preferred reference system (i.e. reference pixel already implemented).
    """
    d = np.array([]);
    diagonals = [];
    date_pairs_used = [];

    for i in range(len(pixel_value)):
        if date_pairs[i][0:7] in datestrs:   # if the interferogram falls within the desired connected component
            if not math.isnan(pixel_value[i]):
                d = np.append(d, pixel_value[i]);  # removes the nans from the computation.
                date_pairs_used.append(date_pairs[i]);
                # might be a slightly shorter array of which interferograms actually got used.
                if coh_value is not None:
                    diagonals.append(np.power(coh_value[i], 2));  # using coherence squared as the weighting.
                else:
                    diagonals.append(1);

    model_num = len(datestrs) - 1;
    # The weighting vector
    W = np.diag(diagonals);

    # print('d:', len(d), d);    # testing code
    # print('date_pairs_used:', len(date_pairs_used));   # testing code

    # More defensive programming for degenerate cases like disconnected networks
    cc_num, num_elements, _ = stacking_utilities.connected_components_search(date_pairs_used, datestrs);
    if num_elements <= 4:
        print("VERY SMALL DATA MATRIX ENCOUNTERED. RETURNING VECTOR OF NANS.");
        empty_vector = np.empty(np.shape(datestrs));  # Length of the TS model
        empty_vector[:] = np.nan;
        return empty_vector;
    if num_elements != len(datestrs):
        print("SINGULAR MATRIX ENCOUNTERED. RETURNING VECTOR OF NANS.");
        empty_vector = np.empty(np.shape(datestrs));  # Length of the TS model
        empty_vector[:] = np.nan;
        return empty_vector;

    # building G matrix line by line.
    G = np.zeros([len(date_pairs_used), model_num]);
    for i in range(len(d)):
        ith_intf = date_pairs_used[i];
        first_image = ith_intf.split('_')[0];  # in format '2017082'
        second_image = ith_intf.split('_')[1];  # in format '2017094'
        first_index = datestrs.index(first_image);
        second_index = datestrs.index(second_image);
        for j in range(second_index - first_index):
            G[i][first_index + j] = 1;

    # solving the SBAS linear least squares equation for displacement between each epoch.
    if coh_value is not None:
        GTWG = np.dot(np.transpose(G), np.dot(W, G))
        GTWd = np.dot(np.transpose(G), np.dot(W, d))
        m = np.dot(np.linalg.inv(GTWG), GTWd)
    else:
        m = np.linalg.lstsq(G, d, rcond=None)[0];  # rcond=None comes from futurewarning

    # Adding up all the displacement.
    m_cumulative = [0];
    for i in range(1, len(m) + 1):
        m_cumulative.append(np.sum(m[0:i]));  # The cumulative phase from start to finish!

    # Conversion from radians to mm
    disp_ts = [i * wavelength / (4 * np.pi) for i in m_cumulative];

    # Convert from range change to subsidence (negative means moving away); set beginning to zero
    disp_ts = [i * -1 for i in disp_ts];  
    disp_ts = [i-disp_ts[0] for i in disp_ts];

    return disp_ts;


def temporal_smoothing_ts(LOS_phase, smoothing):
    """Implementing temporal smoothing after uncorrected timeseries formation.
    I'm doing this as an overconstrained linear inverse problem
    This is similar to a Gaussian smoothing."""
    n_TS = len(LOS_phase);
    G_top = np.eye(n_TS);  # Constructing the G matrix
    alpha_array = np.full((n_TS - 1,), smoothing)
    G_positive_alpha = np.diag(alpha_array);
    G_negative_alpha = np.diag(-alpha_array[0:-1], 1);
    G_bottom = np.add(G_positive_alpha, G_negative_alpha);
    G_last_column = np.zeros((n_TS - 1, 1));
    G_last_column[-1] = -smoothing;
    G_bottom = np.hstack((G_bottom, G_last_column))
    G = np.vstack((G_top, G_bottom))

    # Constructing the data vector
    d = np.hstack((LOS_phase, np.zeros(n_TS - 1, )));

    # Implementing the smoothing
    smoothed_los = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T), d);
    return smoothed_los;


# ------------ OUTPUTS ------------ #

def nsbas_ts_points_outputs(dts, m_cumulative, row, col, name, lon, lat, outdir):
    """ This outdir is expected to be something like "stacking/nsbas/ts". """

    mean_disp = np.nanmean(m_cumulative);
    plotting_ts = [i - mean_disp for i in m_cumulative];

    plt.figure();
    plt.plot(dts, plotting_ts, 'b.');
    plt.xlabel("Time");
    plt.ylabel("Displacement (mm)");
    plt.title(str(row) + ' ' + str(col) + ' ' + str(lon) + ' ' + str(lat) + ' ' + str(name));
    # plt.ylim([-40,50]);
    plt.savefig(outdir + '/' + str(name) + '_' + str(lon) + '_' + str(lat) + '_disp.eps');

    ofile = open(outdir + '/' + str(name) + '_' + str(row) + '_' + str(col) + '_record.txt', 'w');
    print("Writing file %s " % outdir + '/' + str(name) + '_' + str(row) + '_' + str(col) + '_record.txt');
    for i in range(len(dts)):
        ofile.write("%s %f %f %d %d " % (name, lon, lat, row, col));
        ofile.write(dt.datetime.strftime(dts[i], "%Y-%m-%d"));
        ofile.write(" %f\n" % (m_cumulative[i]));
    ofile.close();
    return;
