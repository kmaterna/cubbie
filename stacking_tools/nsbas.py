# This is in Python

import numpy as np
import matplotlib.pyplot as plt
import sys, math
import datetime as dt
import sentinel_utilities
import stacking_utilities
import netcdf_read_write as rwr
import stack_corr


# ------------ UTILITY FUNCTIONS ------------ #
def get_TS_dates(date_julstrings):
    # Get me the x axis associated with a certain set of interferograms
    # Takes a list of N date_julstrings in format [YYYYJJJ_YYYJJJ,...]
    # Returns lists of N long associated with each acquisition: a string in format YYYYJJJ,
    #   a dt object, and the number of days since first image.
    dates_total = [];
    for i in range(len(date_julstrings)):
        dates_total.append(date_julstrings[i][0:7])
        dates_total.append(date_julstrings[i][8:15])
    datestrs = sorted(set(dates_total));
    x_axis_datetimes = [dt.datetime.strptime(x, "%Y%j") for x in datestrs];
    x_axis_days = [(x - x_axis_datetimes[0]).days for x in
                   x_axis_datetimes];  # number of days since first acquisition.
    return datestrs, x_axis_datetimes, x_axis_days;


# Before velocities, might want to make stack_corr_custom if you have excluded significant data.
def make_stack_corr_custom(mydata, signal_spread_file):
    # Stack corr for this exact calculation (given right inputs and outputs). 
    a = stack_corr.stack_corr(mydata, np.nan);
    titlestr = "Signal Spread of " + str(len(mydata.date_deltas)) + " files";
    rwr.produce_output_netcdf(mydata.xvalues, mydata.yvalues, a, 'Percentage', signal_spread_file);
    rwr.produce_output_plot(signal_spread_file, titlestr, signal_spread_file + '.png',
                            'Percentage of coherence', aspect=1 / 4, invert_yaxis=False)
    return;


# ------------ COMPUTE ------------ #
# The point here is to loop through each pixel, determine if there's enough data to use, and then
# make an NSBAS matrix describing each image that's a real number (not nan).

def Velocities(intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data, coh_tuple=None):
    # This is how you access velocity solutions from NSBAS
    retval = np.zeros([len(intf_tuple.yvalues), len(intf_tuple.xvalues)]);
    datestrs, x_dts, x_axis_days = get_TS_dates(intf_tuple.date_pairs_julian);

    def packager_function(i, j, intf_tuple):
        # Giving access to all these variables
        return compute_vel(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data,
                           datestrs, x_axis_days, coh_tuple);

    retval = iterator_func(intf_tuple, packager_function, retval);
    return retval


def Full_TS(intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data, start_index=0,
            end_index=10e6, coh_tuple=None):
    # This is how you access Time Series solutions from NSBAS
    datestrs, x_dts, x_axis_days = get_TS_dates(intf_tuple.date_pairs_julian);
    # Establishing the return array
    empty_vector = [np.empty(np.shape(x_axis_days))];
    retval = [[empty_vector for i in range(len(intf_tuple.xvalues))] for j in range(len(intf_tuple.yvalues))];

    def packager_function(i, j, intf_tuple):
        # Giving access to all these variables.
        return compute_TS(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data,
                          datestrs, x_axis_days, coh_tuple);

    retval = iterator_func(intf_tuple, packager_function, retval, start_index, end_index);
    return retval;


def Velocities_from_TS(ts_tuple):
    # The easy function to take a timeseries saved on disk and construct velocities
    # This one doesn't have a memory leak.
    retval = np.zeros([len(ts_tuple.yvalues), len(ts_tuple.xvalues)]);
    x_axis_days = [(i - ts_tuple.ts_dates[0]).days for i in ts_tuple.ts_dates];

    def packager_function(i, j, ts_tuple):
        return compute_velocity_from_ts(i, j, ts_tuple, x_axis_days);

    retval = iterator_func(ts_tuple, packager_function, retval);
    return retval;


def iterator_func(intf_tuple, func, retval, start_index=0, end_index=None):
    # This iterator performs a for loop. It assumes the return value can be stored in an array of ixj
    # if np.shape(retval) != np.shape(signal_spread_data):
    # 	print("ERROR: signal spread does not match input data. Stopping immediately. ");
    # 	print("Shape of signal spread:", np.shape(signal_spread_data));
    # 	print("Shape of data array:", np.shape(intf_tuple.zvalues[0]));
    # 	sys.exit(0);
    print("Performing NSBAS on %d files" % (len(intf_tuple.zvalues)));
    print("Started at: ");
    print(dt.datetime.now());
    previous_time = dt.datetime.now();
    c = 0;
    true_count = 1;
    if end_index is None:
        end_index = len(intf_tuple.yvalues) * len(intf_tuple.xvalues);
    it = np.nditer(intf_tuple.zvalues[0, :, :], flags=['multi_index'],
                   order='F');  # iterate through the 3D array of data
    while not it.finished:
        i = it.multi_index[0];
        j = it.multi_index[1];
        if c >= start_index:
            if c == end_index:
                break;
            retval[i][j], nanflag = func(i, j, intf_tuple);
            if np.mod(c, 10000) == 0:
                print('Done with ' + str(c) + ' out of ' + str(
                    len(intf_tuple.xvalues) * len(intf_tuple.yvalues)) + ' pixels')
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
    return retval;


# ---------- LOWER LEVEL COMPUTE FUNCTIONS ---------- #

def compute_vel(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data, datestrs,
                x_axis_days, coh_tuple=None):
    TS, nanflag = compute_TS(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref,
                             signal_spread_data, datestrs, x_axis_days, coh_tuple);
    if not nanflag:
        vel = np.polyfit(x_axis_days, TS, 1);
        vel = vel[0] * 365.24;
    else:
        vel = np.nan;
    return vel, nanflag;


def compute_TS(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, rowref, colref, signal_spread_data, datestrs,
               x_axis_days, coh_tuple=None):
    # For a given iteration, what's the right time series?
    empty_vector = np.empty(np.shape(x_axis_days));  # Length of the TS model
    empty_vector[:] = np.nan;
    signal_spread = signal_spread_data[i, j];
    pixel_value = intf_tuple.zvalues[:, i, j];
    reference_pixel_value = intf_tuple.zvalues[:, rowref, colref];
    if coh_tuple is None:
        coh_value = None;
    else:
        coh_value = coh_tuple.zvalues[:, i, j];
    if signal_spread > nsbas_good_perc:
        pixel_value = np.subtract(pixel_value, reference_pixel_value);  # with respect to the reference pixel.
        vel, vector = do_nsbas_pixel(pixel_value, intf_tuple.date_pairs_julian, smoothing, wavelength, datestrs,
                                     x_axis_days, coh_value);
        TS = [vector];
        nanflag = False;
    else:
        TS = [empty_vector];
        nanflag = True;
    return TS, nanflag;


def compute_velocity_from_ts(i, j, ts_tuple, x_axis_days):
    # What is the velocity of bunch of dates and displacements? in mm/yr.
    ts_values = ts_tuple.zvalues[:, i, j];
    if sum(np.isnan(ts_values)) > 30:
        nanflag = 1;
        vel = np.nan;
    else:
        vel = np.polyfit(x_axis_days, ts_values, 1);
        vel = vel[0] * 365.24;
        nanflag = 0;
    return vel, nanflag;


def do_nsbas_pixel(pixel_value, date_pairs, smoothing, wavelength, datestrs, x_axis_days, coh_value=None):
    # pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram, already referenced to reference.
    # date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image, in format 2015157_2018177 (real julian day, 1-offset corrected)
    # This solves Gm = d for the movement of the pixel with smoothing.
    # If coh_value is an array, we do weighted least squares.
    # This function expects the values in the preferred reference system (i.e. reference pixel already implemented).

    # Defensive programming for degenerate cases (happened on coastlines where the water was just coherent enough)
    empty_vector = np.empty(np.shape(x_axis_days));  # Length of the TS model
    empty_vector[:] = np.nan;
    if sum(np.isnan(pixel_value)) > len(pixel_value) * 0.5:
        return np.nan, empty_vector;

    d = np.array([]);
    diagonals = [];
    date_pairs_used = [];

    for i in range(len(pixel_value)):
        if not math.isnan(pixel_value[i]):
            d = np.append(d, pixel_value[i]);  # removes the nans from the computation.
            date_pairs_used.append(
                date_pairs[i]);  # might be a slightly shorter array of which interferograms actually got used.
            if coh_value is not None:
                diagonals.append(np.power(coh_value[i], 2));  # using coherence squared as the weighting.
            else:
                diagonals.append(1);
    model_num = len(datestrs) - 1;

    # The weighting vector
    W = np.diag(diagonals);

    G = np.zeros([len(date_pairs_used), model_num]);

    # building G matrix line by line.
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
        m = np.linalg.lstsq(G, d)[0];

    # modeled_data=np.dot(G,m);
    # plt.figure();
    # plt.plot(d,'.b');
    # plt.plot(modeled_data,'.--g');
    # plt.savefig('d_vs_m.eps')
    # plt.close();

    # Adding up all the displacement.
    m_cumulative = [];
    m_cumulative.append(0);
    for i in range(1, len(m) + 1):
        m_cumulative.append(np.sum(m[0:i]));  # The cumulative phase from start to finish!

    # Smoothing after the time series has been created
    if smoothing > 0:
        m_cumulative = temporal_smoothing_ts(m_cumulative, smoothing);

    # Solving for linear velocity
    x = np.zeros([len(x_axis_days), 2]);
    y = np.array([]);
    for i in range(len(x_axis_days)):
        x[i][0] = x_axis_days[i];
        x[i][1] = 1;
        y = np.append(y, [m_cumulative[i]]);
    model_slopes = np.linalg.lstsq(x, y)[0];  # units: phase per day.

    # Velocity conversion: units in mm / year
    vel = model_slopes[0];  # in radians per day
    vel = vel * wavelength * 365.24 / 2.0 / (2 * np.pi);

    disp_ts = [i * wavelength / (4 * np.pi) for i in m_cumulative];
    # plt.figure();
    # model_line = [model_slopes[1] + x * model_slopes[0] for x in x_axis_days];
    # plt.plot(x_axis_days[0:-1],m,'b.');
    # plt.plot(x_axis_days,m_cumulative,'g.');
    # plt.plot(x_axis_days, model_line,'--g');
    # plt.xlabel("days");
    # plt.ylabel("cumulative phase");
    # plt.text(0,0,str(vel)+"mm/yr slope");
    # plt.savefig('m_model.eps');

    return vel, disp_ts;


def temporal_smoothing_ts(LOS_phase, smoothing):
    # Implementing temporal smoothing after uncorrected timeseries formation.
    # I'm doing this as an overconstrained linear inverse problem
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

def nsbas_ts_outputs(dts, m_cumulative, row, col, name, lon, lat, outdir):
    # This outdir is expected to be something like "stacking/nsbas/ts".
    # It must exist before the function is called.

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
    for i in range(len(dts)):
        ofile.write("%s %f %f %d %d " % (name, lon, lat, row, col));
        ofile.write(dt.datetime.strftime(dts[i], "%Y-%m-%d"));
        ofile.write(" %f\n" % (m_cumulative[i]));
    ofile.close();
    return;
