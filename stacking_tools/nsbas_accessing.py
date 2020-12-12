from subprocess import call
import numpy as np
import glob
import stacking_utilities
import readmytupledata as rmd
from read_write_insar_utilities import netcdf_read_write as rwr
import nsbas
import dem_error_correction
import sbas_testing


def reader_function_gmtsar(intf_files, coh_files, baseline_file, ts_type, dem_error):
    # A massive reader function for SBAS analysis
    if ts_type == 'WNSBAS':
        coh_tuple = rmd.reader(coh_files);
    else:
        coh_tuple = None;
    if dem_error:
        baseline_tuple = stacking_utilities.read_baseline_table(baseline_file);
    else:
        baseline_tuple = None;
    intf_tuple = rmd.reader(intf_files);
    return intf_tuple, coh_tuple, baseline_tuple;


def reader_function_isce(intf_files, coh_files, baseline_file, ts_type, dem_error):
    # A massive reader function for SBAS analysis
    intf_tuple = rmd.reader_isce(intf_files);
    if ts_type == 'WNSBAS':
        coh_tuple = rmd.reader_isce(coh_files);
    else:
        coh_tuple = None;
    if dem_error:
        baseline_tuple = stacking_utilities.read_baseline_table(baseline_file);
    else:
        baseline_tuple = None;
    return intf_tuple, coh_tuple, baseline_tuple;


def repack_param_dictionary(config_params):
    rowref = int(config_params.ref_idx.split('/')[0]);
    colref = int(config_params.ref_idx.split('/')[1]);
    if config_params.file_format == 'isce':  # Working with the file formats
        my_reader_function = reader_function_isce;
    else:
        my_reader_function = reader_function_gmtsar;
    param_dictionary = {"nsbas_good_perc": config_params.nsbas_min_intfs,
                        "sbas_smoothing": config_params.sbas_smoothing, "wavelength": config_params.wavelength,
                        "rowref": rowref, "colref": colref, "ts_output_dir": config_params.ts_output_dir,
                        "signal_spread_filename": config_params.ts_output_dir+'/'+config_params.signal_spread_filename,
                        "dem_error": config_params.dem_error, "ts_type": config_params.ts_type,
                        "reader": my_reader_function,
                        "baseline_file": config_params.baseline_file, "geocoded_flag": config_params.geocoded_intfs};
    return param_dictionary;


def write_output_metrics(param_dict, intf_tuple, metrics):
    # Unpack the dictionary that contains output metrics (if any), write into files.
    if param_dict["dem_error"]:
        gridshape = np.shape(metrics);
        Kz_grid = np.zeros(gridshape);
        for i in range(gridshape[0]):
            for j in range(gridshape[1]):
                if "Kz_error" in metrics[i][j].keys():
                    Kz_grid[i][j] = metrics[i][j]["Kz_error"];
        rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, Kz_grid, 'm',
                                  param_dict["ts_output_dir"] + '/kz_error.grd');
        rwr.produce_output_plot(param_dict["ts_output_dir"] + '/kz_error.grd', 'DEM Error',
                                param_dict["ts_output_dir"] + '/kz_error.png', 'DEM Error (m)');
    return;


def nsbas_ts_format_selector(config_params, intf_files, corr_files):
    # This is the function called from the top level coordinator
    param_dictionary = repack_param_dictionary(config_params);
    if config_params.ts_format == 'velocity':
        drive_velocity(param_dictionary, intf_files, corr_files);
    elif config_params.ts_format == 'points':
        drive_point_ts(param_dictionary, intf_files, corr_files, config_params.ts_points_file);
    elif config_params.ts_format == 'timeseries':
        drive_full_TS(param_dictionary, intf_files, corr_files);
    elif config_params.ts_format == 'velocities_from_timeseries':
        make_vels_from_ts_grids(config_params.ts_output_dir, geocoded=config_params.geocoded_intfs);
    else:
        print("Error!");
    return;


# LET'S GET A VELOCITY FIELD FROM GMTSAR FILES
def drive_velocity(param_dict, intf_files, coh_files):
    intf_tuple, coh_tuple, baseline_tuple = param_dict["reader"](intf_files, coh_files, param_dict["baseline_file"],
                                                                 param_dict["ts_type"], param_dict["dem_error"]);
    [_, _, signal_spread_tuple] = rwr.read_any_grd_xyz(param_dict["signal_spread_filename"]);
    velocities = nsbas.Velocities(param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple);
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr',
                              param_dict["ts_output_dir"] + '/velo_nsbas.grd');
    rwr.produce_output_plot(param_dict["ts_output_dir"] + '/velo_nsbas.grd', 'LOS Velocity',
                            param_dict["ts_output_dir"] + '/velo_nsbas.png', 'velocity (mm/yr)');
    return;


# LET'S GET THE FULL TS FOR EVERY PIXEL
def drive_full_TS(param_dict, intf_files, coh_files):
    param_dict["start_index"] = 0;
    param_dict["end_index"] = 7000000;
    intf_tuple, coh_tuple, baseline_tuple = param_dict["reader"](intf_files, coh_files, param_dict["baseline_file"],
                                                                 param_dict["ts_type"], param_dict["dem_error"]);
    [_, _, signal_spread_tuple] = rwr.read_any_grd_xyz(param_dict["signal_spread_filename"]);
    TS, metrics = nsbas.Full_TS(param_dict, intf_tuple, signal_spread_tuple, baseline_tuple, coh_tuple);
    rwr.produce_output_TS_grids(intf_tuple.xvalues, intf_tuple.yvalues, TS, intf_tuple.ts_dates, 'mm',
                                param_dict["ts_output_dir"]);
    write_output_metrics(param_dict, intf_tuple, metrics);
    return;


# LET'S GET SOME PIXELS AND OUTPUT THEIR TS. 
def drive_point_ts(param_dict, intf_files, coh_files, ts_points_file):
    # Replicating what would happen for a single pixel in the main SBAS loop
    # For general use, please provide a file with [lon, lat, row, col, name]
    lons, lats, names, rows, cols = stacking_utilities.drive_cache_ts_points(ts_points_file, intf_files[0],
                                                                             param_dict["geocoded_flag"]);
    outdir = param_dict["ts_output_dir"] + "/ts";
    print("TS OUTPUT DIR IS: " + outdir);
    call(['mkdir', '-p', outdir], shell=False);
    print("Computing TS for %d pixels" % len(lons));
    intf_tuple, coh_tuple, baseline_tuple = param_dict["reader"](intf_files, coh_files, param_dict["baseline_file"],
                                                                 param_dict["ts_type"], param_dict["dem_error"]);
    datestrs, x_dts, x_axis_days = nsbas.get_TS_dates(intf_tuple.date_pairs_julian);
    reference_pixel_vector = intf_tuple.zvalues[:, param_dict["rowref"], param_dict["colref"]];

    for i in range(len(rows)):
        pixel_value = intf_tuple.zvalues[:, rows[i], cols[i]];
        pixel_value = np.subtract(pixel_value, reference_pixel_vector);  # with respect to the reference pixel. 
        coh_value = None;
        if coh_tuple is not None:
            coh_value = coh_tuple.zvalues[:, rows[i], cols[i]];
        sbas_testing.io_functions.write_testing_pixel(intf_tuple, pixel_value, coh_value,
                                                      outdir+'/test_pixel_'+str(i)+'.txt');
        ts_vector = nsbas.do_nsbas_pixel(pixel_value, intf_tuple.date_pairs_julian,
                                         param_dict["wavelength"], datestrs, coh_value=coh_value);

        if param_dict["dem_error"]:  # If we are implementing a DEM error correction.
            # Shouldn't be done on Weighted least squares.
            ts_vector, Kz_error = dem_error_correction.driver(ts_vector, datestrs, baseline_tuple);
        if param_dict["sbas_smoothing"] > 0:  # Smoothing after the time series has been created
            ts_vector = nsbas.temporal_smoothing_ts(ts_vector, param_dict["sbas_smoothing"]);

        nsbas.nsbas_ts_points_outputs(x_dts, ts_vector, rows[i], cols[i], names[i], lons[i], lats[i], outdir);
    return;


def make_vels_from_ts_grids(ts_dir, geocoded=False):
    if geocoded:
        filelist = glob.glob(ts_dir + "/publish/*_ll.grd");
        mydata = rmd.reader_from_ts(filelist, "lon", "lat", "z");  # put these if using geocoded values
    else:
        filelist = glob.glob(ts_dir + "/????????.grd");
        mydata = rmd.reader_from_ts(filelist);
    vel = nsbas.Velocities_from_TS(mydata);
    rwr.produce_output_netcdf(mydata.xvalues, mydata.yvalues, vel, 'mm/yr', ts_dir + '/velo_nsbas.grd');
    rwr.produce_output_plot(ts_dir + '/velo_nsbas.grd', 'LOS Velocity', ts_dir + '/velo_nsbas.png', 'velocity (mm/yr)');
    return;
