#!/usr/bin/python
import numpy as np
from read_write_insar_utilities import netcdf_plots
from Tectonic_Utils.read_write import netcdf_read_write as rwr
import readmytupledata as rmd
import stacking_utilities


def drive_velocity_simple_stack(config_params, intf_files):
    param_dict = get_simple_stack_params(config_params);
    [_, _, signal_spread_data] = rwr.read_any_grd(param_dict["signal_spread_filename"]);
    intf_tuple = param_dict["reader"](intf_files);
    velocities, x, y = velocity_simple_stack(intf_tuple, param_dict["wavelength"], param_dict["rowref"],
                                             param_dict["colref"], signal_spread_data, 25);
    # last argument is signal threshold (< 100%).  lower signal threshold allows for more data into the stack.
    output_manager_simple_stack(x, y, velocities, param_dict["rowref"], param_dict["colref"], signal_spread_data,
                                param_dict["outdir"]);
    return;


def get_simple_stack_params(config_params):
    # repacking the parameter dictionary
    rowref = int(config_params.ref_idx.split('/')[0]);
    colref = int(config_params.ref_idx.split('/')[1]);
    if config_params.file_format == 'isce':  # Working with the file formats
        my_reader_function = rmd.reader_isce;
    else:
        my_reader_function = rmd.reader;
    param_dictionary = {"wavelength": config_params.wavelength,
                        "rowref": rowref, "colref": colref, "outdir": str(config_params.ts_output_dir),
                        "signal_spread_filename": config_params.ts_output_dir+'/'+config_params.signal_spread_filename,
                        "reader": my_reader_function};
    return param_dictionary;


def pixel_velocity_by_stacking(phase_values, time_intervals, wavelength):
    # The math behind the simple stack method. 
    phase_count = 0;
    time_count = 0.0001;  # we put a small number here to avoid div-by-zero during debugging.
    for i in range(len(phase_values)):
        if not np.isnan(phase_values[i]):
            phase_count = phase_count + phase_values[i];
            time_count = time_count + time_intervals[i];
    velocity = (wavelength / (4 * np.pi)) * (phase_count / time_count);
    return velocity;


def velocity_simple_stack(mytuple, wavelength, rowref, colref, signal_spread_data, signal_threshold):
    """This function takes in a list of files that contain arrays of phases and times. 
    It will compute the velocity of each pixel using the satellite's  wavelength. It will return 2D array of velocities.
    The final argument should be a number between 0 and 100 inclusive that tells the function which pixels
    to exclude based on this signal percentage."""
    print('Number of files being stacked: ' + str(len(mytuple.zvalues)));
    velocities = np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)));
    ref_pixel_values = mytuple.zvalues[:, rowref, colref];
    stacking_utilities.check_clean_computation(rowref, colref, mytuple, signal_spread_data);
    c = 0;

    it = np.nditer(mytuple.zvalues[0, :, :], flags=['multi_index'], order='F');  # iterate through the 3D array of data
    while not it.finished:
        i = it.multi_index[0];
        j = it.multi_index[1];
        signal_spread = signal_spread_data[i, j];
        if signal_spread > signal_threshold:  # if we want a calculation for that day...
            pixel_value = np.subtract(mytuple.zvalues[:, i, j], ref_pixel_values);
            velocities[i, j] = pixel_velocity_by_stacking(pixel_value, mytuple.date_deltas, wavelength);
        else:
            velocities[i, j] = np.nan;
        c = c + 1;
        if np.mod(c, 10000) == 0:
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues) * len(mytuple.yvalues)) + ' pixels')
        it.iternext();
    return velocities, mytuple.xvalues, mytuple.yvalues


def output_manager_simple_stack(x, y, velocities, rowref, colref, signal_spread_data, outdir):
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', outdir + '/velo_simple_stack.grd')
    netcdf_plots.produce_output_plot(outdir + '/velo_simple_stack.grd', 'LOS Velocity ', outdir +
                                     '/velo_simple_stack.png', 'velocity (mm/yr)');
    stacking_utilities.report_on_refpixel(rowref, colref, signal_spread_data, outdir);
    return;


# --------- FOR TESTING AND MANUAL CONTROL ONLY ----------- # 

if __name__ == "__main__":
    print("Manual Control of Super Simple Stack.");
