#!/usr/bin/python
import numpy as np
import glob, sys
import netcdf_read_write as rwr
import readmytupledata as rmd


def drive_velocity_simple_stack(intfs, wavelength, rowref, colref, outdir):
    signal_spread_data = rwr.read_grd(outdir + "/signalspread.nc");
    intf_tuple = rmd.reader(intfs);
    velocities, x, y = velocity_simple_stack(intf_tuple, wavelength, rowref, colref, signal_spread_data, 25);
                # signal threshold < 100%.  lower signal threshold allows for more data into the stack.
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', outdir + '/velo_simple_stack.grd')
    rwr.produce_output_plot(outdir + '/velo_simple_stack.grd', 'LOS Velocity ', outdir + '/velo_simple_stack.png',
                            'velocity (mm/yr)');
    return;


def get_velocity_by_stacking_pixel(phase_values, time_intervals, wavelength):
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
    It will compute the velocity of each pixel using the satellite's  wavelength. It will return a 2D array of velocities. 
    The final argument should be a number between 0 and 100 inclusive that tells the function which pixels
    to exclude based on this signal percentage."""
    print('Number of files being stacked: ' + str(len(mytuple.zvalues)));
    velocities = np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)));
    ref_pixel_values = mytuple.zvalues[:, rowref, colref];
    c = 0;

    it = np.nditer(mytuple.zvalues[0, :, :], flags=['multi_index'], order='F');  # iterate through the 3D array of data
    while not it.finished:
        i = it.multi_index[0];
        j = it.multi_index[1];
        signal_spread = signal_spread_data[i, j];
        if signal_spread > signal_threshold:  # if we want a calculation for that day...
            pixel_value = np.subtract(mytuple.zvalues[:, i, j], ref_pixel_values);
            velocities[i, j] = get_velocity_by_stacking_pixel(pixel_value, mytuple.date_deltas, wavelength);
        else:
            velocities[i, j] = np.nan;
        c = c + 1;
        if np.mod(c, 10000) == 0:
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues) * len(mytuple.yvalues)) + ' pixels')
        it.iternext();
    return velocities, mytuple.xvalues, mytuple.yvalues


# --------- FOR TESTING AND MANUAL CONTROL ONLY ----------- # 

if __name__ == "__main__":
    print("Manual Control of Super Simple Stack.");
    # velocities, x, y = velocity_simple_stack(myfiles_new, 56, 50 )
    # rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', outfile_stem + 'velo_prof_reasonable50_remastered.grd')
    # rwr.produce_output_plot(outfile_stem + 'velo_prof_reasonable50_remastered.grd', 'Velocity Profile Reasonable (15 images removed)', outfile_stem + 'velo_prof_reasonable50_remastered.png', 'velocity (mm/yr)')
