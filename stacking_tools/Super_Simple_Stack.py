#!/usr/bin/python
import numpy as np
import glob
import netcdf_read_write as rwr
import readmytupledata as rmd
import scipy.io.netcdf as netcdf
import matplotlib.pyplot as plt
import sys


def drive_velocity_simple_stack(swath, intfs, wavelength, outdir):
    signal_spread_data=rwr.read_grd("F"+swath+"/"+outdir+"/signalspread.nc");
    intf_tuple = rmd.reader(intfs);
    velocities, x, y = velocity_simple_stack(intf_tuple, wavelength, signal_spread_data, 25);  # signal threshold < 100%.  lower signal threshold allows for more data into the stack.  
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', 'F'+swath+'/'+outdir+'/velo_simple_stack.grd')
    rwr.produce_output_plot('F'+swath+'/'+outdir+'/velo_simple_stack.grd', 'Velocity Profile ',
        'F'+swath+'/'+outdir+'/velo_simple_stack.png', 'velocity (mm/yr)');
    return; 

def get_velocity_by_stacking(phase_values, time_intervals, wavelength):
    # The math behind the simple stack method. 
    phase_count=0; time_count=0.0001;   # we put a small number here to avoid div-by-zero during debugging. 
    for i in range(len(phase_values)):
        if not np.isnan(phase_values[i]):
            phase_count=phase_count+phase_values[i];
            time_count=time_count+time_intervals[i];
    velocity = (wavelength/(4*(np.pi)))*(phase_count/time_count);
    return velocity;

def velocity_simple_stack(mytuple, wavelength, signal_spread_data, signal_threshold):
    """This function takes in a list of files that contain arrays of phases and times. 
    It will compute the velocity of each pixel using the satellite's  wavelength. It will return a 2D array of velocities. 
    The final argument should be a number between 0 and 100 inclusive that tells the function which pixels
    to exclude based on this signal percentage."""
    print('Number of files being stacked: ' + str(len(mytuple.zvalues)));
    velocities = np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)));
    c = 0;

    it = np.nditer(mytuple.zvalues[0,:,:], flags=['multi_index'], order='F');  # iterate through the 3D array of data
    while not it.finished:
        i=it.multi_index[0];
        j=it.multi_index[1];
        signal_spread = signal_spread_data[i,j];
        if signal_spread > signal_threshold: # if we want a calculation for that day... 
            velocities[i,j]=get_velocity_by_stacking(mytuple.zvalues[:,i,j], mytuple.date_deltas, wavelength);
        else:
            velocities[i,j]=np.nan;
        c=c+1;
        if np.mod(c,10000)==0:
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues)*len(mytuple.yvalues)) + ' pixels')        
        it.iternext();
    return velocities, mytuple.xvalues, mytuple.yvalues


# --------- FOR TESTING AND MANUAL CONTROL ONLY ----------- # 

def configure():
    ramps = "Metadata/Ramp_need_fix.txt"
    outfile_stem = "Stacking/Simple_Stack/Ionosphere_corrected/"
    remove_ramp = 1
    myfiles = glob.glob("intf_all_remote/???????_???????/unwrap_ref.grd")
    myfiles_no_ramp = glob.glob("intf_all_remote/???????_???????/unwrap_ref_corrected.grd")
    return ramps, outfile_stem, myfiles, myfiles_no_ramp, remove_ramp

def inputs(ramps, myfiles, myfiles_no_ramp, remove_ramp, manual_exclude):
    f = open(ramps, 'r')
    raw, content = f.readlines()[:], []
    for i in range(len(raw)):
        content.append(raw[i].strip('\n'))
    myfiles_new = []
    for i in range(len(myfiles)):
        if remove_ramp != 0:
            test = myfiles[i].replace("ref", "ref_corrected")
            if test in myfiles_no_ramp:
                myfiles_new.append(test)
            if myfiles[i][16:31] not in content:
                myfiles_new.append(myfiles[i])
        else:
            myfiles_new.append(myfiles[i])
    return myfiles_new


if __name__ == "__main__":
    ramps, outfile_stem, myfiles, myfiles_no_ramp, remove_ramp = configure()
    myfiles_new = inputs(ramps,myfiles, myfiles_no_ramp, remove_ramp, 1)
    velocities, x, y = velocity_simple_stack(myfiles_new, 56, 50 )
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', outfile_stem + 'velo_prof_reasonable50_remastered.grd')
    rwr.produce_output_plot(outfile_stem + 'velo_prof_reasonable50_remastered.grd', 'Velocity Profile Reasonable (15 images removed)', outfile_stem + 'velo_prof_reasonable50_remastered.png', 'velocity (mm/yr)')
