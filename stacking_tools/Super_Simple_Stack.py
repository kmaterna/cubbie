#!/usr/bin/python
import numpy as np
import glob
import netcdf_read_write as rwr
import readmytupledata as rmd
import scipy.io.netcdf as netcdf
import matplotlib.pyplot as plt
import sys

def configure():
    ramps = "Metadata/Ramp_need_fix.txt"
    outfile_stem = "Stacking/Simple_Stack/Ionosphere_corrected/"
    remove_ramp = 1
    myfiles = glob.glob("intf_all_remote/???????_???????/unwrap_ref.grd")
    myfiles_no_ramp = glob.glob("intf_all_remote/???????_???????/unwrap_ref_corrected.grd")
    return ramps, outfile_stem, myfiles, myfiles_no_ramp, remove_ramp

def inputs(ramps, myfiles, myfiles_no_ramp, remove_ramp):
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

def velocity_simple_stack(filepathslist, wavelength, manual_exclude, signal_threshold):
    """This function takes in a list of files that contain arrays of phases and times. It
    will compute the velocity of each pixel using the given wavelength of the satellite.
    Finally, it will return a 2D array of velocities, ready to be plotted. For the manual exclude
    argument, enter either 0 (no images excluded), 1 (15 images excluded), or 2 (40 images excluded). The
    final argument should be a number between 0 and 100 inclusive that tells the function which pixels
    to exclude based on this signal percentage."""

    print(signal_threshold)
    if manual_exclude != 0:
        f = open('Metadata/manual_remove.txt', 'r')
        if manual_exclude == 1:
            content, x = f.readlines()[0:15], []
            for i in range(len(content)):
                x.append(content[i].strip('\n'))
        if manual_exclude == 2:
            content = f.read()
            x = content.split('\n')
        f.close()
        filesmodified = []
        filepathslist = filesmodified
        for i in range(len(myfiles_new)):
            if myfiles_new[i][16:31] not in x:
                filesmodified.append(myfiles_new[i])
    print('Number of files being stacked: ' + str(len(filepathslist)))
    signal_spread_data=rwr.read_grd("signalspread_please_test.nc")  # This will have to be refactored.... 
    mytuple = rmd.reader(filepathslist)
    phases, times = [], []
    velocities = np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)))
    i,j,f,c = 0,0,0,0
    for z in np.nditer(mytuple.zvalues, order='F'):
        if np.isnan(z) == False:
            if signal_spread_data[i,j] < signal_threshold:
                times.append(np.nan)
            else:
                phases.append(mytuple.zvalues[f][i][j])
                times.append(mytuple.date_deltas[f])
        if np.isnan(z) == True:
            times.append(np.nan)
        f+=1
        if f == len(mytuple.zvalues):
            velocities[i,j] = (wavelength/(4*(np.pi)))*((np.sum(phases))/(np.sum(times)))
            phases, times = [], []
            c+=1
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues)*len(mytuple.yvalues)) + ' pixels')
            f=0
            j+=1
            if j==len(mytuple.xvalues):
                j=0
                i+=1
                if i == len(mytuple.yvalues):
                    i=0
    return velocities, mytuple.xvalues, mytuple.yvalues



if __name__ == "__main__":
    ramps, outfile_stem, myfiles, myfiles_no_ramp, remove_ramp = configure()
    myfiles_new = inputs(ramps,myfiles, myfiles_no_ramp, remove_ramp)
    velocities, x, y = velocity_simple_stack(myfiles_new, 56, 1, 50 )
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', outfile_stem + 'velo_prof_reasonable50_remastered.grd')
    rwr.produce_output_plot(outfile_stem + 'velo_prof_reasonable50_remastered.grd', 'Velocity Profile Reasonable (15 images removed)', outfile_stem + 'velo_prof_reasonable50_remastered.png', 'velocity (mm/yr)')
