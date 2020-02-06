#!/usr/bin/python
import numpy as np
import glob, sys
import netcdf_read_write as rwr
import readmytupledata as rmd

def stack_corr(mytuple, cutoff):
    """This function takes in a mytuple of data (argument 1) and counts how many times a certain
    piece of data is above a specified cutoff value (argument 2) in each 2-D array stored in mytuple.
    It returns a 2-D array of percentages, showing how much certain pieces of data satisfy the given cutoff
    condition. You can use cutoff=np.nan to do number of non-nans"""
    print('Number of files being stacked: ' + str(len(mytuple.filepaths)));
    a=np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)))
    c = 0;
    it = np.nditer(mytuple.zvalues[0,:,:], flags=['multi_index'], order='F');  # iterate through the 3D array of data
    while not it.finished:
        i=it.multi_index[0];
        j=it.multi_index[1];
        data_vector = mytuple.zvalues[:,i,j]
        a[i][j] = get_signal_spread(data_vector, cutoff);
        c=c+1;
        if np.mod(c,20000)==0:
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues)*len(mytuple.yvalues)) + ' pixels')
        it.iternext();
    return a;

def get_signal_spread(data_vector, cutoff):
    if np.isnan(cutoff):
        a = 100 * np.sum(~np.isnan(data_vector))/len(data_vector);
    else:
        a = 100 * np.sum(np.where(data_vector>cutoff))/len(data_vector);
    return a;

def drive_signal_spread_calculation(swath, unwrap_ref_dir, intfs, output_dir):
    print("Making stack_corr")
    mytuple=rmd.reader(intfs)
    a=stack_corr(mytuple, np.nan)
    rwr.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', 'F'+swath+'/'+output_dir+'/signalspread.nc')
    rwr.produce_output_plot('F'+swath+'/'+output_dir+'/signalspread.nc', 'Signal Spread', 
        'F'+swath+'/'+output_dir+'/signalspread.png', 'Percentage of coherence (out of '+str(len(intfs))+' images)' )
    return;


if __name__ == "__main__":
    myfiles = glob.glob("intf_all_remote/???????_???????/corr.grd")
    mytuple=rmd.reader(myfiles)
    a=stack_corr(mytuple, 0.1)
    rwr.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', 'signalspread.nc')
    rwr.produce_output_plot('signalspread.nc', 'Signal Spread', 'signalspread.png', 'Percentage of coherence' )
