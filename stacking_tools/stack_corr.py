#!/usr/bin/python
import numpy as np
import glob, sys
import netcdf_read_write as rwr
import readmytupledata as rmd
import netcdf_read_write


def stack_corr(mytuple, cutoff):
    """This function takes in a mytuple of data (argument 1) and counts how many times a certain
    piece of data is above a specified cutoff value (argument 2) in each 2-D array stored in mytuple.
    It returns a 2-D array of percentages, showing how much certain pieces of data satisfy the given cutoff
    condition. You can use cutoff=np.nan to do number of non-nans"""
    print('Number of files being stacked: ' + str(len(mytuple.filepaths)));
    a = np.zeros((len(mytuple.yvalues), len(mytuple.xvalues)))
    c = 0;
    it = np.nditer(mytuple.zvalues[0, :, :], flags=['multi_index'], order='F');  # iterate through the 3D array of data
    while not it.finished:
        i = it.multi_index[0];
        j = it.multi_index[1];
        data_vector = mytuple.zvalues[:, i, j]
        a[i][j] = get_signal_spread(data_vector, cutoff);
        c = c + 1;
        if np.mod(c, 20000) == 0:
            print('Done with ' + str(c) + ' out of ' + str(len(mytuple.xvalues) * len(mytuple.yvalues)) + ' pixels')
        it.iternext();
    return a;


def get_signal_spread(data_vector, cutoff):
    # For a pixel, what is the percentage of good images? 
    if np.isnan(
            cutoff):  # for GMTSAR, we usually use this criterion (cutoff has been imposed during unwrapping, and the bad pixels are already nans)
        a = 100 * np.sum(~np.isnan(data_vector)) / len(data_vector);
    else:
        a = 100 * len(np.where(data_vector > cutoff)[0]) / len(data_vector);  # This has been tested.
    return a;


def drive_signal_spread_calculation(intfs, cutoff, output_dir, output_filename):
    print("Making stack_corr")
    output_file = output_dir + "/" + output_filename
    mytuple = rmd.reader(intfs[0:2])  # get rid of [0:2] after testing
    # a = stack_corr(mytuple, cutoff)  # if unwrapped files, we use Nan to show when it was unwrapped successfully.
    intf_file = intfs[1];
    [x, y, z] = netcdf_read_write.read_netcdf4_xyz(intf_file);
    # netcdf_read_write.produce_output_netcdf(x, y, z, 'mm', output_dir+'/python_written_file.nc');
    netcdf_read_write.write_netcdf4(x, y, z, output_dir+'/python_written_file.nc');
    sys.exit(0)

    a = np.zeros(np.shape(mytuple.zvalues[0])); # get rid of this after testing
    rwr.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', output_file)
    rwr.produce_output_plot(output_file, 'Signal Spread', output_dir + '/signalspread.png',
                            'Percentage of coherence (out of ' + str(len(intfs)) + ' images)', aspect=1.2);
    return;


if __name__ == "__main__":
    myfiles = glob.glob("intf_all_remote/???????_???????/corr.grd")
    mytuple = rmd.reader(myfiles)
    a = stack_corr(mytuple, 0.1)
    rwr.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', 'signalspread.nc')
    rwr.produce_output_plot('signalspread.nc', 'Signal Spread', 'signalspread.png', 'Percentage of coherence')
