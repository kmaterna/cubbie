#!/usr/bin/python
import numpy as np
import glob

import read_write_insar_utilities.netcdf_plots as netcdf_plots
import readmytupledata as rmd
from Tectonic_Utils.read_write import netcdf_read_write


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
    if np.isnan(cutoff):
        # for GMTSAR, we usually use this criterion
        # (cutoff has been imposed during unwrapping, and the bad pixels are already nans)
        a = 100 * np.sum(~np.isnan(data_vector)) / len(data_vector);
    else:
        a = 100 * len(np.where(data_vector > cutoff)[0]) / len(data_vector);  # This has been tested.
    return a;


def dummy_signal_spread(intfs, output_dir, output_filename):
    # Make a perfect signal spread for passing to other applications
    print("Making a dummy signal spread that matches interferograms' dimensions (perfect 100).");
    output_filename = output_dir + "/" + output_filename;
    [xdata, ydata, zdata] = netcdf_read_write.read_netcdf4(intfs[0]);
    a = np.add(np.zeros(np.shape(zdata)), 100);
    netcdf_read_write.produce_output_netcdf(xdata, ydata, a, 'Percentage', output_filename, dtype=np.float32)
    netcdf_plots.produce_output_plot(output_filename, 'Signal Spread', output_dir + '/signalspread.png',
                                     'Percentage of coherence (out of ' + str(len(intfs)) + ' images)', aspect=1.2);
    return;


def signal_spread_to_mask(ss_file, cutoff, mask_file):
    # Given a signal spread file, make a nice mask that we can use for plotting.
    [xdata, ydata, zdata] = netcdf_read_write.read_netcdf3(ss_file);
    mask_response = np.zeros(np.shape(zdata));
    for i in range(len(ydata)):
        for j in range(len(xdata)):
            if zdata[i][j] >= cutoff:
                mask_response[i][j] = 1;
            else:
                mask_response[i][j] = np.nan;
    netcdf_read_write.produce_output_netcdf(xdata, ydata, mask_response, 'unitless', mask_file);
    return;


def drive_signal_spread_calculation(corr_files, cutoff, output_dir, output_filename):
    print("Making stack_corr")
    output_file = output_dir + "/" + output_filename
    mytuple = rmd.reader(corr_files)  
    a = stack_corr(mytuple, cutoff)  # if unwrapped files, we use Nan to show when it was unwrapped successfully.
    netcdf_read_write.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', output_file)
    netcdf_plots.produce_output_plot(output_file, 'Signal Spread', output_dir + '/signalspread.png',
                                     'Percentage of coherence (out of ' + str(len(corr_files)) + ' images)',
                                     aspect=1.2);
    return;


def drive_signal_spread_isce(corr_files, cutoff, output_dir, output_filename):
    cor_data = rmd.reader_isce(corr_files);
    a = stack_corr(cor_data, cutoff);
    netcdf_read_write.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', output_dir+'/' +
                                            output_filename);
    netcdf_plots.produce_output_plot(output_dir + '/' + output_filename, 'Signal Spread above cor=' + str(cutoff),
                                     output_dir + '/signalspread_full.png', 'Percentage of coherence', aspect=1 / 4,
                                     invert_yaxis=False);
    return;


if __name__ == "__main__":
    myfiles = glob.glob("intf_all_remote/???????_???????/corr.grd")
    mytuple = rmd.reader(myfiles)
    a = stack_corr(mytuple, 0.1)
    netcdf_read_write.produce_output_netcdf(mytuple.xvalues, mytuple.yvalues, a, 'Percentage', 'signalspread.nc')
    netcdf_plots.produce_output_plot('signalspread.nc', 'Signal Spread', 'signalspread.png', 'Percentage of coherence')
