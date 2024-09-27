#!/usr/bin/env python

"""
This is for when you've run a large SBAS in chunks of several million pixels each
Because it saves time to run in parallel.
"""

import numpy as np
import glob
import os
from subprocess import call
from ...read_write_insar_utilities import netcdf_plots
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3, produce_output_netcdf


def get_input_dirs():
    input_files = ["/Volumes/Ironwolf/Track_71/stacking/no_smoothing/0_3500000",
                   "/Volumes/Ironwolf/Track_71/stacking/no_smoothing/3500000_7000000"]
    return input_files


def get_datestrs():
    files = glob.glob("/Volumes/Ironwolf/Track_71/stacking/no_smoothing/0_3500000/*.grd")
    datestrs = [os.path.split(file)[1][0:8] for file in files]
    print(datestrs)
    return datestrs


def combine_all_files(datestr, input_dirs, output_dir):
    print("\nCombining files for date %s" % datestr)

    filename = str(os.path.join(input_dirs[0], datestr + ".grd"))
    xdata, ydata, zdata0 = read_netcdf3(filename)
    filename1 = str(os.path.join(input_dirs[1], datestr + ".grd"))
    xdata, ydata, zdata1 = read_netcdf3(filename1)
    zdata_total = np.zeros(np.shape(zdata0))

    for j in range(len(ydata)):
        if np.mod(j, 200) == 0:
            print(j)
        for k in range(len(xdata)):
            vector = [zdata0[j][k],
                      zdata1[j][k]]  # , zdata2[j][k], zdata3[j][k], zdata4[j][k], zdata5[j][k], zdata6[j][k] ]
            zdata_total[j][k] = np.sum(vector)
    output_file = os.path.join(output_dir, datestr + ".grd")
    output_plot = os.path.join(output_dir, datestr + ".png")
    produce_output_netcdf(xdata, ydata, zdata_total, "mm", output_file)
    netcdf_plots.produce_output_plot(output_file, datestr, output_plot, "mm", aspect=1.0,
                                     invert_yaxis=True, vmin=-50, vmax=100)
    return


if __name__ == "__main__":
    my_output_dir = "/Volumes/Ironwolf/Track_71/stacking/no_smoothing/combined/"
    os.makedirs(my_output_dir, exist_ok=True)
    my_input_dirs = get_input_dirs()
    my_datestrs = get_datestrs()
    for i in range(len(my_datestrs)):
        combine_all_files(my_datestrs[i], my_input_dirs, my_output_dir)

    # Then, quickly geocode all the time series files.
    for one_date in my_datestrs:
        call(["quick_geocode.csh", "stacking/no_smoothing/combined", "merged", one_date + ".grd",
              one_date + "_ll"], shell=False)
