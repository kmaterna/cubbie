#!/usr/bin/env python
"""
The equivalent of gmt grdinfo for isce files, including displaying the number of nans.
"""

import numpy as np
from s1_batches.read_write_insar_utilities import isce_read_write
import matplotlib.pyplot as plt
import argparse

help_message = "Get information about an ISCE binary data file. \nUsage: " \
               "isceinfo.py --data_file phase.int"


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--data_file', type=str,
                   help='''filename for raster information, REQUIRED''', required=True)
    p.add_argument('-r', '--region',
                   help='''Analyze a region or subset; default is the whole image.''')
    p.add_argument('-v', '--vmin',
                   help='''Minimum color for quick-and-dirty plot.''', default=None)
    p.add_argument('-w', '--vmax',
                   help='''Maximum color for quick-and-dirty plot.''', default=None)
    p.add_argument('-p', '--plot_data',
                   help='''Whether to plot data or not''', default=True)
    exp_dict = vars(p.parse_args())
    return exp_dict


def main_body(args):
    sourcefile = args['data_file']

    if '.unw' in sourcefile:
        xdata, ydata, data = isce_read_write.read_scalar_data(sourcefile)  # reading second band
    else:
        xdata, ydata, data = isce_read_write.read_scalar_data(sourcefile)

    print("Range of X axis: xmin:", np.nanmin(xdata), "xmax:", np.nanmax(xdata),
          ", xlen:", len(xdata), ", xinc:", xdata[1]-xdata[0])
    print("Range of Y axis: ymin:", np.nanmin(ydata), "ymax:", np.nanmax(ydata),
          ", ylen:", len(ydata), ", yinc:", ydata[1]-ydata[0])
    print("Type of data: ", type(data[0][0]))
    print("Shape of data: ", np.shape(data))
    num_values = np.size(data)
    print("Number of values: ", num_values)
    print("Spread of data: zmin:", np.nanmin(data), ", zmax:", np.nanmax(data))
    num_nans = np.sum(np.isnan(data))
    percent_nans = (num_nans / num_values) * 100
    print("Number of Nans: ", num_nans, " out of ", num_values, "(", percent_nans, " percent)")

    if args['plot_data']:
        plt.figure()
        vmin = args['vmin'] if args['vmin'] is not None else np.nanmin(data)
        vmax = args['vmax'] if args['vmax'] is not None else np.nanmin(data)
        plt.imshow(data, vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.show()
    return


if __name__ == "__main__":
    args_dict = parse_arguments()
    main_body(args_dict)
