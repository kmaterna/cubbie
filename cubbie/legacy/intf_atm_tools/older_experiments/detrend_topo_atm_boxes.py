"""
A more sophisticated topo-vs.-atmosphere software
Splits the area into chunks, and solves for a linear trend in each box.
Plots the variance and slope.
These may change with distance from the coast or elevation.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from tectonic_utils.read_write.netcdf_read_write import read_netcdf3


def main_function(input_dir='intf_all', outdir='atm_topo'):
    [filename, demfile, rowsample, colsample] = configure(input_dir, outdir)
    [topo, _, _, zdata] = inputs(filename, demfile)
    [slope_array] = local_compute(topo, zdata, rowsample, colsample)
    outputs(zdata, topo, slope_array, outdir)
    return


# ------ CONFIGURE THE MAJOR LOOP ------------ #
def configure(input_dir, outdir):
    # input_file=input_dir+'2016196_2016220'+'/'+'phasefilt.grd'
    input_file = os.path.join(input_dir, '2016244_2016268', 'phasefilt.grd')
    print("Detrending atm/topo file %s" % input_file)
    demfile = os.path.join('topo', 'topo_ra.grd')
    os.makedirs(outdir, exist_ok=True)
    rowsample = 50
    colsample = 30
    return [input_file, demfile, rowsample, colsample]


def inputs(inputfile, demfile):
    [_, _, topo] = read_netcdf3(demfile)
    [xdata, ydata, zdata] = read_netcdf3(inputfile)
    topo = np.flipud(topo)
    return [topo, xdata, ydata, zdata]


# ------ COMPUTE ------------ #
def local_compute(topo, zdata, rowsample, colsample):
    rowdim, coldim = np.shape(zdata)  # rowdim = 1524. coldim = 661.
    num_row_iterations = int(np.ceil(rowdim / float(rowsample)))
    num_col_iterations = int(np.ceil(coldim / float(colsample)))

    slope_array = np.zeros((num_row_iterations, num_col_iterations))

    for i in range(num_row_iterations):
        for j in range(num_col_iterations):

            zarray = []  # the 1D array with original phase values
            demarray = []  # the 1D array with topography
            startrow = i * rowsample
            startcol = j * colsample

            # Collect valid phase values
            for m in range(startrow, min(startrow + rowsample, rowdim)):
                for n in range(startcol, min(startcol + colsample, coldim)):
                    if ~np.isnan(zdata[m][n]) and zdata[m][n] != 0.000:
                        zarray.append(zdata[m][n])
                        demarray.append(topo[m][n])

            # Now generate a best-fitting slope between phase and topography

            # Here we need an adjustmenet for phase jumps. What is the appropriate slope?

            sea_level = 20  # meters  If the span of elevation is less than this, we call it ocean.
            if not demarray:
                coef = [np.nan, np.nan]
            elif max(demarray) - min(demarray) <= sea_level:
                print("we found an element with no topography: %d, %d " % (i, j))
                coef = [np.nan, np.nan]
            else:
                coef = np.polyfit(demarray, zarray, 1)

            slope_array[i][j] = coef[0]

            # Making a plot
            if j == num_col_iterations / 2:  # and i<min(num_row_iterations, num_col_iterations):
                print('editing slope_array %d' % i)
                if not demarray:
                    continue  # cannot make plot for empty array.
                f, axarr = plt.subplots(1, 3)
                axarr[0].plot(demarray, zarray, '.')
                xvals = np.arange(min(demarray), max(demarray), 1)
                axarr[0].plot(xvals, [(x - min(demarray)) * coef[0] for x in xvals], '--r')
                axarr[0].set_ylim([-np.pi, np.pi])
                axarr[1].imshow(zdata[startrow:startrow + rowsample, startcol:startcol + colsample], cmap='jet')
                axarr[2].imshow(topo[startrow:startrow + rowsample, startcol:startcol + colsample], cmap='gray')
                plt.savefig(os.path.join('atm_topo', 'testbox' + str(i) + '.eps'))
                plt.close()

    return [slope_array]


def outputs(zdata, topo, slope_array, outdir):
    f, axarr = plt.subplots(1, 3)
    axarr[0].imshow(zdata, cmap='jet')
    axarr[1].imshow(topo, cmap='gray')
    ax2 = axarr[2].imshow(slope_array, cmap='hsv')
    plt.colorbar(ax2)
    plt.savefig(os.path.join(outdir, 'testbox.eps'))
    plt.close()
    return


if __name__ == "__main__":
    main_function()
