import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import sys
import os
import datetime as dt
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3


def top_level_driver(skip_file=None):
    [file_names, outdir, num_plots_x, num_plots_y] = configure()
    [xdata, ydata, data_all, date_pairs, skip_intfs] = inputs(file_names, skip_file)
    make_plots(xdata, ydata, data_all, date_pairs, outdir, num_plots_x, num_plots_y, skip_intfs)
    return


# ------------- CONFIGURE ------------ # 
def configure():
    file_dir = "intf_all"
    file_type = "phasefilt.grd"
    outdir = os.path.join('phasefilt', '')

    os.makedirs(outdir, exist_ok=True)

    file_names = glob.glob(os.path.join(file_dir, "*", file_type))
    if len(file_names) == 0:
        print("Error! No files matching search pattern.")
        sys.exit(1)
    print("Reading "+str(len(file_names))+" files.")
    num_plots_x = 4
    num_plots_y = 3
    return [file_names, outdir, num_plots_x, num_plots_y]


# ------------- INPUTS ------------ # 
def inputs(file_names, skip_file):
    try:
        filename = file_names[0]
        [xdata, ydata] = read_netcdf3(filename)[0:2]  # can read either netcdf3 or netcdf4.
    except TypeError:
        [xdata, ydata, _] = netcdf_read_write.read_netcdf4(file_names[0])
    data_all = []
    date_pairs = []

    file_names = sorted(file_names)  # To force into date-ascending order.

    for ifile in file_names:  # Read the data
        try:
            data = read_netcdf3(ifile)[2]
        except TypeError:
            [_, _, data] = netcdf_read_write.read_netcdf4(ifile)
        data_all.append(data)
        pairname = ifile.split('/')[-2][0:15]
        date_pairs.append(pairname)  # returning something like '2016292_2016316' for each intf
        print(pairname)

    skip_intfs = []
    if skip_file is not None:
        ifile = open(skip_file, 'r')
        for line in ifile:
            skip_intfs.append(line.split()[0])
        ifile.close()

    return [xdata, ydata, data_all, date_pairs, skip_intfs]


def make_plots(_xdata, _ydata, data_all, date_pairs, outdir, num_plots_x, num_plots_y, skip_intfs):

    for i in range(len(data_all)):
        if np.mod(i, num_plots_y*num_plots_x) == 0:
            count = i

            fignum = i/(num_plots_y*num_plots_x)  # counting figures up 0 to 1 to 2....

            # Looping forward and plotting the next 12 plots...
            f, axarr = plt.subplots(num_plots_y, num_plots_x, figsize=(10, 10))
            for k in range(num_plots_y):
                for m in range(num_plots_x):
                    if count == len(data_all):
                        break

                    # How many days separate this interferogram?
                    day1 = date_pairs[count].split('_')[0]
                    day2 = date_pairs[count].split('_')[1]
                    if day1[4:7] == "000":
                        day1 = day1[0:6]+"1"
                    if day2[4:7] == "000":
                        day2 = day2[0:6]+"1"
                    dt1 = dt.datetime.strptime(day1, '%Y%j')
                    dt2 = dt.datetime.strptime(day2, '%Y%j')
                    deltat = dt2-dt1
                    daysdiff = deltat.days

                    # The actual plotting
                    axarr[k][m].imshow(data_all[count], cmap='jet', aspect=0.5)
                    axarr[k][m].invert_yaxis()
                    axarr[k][m].invert_xaxis()
                    axarr[k][m].get_xaxis().set_ticks([])
                    axarr[k][m].get_yaxis().set_ticks([])
                    if str(date_pairs[count]) in skip_intfs:
                        axarr[k][m].set_title(str(date_pairs[count])+'   '+str(daysdiff)+' days', fontsize=8,
                                              color='red', fontweight='bold')
                    else:
                        axarr[k][m].set_title(str(date_pairs[count])+'   '+str(daysdiff)+' days', fontsize=8,
                                              color='black')

                    count = count+1
            plt.savefig(outdir+"selected_data_"+str(int(fignum))+".eps")
            plt.close()
    return


if __name__ == "__main__":
    top_level_driver()
