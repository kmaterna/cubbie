import numpy as np
import glob
from subprocess import call
import read_write_insar_utilities.netcdf_plots
from Tectonic_Utils.read_write import netcdf_read_write

# This is for when you've run a large SBAS in chunks of several million pixels each
# Because it saves time to run in parallel.
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3


def get_input_dirs():
    input_files = [];
    input_files.append("/Volumes/Ironwolf/Track_71/stacking/no_smoothing/0_3500000");
    input_files.append("/Volumes/Ironwolf/Track_71/stacking/no_smoothing/3500000_7000000");
    return input_files;


def get_datestrs():
    files = glob.glob("/Volumes/Ironwolf/Track_71/stacking/no_smoothing/0_3500000/*.grd");
    datestrs = [i.split('/')[-1][0:8] for i in files];
    print(datestrs);
    return datestrs;


def combine_all_files(datestr, input_dirs, output_dir):
    print("\nCombining files for date %s" % datestr);

    filename = input_dirs[0] + "/" + datestr + ".grd"
    xdata, ydata, zdata0 = read_netcdf3(filename);
    filename1 = input_dirs[1] + "/" + datestr + ".grd"
    xdata, ydata, zdata1 = read_netcdf3(filename1);
    zdata_total = np.zeros(np.shape(zdata0));

    for j in range(len(ydata)):
        if np.mod(j, 200) == 0:
            print(j)
        for k in range(len(xdata)):
            vector = [zdata0[j][k],
                      zdata1[j][k]];  # , zdata2[j][k], zdata3[j][k], zdata4[j][k], zdata5[j][k], zdata6[j][k] ];
            zdata_total[j][k] = np.sum(vector);
    output_file = output_dir + "/" + datestr + ".grd";
    output_plot = output_dir + "/" + datestr + ".png";
    netcdf_read_write.produce_output_netcdf(xdata, ydata, zdata_total, "mm", output_file);
    read_write_insar_utilities.netcdf_plots.produce_output_plot(output_file, datestr, output_plot, "mm", aspect=1.0, invert_yaxis=True,
                                                                vmin=-50, vmax=100);
    return;


if __name__ == "__main__":
    output_dir = "/Volumes/Ironwolf/Track_71/stacking/no_smoothing/combined/"
    call(["mkdir", "-p", output_dir], shell=False);
    input_dirs = get_input_dirs();
    datestrs = get_datestrs();
    for i in range(len(datestrs)):
        combine_all_files(datestrs[i], input_dirs, output_dir);

    # Then, quickly geocode all the time series files.
    for i in range(len(datestrs)):
        call(["quick_geocode.csh", "stacking/no_smoothing/combined", "merged", datestrs[i] + ".grd",
              datestrs[i] + "_ll"], shell=False);
