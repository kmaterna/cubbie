# A script to take a stack of images plus a DEM
# Solve for best-fitting linear trend
# Right now this is a global trend across the whole scene. 
# Remove trend and save the adjusted stack in out_dir
# In order to work correctly, this script needs a specified reference pixel. 

import numpy as np
import matplotlib.pyplot as plt
import glob, sys, os
import subprocess
import datetime as dt
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3


def main_function(staging_directory, outdir, rowref, colref, starttime, endtime):
    [filenames, demfile] = configure(staging_directory, outdir, starttime, endtime);
    demdata = subsample_read_dem(filenames[0], demfile);

    for item in filenames:
        [xdata, ydata, zdata] = read_netcdf3(item);
        [corrected_zdata, zarray, corrarray, demarray] = global_compute_item(zdata, demdata, rowref, colref);
        output_item(xdata, ydata, zdata, corrected_zdata, zarray, corrarray, demarray, item, outdir);
    return;


# ------ CONFIGURE THE MAJOR LOOP ------------ # 
def configure(staging_directory, outdir, starttime, endtime):
    print("Detrending atm/topo for all files in %s and storing the result in %s " % (staging_directory, outdir))
    file_names = glob.glob(staging_directory + "/*_*_unwrap.grd");
    valid_file_names = [];
    if len(file_names) == 0:
        print("Error! No files matching search pattern within " + staging_directory);
        sys.exit(1);

    # Only read in the files within the starttime and endtime boundaries.
    start_dt = dt.datetime.strptime(str(starttime), "%Y%m%d");
    end_dt = dt.datetime.strptime(str(endtime), "%Y%m%d");

    for ifile in file_names:
        pairname = ifile.split('/')[-1][0:15];
        image1 = pairname.split('_')[0];
        image2 = pairname.split('_')[1];
        image1_dt = dt.datetime.strptime(image1, "%Y%j");
        image2_dt = dt.datetime.strptime(image2, "%Y%j");

        if start_dt <= image1_dt <= end_dt:
            if start_dt <= image2_dt <= end_dt:
                valid_file_names.append(ifile);

    subprocess.call(['mkdir', '-p', outdir], shell=False);
    demfile = 'topo/topo_ra.grd';
    return [valid_file_names, demfile];


# -------- PREPROCESS I/O FUNCTIONS ---------- # 
def subsample_read_dem(samplefile, demfile):
    # Take an example interferogram and subsample the topo_ra to exactly match this file size.
    # You may have to play with -T/-r options in GMT grdsample to force the same gridcell/gridline registration
    # You also may have to force the netcdf4 file into netcdf3 for later reading in python.

    if not os.path.isfile(demfile):
        print("ERROR! %s does not exist- exiting " % demfile);
        sys.exit(1);

    subsampled_file = 'topo/topo_ra_subsampled.grd';
    intervals = subprocess.check_output(['gmt', 'grdinfo', '-I', samplefile], shell=False);
    intervals = intervals.split('\n')[0];
    ranges = subprocess.check_output(['gmt', 'grdinfo', '-I-', samplefile], shell=False);
    ranges = ranges.split('\n')[0];
    command = 'gmt grdsample ' + demfile + ' -Gtopo/temp.grd -T ' + intervals + ' ' + ranges;
    print(command);
    subprocess.call(command, shell=True);
    subprocess.call('nccopy -k classic topo/temp.grd ' + subsampled_file, shell=True);
    subprocess.call(['rm', 'topo/temp.grd'], shell=False);

    [_, _, z] = read_netcdf3(subsampled_file);
    return z;


# ---------- COMPUTE FUNCTIONS ------------ #

def global_compute_item(zdata, demdata, rowref, colref):
    zarray = [];  # the 1D array with original phase values
    demarray = [];  # the 1D array with topography
    corrarray = [];  # the 1D array with corrected phase

    # Collect valid phase values
    rowdim, coldim = np.shape(zdata);
    for i in range(rowdim):
        for j in range(coldim):
            if ~np.isnan(zdata[i][j]):
                zarray.append(zdata[i][j]);
                demarray.append(demdata[i][j]);

    # Now generate a best-fitting slope between phase and topography
    coef = np.polyfit(demarray, zarray, 1);

    # Now remove the slope (pinning everything so that the reference pixel has phase 0)
    corrected_zdata = np.zeros(np.shape(zdata));
    rowdim, coldim = np.shape(zdata);
    reference_pixel_offset = zdata[rowref][colref] - coef[0] * demdata[rowref][colref];
    for i in range(rowdim):
        for j in range(coldim):
            if ~np.isnan(zdata[i][j]):
                corrected_zdata[i][j] = zdata[i][j] - coef[0] * demdata[i][j] - reference_pixel_offset;
                corrarray.append(zdata[i][j] - coef[0] * demdata[i][j] - reference_pixel_offset);
            else:
                corrected_zdata[i][j] = np.nan;

    # print(corrected_zdata[rowref][colref]);

    return [corrected_zdata, zarray, corrarray, demarray];


# --------------- OUTPUTS ---------------- #
def output_item(xdata, ydata, zdata, corrected_zdata, zarray, corrarray, demarray, item, outdir):
    item_name = item.split('/')[-1];
    item_name_short = item_name.split('unwrap.grd')[0];
    netcdfname = outdir + '/' + item_name;
    netcdf_read_write.produce_output_netcdf(xdata, ydata, corrected_zdata, 'unwrapped_phase', netcdfname);
    netcdf_read_write.flip_if_necessary(netcdfname);

    plt.figure();
    plt.plot(demarray, zarray, '.');
    plt.plot(demarray, corrarray, '.r', alpha=0.15)
    plt.ylabel('phase');
    plt.xlabel('topo');
    plt.title('Initial and Corrected Phase vs. Topography')
    plt.savefig('intf_all/atm_topo_corrected.grd' + '/phase_topo_' + item_name_short + '.png');
    plt.close();

    plt.figure();
    plt.subplot(1, 2, 1);
    plt.imshow(zdata, cmap='hsv', vmin=np.min(zarray), vmax=np.max(zarray));
    plt.title('Before Correction');
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.subplot(1, 2, 2);
    plt.imshow(corrected_zdata, cmap='hsv', vmin=np.min(zarray), vmax=np.max(zarray));
    plt.title('After Correction');
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    cb = plt.colorbar();
    cb.set_label("Unwrapped Phase (radians)", size=12);
    plt.savefig(outdir + '/imshow_' + item_name_short + '.png');
    plt.close();

    return;


if __name__ == "__main__":
    main_function('intf_all/unwrap.grd', 'intf_all/atm_topo_corrected.grd', 241, 175);
