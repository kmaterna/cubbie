# Code to take a set of interferograms that span a particular earthquake and generate an average. 
# The average should contain less noise than the original interferograms.

import numpy as np
import sys
import readmytupledata as rmd
import netcdf_read_write as rwr


def drive_coseismic_stack_gmtsar(intf_files, wavelength, rowref, colref, outdir):
    intf_tuple = rmd.reader(intf_files);
    average_coseismic = get_avg_coseismic(intf_tuple, rowref, colref, wavelength);
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, 'mm',
                              outdir + '/coseismic.grd');
    rwr.produce_output_plot(outdir + '/coseismic.grd', 'LOS Displacement', outdir + '/coseismic.png',
                            'displacement (mm)');
    return;


def drive_coseismic_stack_isce(intf_files, wavelength, rowref, colref, outdir):
    intf_tuple = rmd.reader_isce(intf_files);
    average_coseismic = get_avg_coseismic(intf_tuple, rowref, colref, wavelength);
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, 'mm',
                              outdir + '/coseismic.grd');
    rwr.produce_output_plot(outdir + '/coseismic.grd', 'LOS Displacement', outdir + '/coseismic.png',
                            'displacement (mm)',
                            aspect=1 / 8, invert_yaxis=False, vmin=-50, vmax=200);
    return;


def get_avg_coseismic(intf_tuple, rowref, colref, wavelength):
    # I could send this into the iterator_func in NSBAS if I wanted to.
    # Negative sign matches the NSBAS code
    disp = np.zeros(np.shape(intf_tuple.zvalues[0, :, :]));
    for i in range(len(intf_tuple.yvalues)):
        for j in range(len(intf_tuple.xvalues)):
            pixel_value = np.subtract(intf_tuple.zvalues[:, i, j], intf_tuple.zvalues[:, rowref, colref]);
            disp[i][j] = np.nanmean(pixel_value) * -wavelength / (4 * np.pi);
    return disp;
