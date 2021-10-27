# Code to take a set of interferograms that span a particular earthquake and generate an average. 
# The average should contain less noise than the original interferograms.
# Can be used with gmtsar or isce

import numpy as np
from read_write_insar_utilities import netcdf_plots
from . import readmytupledata as rmd
from Tectonic_Utils.read_write import netcdf_read_write as rwr


def drive_coseismic_stack(config_params, intf_files):
    param_dict = get_coseismic_params(config_params);
    intf_tuple = param_dict["reader"](intf_files);
    average_coseismic = get_avg_coseismic(intf_tuple, param_dict["rowref"], param_dict["colref"],
                                          param_dict["wavelength"]);
    output_manager_coseismic(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, param_dict["outdir"]);
    return;


def get_avg_coseismic(intf_tuple, rowref, colref, wavelength):
    """I could send this into the iterator_func in NSBAS if I wanted to.
    Negative sign matches the NSBAS code """
    disp = np.zeros(np.shape(intf_tuple.zvalues[0, :, :]));
    for i in range(len(intf_tuple.yvalues)):
        for j in range(len(intf_tuple.xvalues)):
            pixel_value = np.subtract(intf_tuple.zvalues[:, i, j], intf_tuple.zvalues[:, rowref, colref]);
            disp[i][j] = np.nanmean(pixel_value) * -wavelength / (4 * np.pi);
    return disp;


def get_coseismic_params(config_params):
    """Unpack the parameter object into a new parameter dictionary"""
    rowref = int(config_params.ref_idx.split('/')[0]);
    colref = int(config_params.ref_idx.split('/')[1]);
    if config_params.file_format == 'isce':  # Working with the file formats
        my_reader_function = rmd.reader_isce;
    else:
        my_reader_function = rmd.reader;
    param_dictionary = {"wavelength": config_params.wavelength,
                        "rowref": rowref, "colref": colref, "outdir": str(config_params.ts_output_dir),
                        "signal_spread_filename": config_params.ts_output_dir+'/'+config_params.signal_spread_filename,
                        "reader": my_reader_function};
    return param_dictionary;


def output_manager_coseismic(x, y, average_coseismic, outdir):
    rwr.produce_output_netcdf(x, y, average_coseismic, 'mm', outdir+'/coseismic.grd');
    netcdf_plots.produce_output_plot(outdir+'/coseismic.grd', 'LOS Displacement', outdir+'/coseismic.png', 'disp (mm)');
    return;
