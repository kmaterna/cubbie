# Code to take a set of interferograms that span a particular earthquake and generate an average. 
# The average should contain less noise than the original interferograms.

import numpy as np
import readmytupledata as rmd
import netcdf_read_write as rwr


def drive_coseismic_stack(config_params, intf_files):
    # Can be used for gmtsar or isce
    param_dict = get_coseismic_params(config_params);
    intf_tuple = param_dict["reader"](intf_files);
    average_coseismic = get_avg_coseismic(intf_tuple, param_dict["rowref"], param_dict["colref"],
                                          param_dict["wavelength"]);
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, 'mm',
                              param_dict["outdir"]+'/coseismic.grd');
    rwr.produce_output_plot(param_dict["outdir"]+'/coseismic.grd', 'LOS Displacement',
                            param_dict["outdir"]+'/coseismic.png', 'displacement (mm)');
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


def get_coseismic_params(config_params):
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
