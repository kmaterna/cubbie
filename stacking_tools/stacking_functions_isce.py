import sys
import numpy as np
import read_write_insar_utilities.netcdf_plots
import readmytupledata
import stacking_utilities
import stack_corr
from Tectonic_Utils.read_write import netcdf_read_write as rwr
from intf_generating import isce_geocode_tools, unwrapping_isce_custom
import stacking_configparser


def custom_isce_unwrapping(config_params):
    custom_params = stacking_configparser.read_config_isce(config_params.config_file)
    unwrapping_isce_custom.main_function(custom_params.rlks, custom_params.alks, custom_params.filt,
                                         custom_params.xbounds, custom_params.ybounds,
                                         custom_params.cor_cutoff_mask);
    intf_file_tuples = stacking_utilities.get_list_of_intf_all(config_params);
    corr_files = [x[3] for x in intf_file_tuples];
    stack_corr.drive_signal_spread_isce(corr_files, 0.5, config_params.ts_output_dir, "signalspread_full.nc");
    return;


def stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, ts_output_dir, label=""):
    # WE MAKE THE SIGNAL SPREAD FOR THE CUT IMAGES
    cor_files = [i.replace("fully_processed.unwrappedphase", "cut.cor") for i in intf_files];  # get for isce
    netcdfname = ts_output_dir + '/signalspread_cut_ref' + label + '.nc'
    cor_value = np.nan;
    cor_data = readmytupledata.reader_isce(cor_files);
    a = stack_corr.stack_corr(cor_data, cor_value);
    rwr.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', netcdfname)
    read_write_insar_utilities.netcdf_plots.produce_output_plot(netcdfname, 'Signal Spread above cor=' + str(cor_value),
                                                                ts_output_dir + '/signalspread_cut_ref' + label + '.png', 'Percentage of coherence',
                                                                aspect=1 / 4, invert_yaxis=False, dot_points=[[colref], [rowref]]);
    signal_spread_ref = a[rowref, colref];
    print("Signal Spread of the reference pixel = %.2f " % signal_spread_ref);
    if signal_spread_ref < 50:
        print("WARNING: Your reference pixel has very low coherence. Consider picking a different one.");
        print("STOPPING ON PURPOSE.");
        sys.exit(0);
    return;


def geocode_isce_uavsar(config_params):
    geocode_directory = config_params.ts_output_dir + "/isce_geocode";
    # Deleting the contents of this folder would be a good automatic step in the future.
    isce_geocode_tools.gmtsar_nc_stack_2_isce_stack(config_params.ts_output_dir + "/TS.nc", geocode_directory,
                                                    bands=2);  # write the TS data into isce binaries
    W, E, S, N = isce_geocode_tools.geocode_UAVSAR_stack(config_params,
                                                         geocode_directory);  # do this once or more than once
    isce_geocode_tools.create_isce_stack_unw_geo(geocode_directory, W, E, S, N);
    isce_geocode_tools.create_isce_stack_rdr_geo(geocode_directory, W, E, S, N);
    isce_geocode_tools.inspect_isce(geocode_directory);
    return;
