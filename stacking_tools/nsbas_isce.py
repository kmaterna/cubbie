from subprocess import call
import numpy as np
import sys
import stacking_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 
import nsbas


def drive_nsbas(swath, intf_files, nsbas_min_intfs, sbas_smoothing, wavelength, outdir, coh_files=[]):

    # SETUP. 
    signal_spread_file=outdir+"/signalspread_cut.nc"
    intf_tuple = rmd.reader_isce(intf_files); 
    if coh_files != []:
        coh_tuple = rmd.reader_isce(coh_files);
    else:
        coh_tuple = [];
    xdates = stacking_utilities.get_xdates_from_intf_tuple(intf_tuple);
    # nsbas.make_stack_corr_custom(intf_tuple, signal_spread_file);  # for safety, let's make signalspread again. 
    signal_spread_data=rwr.read_grd(signal_spread_file);

    # Select a single pixel and look at its time series. 
    pixel_coords = [550, 150];
    nsbas.single_pixel_ts(intf_tuple, pixel_coords, sbas_smoothing, wavelength, signal_spread_file, outdir, coh_tuple);

    # # TIME SERIES
    TS = nsbas.driver_Full_TS(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data, coh_tuple);
    TS_NC_file = outdir+"/TS.nc";
    TS_image_file = outdir+"/TS.png";
    rwr.produce_output_timeseries(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', TS_NC_file);
    stacking_utilities.plot_full_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1/4);

    return;
