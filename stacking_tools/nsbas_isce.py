from subprocess import call
import numpy as np
import sys
import stacking_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 
import nsbas


def drive_velocity_nsbas(swath, intfs, nsbas_min_intfs, sbas_smoothing, wavelength, outdir):

    # SETUP. 
    signal_spread_file='F'+swath+'/'+outdir+"/signalspread_cut.nc"
    intf_tuple = rmd.reader_isce(intfs); 
    xdates = stacking_utilities.get_xdates_from_intf_tuple(intf_tuple);
    nsbas.make_stack_corr_custom(intf_tuple, signal_spread_file);  # for safety, let's make signalspread again. 
    signal_spread_data=rwr.read_grd(signal_spread_file);

    # Select a single pixel and look at its time series. 
    pixel_coords = [550, 150];
    nsbas.single_pixel_ts(intf_tuple, pixel_coords, sbas_smoothing, wavelength, signal_spread_file);

    # TIME SERIES
    TS, xdates = nsbas.compute_fullTS_nsbas(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data);
    TS_NC_file = 'F'+swath+'/'+outdir+"/TS.nc";
    TS_image_file = 'F'+swath+'/'+outdir+"/TS.png";
    rwr.produce_output_timeseries(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', TS_NC_file);
    stacking_utilities.plot_full_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1/4);

    # AVERAGE VELOCITIES
    # velocities = nsbas.compute_nsbas(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data); 
    # rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr', 'F'+swath+'/'+outdir+'/velo_nsbas.grd');
    # rwr.produce_output_plot('F'+swath+'/'+outdir+'/velo_nsbas.grd', 'LOS Velocity',
    #     'F'+swath+'/'+outdir+'/velo_nsbas.png', 'velocity (mm/yr)',aspect=1/4, invert_yaxis=False);

    return;
