from subprocess import call
import numpy as np
import sys, glob
import datetime as dt 
import stacking_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 
import nsbas


# LET'S GET A VELOCITY FIELD
def drive_velocity_gmtsar(intf_files, nsbas_min_intfs, smoothing, wavelength, rowref, colref, outdir, coh_files=[]):
    # GMTSAR DRIVING VELOCITIES
    signal_spread_file=outdir+"/signalspread.nc"
    intf_tuple = rmd.reader(intf_files); 
    coh_tuple = [];
    if coh_files != []:
        coh_tuple = rmd.reader(coh_files);
    signal_spread_data=rwr.read_grd(signal_spread_file);
    velocities = nsbas.Velocities(intf_tuple, nsbas_min_intfs, smoothing, wavelength, rowref, colref, signal_spread_data, coh_tuple); 
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr', outdir+'/velo_nsbas.grd');
    rwr.produce_output_plot(outdir+'/velo_nsbas.grd', 'LOS Velocity', outdir+'/velo_nsbas.png', 'velocity (mm/yr)');
    return;


# LET'S GET SOME PIXELS AND OUTPUT THEIR TS. 
def drive_point_ts_gmtsar(intf_files, ts_points_file, smoothing, wavelength, rowref, colref, outdir, coh_files=[]):
    # For general use, please provide a file with [lon, lat, row, col, name]
    lons, lats, names, rows, cols = stacking_utilities.drive_cache_ts_points(ts_points_file);
    if len(lons)==0:
        return; 
    outdir=outdir+"/ts";
    print("TS OUTPUT DIR IS: " + outdir);
    call(['mkdir','-p',outdir],shell=False);
    intf_tuple = rmd.reader(intf_files);
    datestrs, x_dts, x_axis_days = nsbas.get_TS_dates(intf_tuple.dates_correct); 
    reference_pixel_vector = intf_tuple.zvalues[:,rowref,colref];
    for i in range(len(rows)):
        pixel_value = intf_tuple.zvalues[:,rows[i],cols[i]];
        pixel_value = np.subtract(pixel_value, reference_pixel_vector);  # with respect to the reference pixel. 
        vel, m_cumulative = nsbas.do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, smoothing, wavelength, datestrs, x_axis_days); 
        m_cumulative=[i*-1 for i in m_cumulative];  # My sign convention seems to be opposite to Katia's
        nsbas.nsbas_ts_outputs(x_dts, m_cumulative, rows[i], cols[i], names[i], lons[i], lats[i], outdir);
    return;


# LET'S GET THE FULL TS FOR EVERY PIXEL
def drive_full_TS_gmtsar(intf_files, nsbas_min_intfs, sbas_smoothing, wavelength, rowref, colref, outdir, coh_files=[]):
    # SETUP. 
    start_index = 0;
    end_index   = 7000000;
    signal_spread_file=outdir+"/signalspread.nc"

    intf_tuple = rmd.reader(intf_files);
    coh_tuple=[];
    if coh_files != []:
        coh_tuple = rmd.reader(coh_files);
    xdates = stacking_utilities.get_xdates_from_intf_tuple(intf_tuple);
    signal_spread_data=rwr.read_grd(signal_spread_file);

    # TIME SERIES
    TS = nsbas.Full_TS(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, rowref, colref, signal_spread_data, 
        start_index=start_index, end_index=end_index, coh_tuple=coh_tuple);
    # TS_NC_file = outdir+"/TS.nc";
    # rwr.produce_output_timeseries(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', TS_NC_file);
    rwr.produce_output_TS_grids(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', outdir);
    return;



def make_vels_from_ts(ts_dir,geocoded=False):
    if geocoded:
        filelist = glob.glob(ts_dir+"/????????_ll.grd");
        mydata = rmd.reader_from_ts(filelist,"lon","lat","z");  # put these if using geocoded values
    else:
        filelist = glob.glob(ts_dir+"/????????.grd");
        mydata = rmd.reader_from_ts(filelist); 
    vel = nsbas.Velocities_from_TS(mydata);
    rwr.produce_output_netcdf(mydata.xvalues, mydata.yvalues, vel, 'mm/yr', ts_dir+'/velo_nsbas.grd');
    rwr.produce_output_plot(ts_dir+'/velo_nsbas.grd', 'LOS Velocity', ts_dir+'/velo_nsbas.png', 'velocity (mm/yr)');
    return;




# LET'S GET THE FULL TS FOR UAVSAR/ISCE FILES. 
def drive_full_TS_isce(intf_files, nsbas_min_intfs, sbas_smoothing, wavelength, rowref, colref, outdir, coh_files=[]):
    # SETUP. 
    signal_spread_file=outdir+"/signalspread_cut.nc"
    intf_tuple = rmd.reader_isce(intf_files); 
    coh_tuple=[];
    if coh_files != []:
        coh_tuple = rmd.reader_isce(coh_files);
    xdates = stacking_utilities.get_xdates_from_intf_tuple(intf_tuple);
    signal_spread_data=rwr.read_grd(signal_spread_file);

    # TIME SERIES
    TS = nsbas.Full_TS(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, rowref, colref, signal_spread_data, coh_tuple=coh_tuple);

    # OUTPUTS
    TS_NC_file = outdir+"/TS.nc";
    TS_image_file = outdir+"/TS.png";
    rwr.produce_output_timeseries(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', TS_NC_file);
    stacking_utilities.plot_full_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1/8);
    return;



def drive_point_ts_isce():
    # Select a single pixel and look at its time series. 
    # For the future, when I do this, I need a file with pixels and row/col (don't mix code and data). 
    # I haven't written this yet. 
    return;

