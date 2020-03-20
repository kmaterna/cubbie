from subprocess import call
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import datetime as dt 
import sys
import sentinel_utilities
import stacking_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 
import nsbas
import stack_corr



def drive_velocity_nsbas(swath, intfs, nsbas_min_intfs, sbas_smoothing, wavelength, outdir):

    intf_tuple = rmd.reader_isce(intfs); 
    make_stack_corr_custom(intf_tuple, swath, outdir);  # for safety, let's make signalspread again. 
    signal_spread_data=rwr.read_grd('F'+swath+'/'+outdir+"/signalspread_cut.nc");

    # pixel_value = intf_tuple.zvalues[:,300,20];
    # vel, xdates, ts = nsbas.do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, sbas_smoothing, wavelength, full_ts_return=True)
    # print("Random Pixel");
    # fig = plt.figure()
    # plt.plot(xdates, ts,'.-',markersize=10);
    # plt.savefig('test_pixel.png');

    # TIME SERIES
    TS, xdates = nsbas.compute_fullTS_nsbas(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data);
    TS_file = 'F'+swath+'/'+outdir+"/TS.nc";
    TS_image_file = 'F'+swath+'/'+outdir+"/TS.png";
    rwr.produce_output_timeseries(intf_tuple.xvalues, intf_tuple.yvalues, TS, xdates, 'mm', TS_file);
    plot_full_timeseries(TS_file, xdates, TS_image_file);


    # AVERAGE VELOCITIES
    # velocities = nsbas.compute_nsbas(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data); 
    # rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr', 'F'+swath+'/'+outdir+'/velo_nsbas.grd');
    # rwr.produce_output_plot('F'+swath+'/'+outdir+'/velo_nsbas.grd', 'LOS Velocity',
    #     'F'+swath+'/'+outdir+'/velo_nsbas.png', 'velocity (mm/yr)',aspect=1/4, invert_yaxis=False);

    return;

def make_stack_corr_custom(mydata, swath, outdir):
    # Stack corr for this exact calculation (given right inputs and outputs). 
    a=stack_corr.stack_corr(mydata, np.nan);
    rwr.produce_output_netcdf(mydata.xvalues, mydata.yvalues, a, 'Percentage', 'F'+swath+'/'+outdir+'/signalspread_cut.nc');
    rwr.produce_output_plot('F'+swath+'/'+outdir+'/signalspread_cut.nc', 'Signal Spread', 'F'+swath+'/'+outdir+'/signalspread_cut.png', 
        'Percentage of coherence', aspect=1/4, invert_yaxis=False )
    return;


def get_axarr_numbers(rows, cols, idx):
    # Given an incrementally counting idx number and a subplot dimension, where is our plot? 
    total_plots = rows*cols;
    col_num = np.mod(idx, cols);
    row_num = int(np.floor(idx/cols));
    return row_num, col_num;


def plot_full_timeseries(TS_file, xdates, TS_image_file):
    # Make a nice time series plot. 
    tdata, xdata, ydata, TS_array = rwr.read_3D_netcdf(TS_file);
    vmin=-50;
    vmax=200;

    f, axarr = plt.subplots(3,4,figsize=(16,10));
    for i in range(len(xdates)):
        rownum, colnum = get_axarr_numbers(3,4,i);
        axarr[rownum][colnum].imshow(TS_array[i,:,:],aspect=1/4,cmap='rainbow',vmin=vmin,vmax=vmax);
        titlestr = dt.datetime.strftime(xdates[i],"%Y-%m-%d");
        axarr[rownum][colnum].get_xaxis().set_visible(False);
        axarr[rownum][colnum].set_title(titlestr,fontsize=20);

    cbarax = f.add_axes([0.75,0.35,0.2,0.3],visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(TS_image_file);
    return;