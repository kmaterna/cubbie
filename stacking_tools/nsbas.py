# This is in Python

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import glob, sys, math
import datetime as dt 
from subprocess import call
import sentinel_utilities
import stacking_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 
import stack_corr

def drive_velocity_nsbas(swath, intfs, nsbas_min_intfs, sbas_smoothing, wavelength, outdir):
    signal_spread_file=outdir+"/signalspread_cut.nc"
    intf_tuple = rmd.reader(intfs); 
    make_stack_corr_custom(intf_tuple, signal_spread_file);  # for safety, let's make signalspread again. 
    signal_spread_data=rwr.read_grd(signal_spread_file);
    velocities = driver_vels(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data); 
    
    # An issue with the lightswitch at Moffett turning off... 
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr', '/Users/kmaterna/Documents/testvelo/F'+swath+'_'+'velo_nsbas.grd');  # just in case we screw up
    rwr.produce_output_plot('/Users/kmaterna/Documents/testvelo/F'+swath+'_velo_nsbas.grd', 'LOS Velocity',
        '/Users/kmaterna/Documents/testvelo/F'+swath+'_velo_nsbas.png', 'velocity (mm/yr)');

    # The more general case. 
    rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, velocities, 'mm/yr', outdir+'/velo_nsbas.grd');
    rwr.produce_output_plot(outdir+'/velo_nsbas.grd', 'LOS Velocity', outdir+'/velo_nsbas.png', 'velocity (mm/yr)');
    return;

def drive_ts_nsbas(config_params):
	lons, lats, names, swaths, rows, cols = stacking_utilities.get_rows_cols(config_params.ts_points_file);
	if len(rows)==0:
		return;
	drive_ts_nsbas_one_swath(config_params, '1', rows, cols, swaths, names, lons, lats, config_params.sbas_smoothing, config_params.wavelength);
	drive_ts_nsbas_one_swath(config_params, '2', rows, cols, swaths, names, lons, lats, config_params.sbas_smoothing, config_params.wavelength);
	drive_ts_nsbas_one_swath(config_params, '3', rows, cols, swaths, names, lons, lats, config_params.sbas_smoothing, config_params.wavelength);
	return;

# FOR A GIVEN SWATH, LET'S GET SOME PIXELS AND OUTPUT THEIR TS. 
def drive_ts_nsbas_one_swath(config_params, select_swath, rows, cols, swaths, names, lons, lats, smoothing, wavelength):
	outdir = config_params.ts_output_dir+"/ts";
	print("TS OUTPUT DIR IS: " + outdir)
	call(['mkdir','-p',outdir],shell=False);
	
	rows=np.array(rows);
	cols=np.array(cols);
	names=np.array(names);
	lons=np.array(lons);
	lats=np.array(lats);
	select_rows = rows[np.array(swaths)==select_swath];
	select_cols = cols[np.array(swaths)==select_swath];
	select_names = names[np.array(swaths)==select_swath];
	select_lons = lons[np.array(swaths)==select_swath];
	select_lats = lats[np.array(swaths)==select_swath];
	print("For Swath %s, extracting time series for:" % (select_swath) );
	for i in range(len(select_rows)):
		print(select_rows[i], select_cols[i]);
	if len(select_rows)==0:
		return;	
	intfs = stacking_utilities.make_selection_of_intfs(config_params, swath=select_swath);
	intf_tuple = rmd.reader(intfs);
	for i in range(len(select_rows)):
		pixel_value = intf_tuple.zvalues[:,select_rows[i],select_cols[i]];
		vel, dts, m_cumulative = do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, smoothing, wavelength, full_ts_return=True); 
		m_cumulative=[i*-1 for i in m_cumulative];  # My sign convention seems to be opposite to Katia's
		nsbas_ts_outputs(dts, m_cumulative, select_swath, select_rows[i], select_cols[i], select_names[i], select_lons[i], select_lats[i], outdir);
	return;

# This might be redundant, not sure. 
def single_pixel_ts(intf_tuple, pixel_coords, sbas_smoothing, wavelength, signal_spread_file,outdir, coh_tuple=[]):
	# Plot a single pixel's time series. 
	datestrs, x_dts, xdates = get_TS_dates(intf_tuple.dates_correct);
	pixel_value = intf_tuple.zvalues[:,pixel_coords[0],pixel_coords[1]];
	if coh_tuple==[]:
		coh_value=[];
	else:
		coh_value = coh_tuple.zvalues[:,pixel_coords[0], pixel_coords[1]];
	ts = do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, sbas_smoothing, wavelength, datestrs, xdates, coh_value, 
		full_ts_return=True);
	output_single_pixel_ts(pixel_coords, xdates, ts, signal_spread_file, outdir);
	return;


# Before velocities, might want to make stack_corr_custom if you have excluded significant data. 
def make_stack_corr_custom(mydata,signal_spread_file):
    # Stack corr for this exact calculation (given right inputs and outputs). 
    a=stack_corr.stack_corr(mydata, np.nan);
    rwr.produce_output_netcdf(mydata.xvalues, mydata.yvalues, a, 'Percentage', signal_spread_file);
    rwr.produce_output_plot(signal_spread_file, 'Signal Spread', signal_spread_file+'.png', 
        'Percentage of coherence', aspect=1/4, invert_yaxis=False )
    return;



# ------------ COMPUTE ------------ #
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an NSBAS matrix describing each image that's a real number (not nan). 

def driver_vels(intf_tuple, nsbas_good_perc, smoothing, wavelegnth, signal_spread_data, coh_tuple = []):
	retval = np.zeros([len(intf_tuple.yvalues), len(intf_tuple.xvalues)]);
	def packager_function(i, j, intf_tuple):
		# Giving access to all these variables
		return compute_vel(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data, datestrs, x_axis_days, coh_tuple);
	retval = iterator_func(intf_tuple, nsbas_good_perc, smoothing, wavelength, 
		signal_spread_data, packager_function, coh_tuple);  # how does it get compute_vels()?
	return retval

def driver_Full_TS(intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data, coh_tuple = []):
	datestrs, x_dts, x_axis_days = get_TS_dates(intf_tuple.dates_correct); 
	# Establishing the return array
	empty_vector=[np.empty(np.shape(x_axis_days))];
	retval = [ [ empty_vector for i in range(len(intf_tuple.xvalues))] for j in range(len(intf_tuple.yvalues))];
	def packager_function(i, j, intf_tuple):
		# Giving access to all these variables. 
		return compute_TS(i,j,intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data, datestrs, x_axis_days, coh_tuple);
	retval = iterator_func(intf_tuple, packager_function, retval);  
	return retval;


def iterator_func(intf_tuple, func, retval):
	# This iterator performs a for loop. It assumes the return value can be stored in an array of ixj
	# if np.shape(retval) != np.shape(signal_spread_data):
	# 	print("ERROR: signal spread does not match input data. Stopping immediately. ");
	# 	print("Shape of signal spread:", np.shape(signal_spread_data));
	# 	print("Shape of data array:", np.shape(intf_tuple.zvalues[0]));
	# 	sys.exit(0);
	print("Performing NSBAS on %d files" % (len(intf_tuple.zvalues)) );
	print("Started at: ");
	print(dt.datetime.now());
	c = 0;
	it = np.nditer(intf_tuple.zvalues[0,:,:], flags=['multi_index'], order='F');  # iterate through the 3D array of data
	while not it.finished:
		i=it.multi_index[0];
		j=it.multi_index[1];
		retval[i][j]=func(i,j,intf_tuple);
		c=c+1;
		if np.mod(c,10000)==0:
			print('Done with ' + str(c) + ' out of ' + str(len(intf_tuple.xvalues)*len(intf_tuple.yvalues)) + ' pixels')        
		# if c==40000:
		# 	break;
		it.iternext();
	print("Finished at: ");
	print(dt.datetime.now());
	return retval;


def compute_vel(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data, datestrs, x_axis_days, coh_tuple=[]):
	signal_spread = signal_spread_data[i,j];
	pixel_value = intf_tuple.zvalues[:,i,j];
	if coh_tuple==[]:
		coh_value=[];
	else:
		coh_value = coh_tuple.zvalues[:,i,j];
	if signal_spread > nsbas_good_perc: 
		vel = do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, smoothing, wavelength, datestrs, x_axis_days, coh_value, full_ts_return=False); 
	else:
		vel = np.nan;
	return vel;

def compute_TS(i, j, intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data, datestrs, x_axis_days, coh_tuple=[]):
	# For a given iteration, what's the right time series? 
	empty_vector = np.empty(np.shape(x_axis_days));  # Length of the TS model
	empty_vector[:] = np.nan;
	signal_spread = signal_spread_data[i,j];
	pixel_value = intf_tuple.zvalues[:,i,j];
	if coh_tuple==[]:
		coh_value=[];
	else:
		coh_value = coh_tuple.zvalues[:,i,j];
	if signal_spread > nsbas_good_perc:
		vector = do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, smoothing, wavelength, datestrs, x_axis_days, coh_value, full_ts_return=True); 
		TS = [vector];
	else:
		TS = [empty_vector];
	return TS;


def get_TS_dates(date_pairs):
	# Get me the x axis associated with a certain set of interferograms
	# Returns a list of n strings for each satellite pass
	# Also the dt objects
	# Also the days since the start of the track. 
	dates_total=[];
	for i in range(len(date_pairs)):
		dates_total.append(date_pairs[i][0:7])
		dates_total.append(date_pairs[i][8:15])
	dates_total = set(dates_total);
	datestrs=sorted(dates_total);
	x_axis_datetimes=[dt.datetime.strptime(x,"%Y%j") for x in datestrs];
	x_axis_days=[(x - x_axis_datetimes[0]).days for x in x_axis_datetimes];  # number of days since first acquisition. 
	return datestrs, x_axis_datetimes, x_axis_days;


def do_nsbas_pixel(pixel_value, date_pairs, smoothing, wavelength, datestrs, x_axis_days, coh_value=[], full_ts_return=False):
	# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
	# dates: if we have 35 images, this is the date of each image, in format 
	# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image, in format 2015157_2018177 (real julian day, 1-offset corrected)
	# This solves Gm = d for the movement of the pixel with smoothing. 
	# If coh_value is an array, we do weighted least squares. 
	# Need to implement WLS next. 

	d = np.array([]);
	diagonals = [];
	date_pairs_used=[];

	for i in range(len(pixel_value)):
		if not math.isnan(pixel_value[i]):
			d = np.append(d, pixel_value[i]);  # removes the nans from the computation. 
			date_pairs_used.append(date_pairs[i]);  # might be a slightly shorter array of which interferograms actually got used. 
			if coh_value != []:
				diagonals.append(np.power(coh_value[i],2));  # using coherence squared as the weighting. 
			else:
				diagonals.append(1);
	model_num=len(datestrs)-1;

	# Average weights at the end of the weigting vector
	Wavg = np.nanmean(diagonals)
	for i in range(model_num-1):
		diagonals.append(Wavg)
	W = np.diag(diagonals);

	G = np.zeros([len(date_pairs_used)+model_num-1, model_num]);  # in one case, 91x35
	# print(np.shape(G));
	
	# building G matrix line by line. 
	for i in range(len(d)):  
		ith_intf = date_pairs_used[i];
		first_image=ith_intf.split('_')[0]; # in format '2017082'
		second_image=ith_intf.split('_')[1]; # in format '2017094'
		first_index=datestrs.index(first_image);
		second_index=datestrs.index(second_image);
		for j in range(second_index-first_index):
			G[i][first_index+j]=1;

	# Building the smoothing matrix with 1, -1 pairs
	for i in range(len(date_pairs_used),len(date_pairs_used)+model_num-1):
		position=i-len(date_pairs_used);
		G[i][position]=1*smoothing;
		G[i][position+1]=-1*smoothing;
		d = np.append(d,0);
  
	# solving the SBAS linear least squares equation for displacement between each epoch.
	if coh_value != []:
		GTWG = np.dot(np.transpose(G), np.dot(W,G))
		GTWd = np.dot(np.transpose(G), np.dot(W,d))
		m = np.dot(np.linalg.inv(GTWG), GTWd )
	else:
		m = np.linalg.lstsq(G,d)[0];

	# modeled_data=np.dot(G,m);
	# plt.figure();
	# plt.plot(d,'.b');
	# plt.plot(modeled_data,'.--g');
	# plt.savefig('d_vs_m.eps')
	# plt.close();
	# sys.exit(0);

	# Adding up all the displacement. 
	m_cumulative=[];
	m_cumulative.append(0);
	for i in range(1,len(m)+1):
		m_cumulative.append(np.sum(m[0:i]));  # The cumulative phase from start to finish! 

	# Solving for linear velocity
	x=np.zeros([len(x_axis_days),2]);
	y=np.array([]);
	for i in range(len(x_axis_days)):
		x[i][0]=x_axis_days[i];
		x[i][1]=1;  
		y=np.append(y,[m_cumulative[i]]);
	model_slopes = np.linalg.lstsq(x,y)[0];  # units: phase per day. 
	model_line = [model_slopes[1]+ x*model_slopes[0] for x in x_axis_days];

	# Velocity conversion: units in mm / year
	vel=model_slopes[0];  # in radians per day
	vel=vel*wavelength*365.24/2.0/(2*np.pi);

	disp_ts = [i*wavelength/(4*np.pi) for i in m_cumulative];
	# plt.figure();
	# plt.plot(x_axis_days[0:-1],m,'b.');
	# plt.plot(x_axis_days,m_cumulative,'g.');
	# plt.plot(x_axis_days, model_line,'--g');
	# plt.xlabel("days");
	# plt.ylabel("cumulative phase");
	# plt.text(0,0,str(vel)+"mm/yr slope");
	# plt.savefig('m_model.eps');

	if full_ts_return:
		return disp_ts;
	else:
		return vel;



# ------------ OUTPUT ------------ #




def output_single_pixel_ts(pixel_coords, xdates, ts, signal_spread_file, outdir):
    print("Random Pixel");
    fig = plt.figure()
    plt.plot(xdates, ts,'.-',markersize=10);
    plt.title("Pixel "+str(pixel_coords[0])+", "+str(pixel_coords[1]));
    plt.ylabel('Displacement (mm)');
    plt.savefig(outdir+'/test_pixel.png');
    rwr.produce_output_plot(signal_spread_file, 'Signal Spread', outdir+'/show_test_pixels.png', 
        'Percentage of coherence', aspect=1/4, invert_yaxis=False, dot_points = [[pixel_coords[1]],[pixel_coords[0]]] );
    return;



def nsbas_ts_outputs(dts, m_cumulative, swath, row, col, name, lon, lat, outdir):
	# This outdir is expected to be something like "F2/stacking/nsbas/ts". 
	# It must exist before the function is called. 

	mean_disp = np.nanmean(m_cumulative);
	plotting_ts = [i-mean_disp for i in m_cumulative];

	plt.figure();
	plt.plot(dts,plotting_ts,'b.');
	plt.xlabel("Time");
	plt.ylabel("Displacement (mm)");
	plt.title(str(swath)+' '+str(row)+' '+str(col)+' '+str(lon)+' '+str(lat)+' '+str(name));
	plt.ylim([-40,50]);
	plt.savefig(outdir+'/'+str(name)+'_'+str(lon)+'_'+str(lat)+'_disp.eps');

	ofile=open(outdir+'/'+str(name)+'_'+str(row)+'_'+str(col)+'_record.txt','w');
	for i in range(len(dts)):
		ofile.write("%s %f %f %s %d %d " % (name, lon, lat, swath, row, col) );
		ofile.write(dt.datetime.strftime(dts[i],"%Y-%m-%d"));
		ofile.write(" %f\n" % (m_cumulative[i]) );
	ofile.close();
	return;

def outputs(xdata, ydata, number_of_datas, zdim, vel, out_dir):

	# Visualizing the velocity field in a few different ways. 
	zdata2=np.reshape(vel, [len(xdata)*len(ydata), 1])
	zdata2=sentinel_utilities.remove_nans_array(zdata2);
	plt.figure();
	plt.hist(zdata2,bins=80);
	plt.gca().set_yscale('log');
	plt.title('Pixels by Velocity: mean=%.2fmm/yr, sdev=%.2fmm/yr' % (np.mean(zdata2), np.std(zdata2)) )
	plt.ylabel('Number of Pixels');
	plt.xlabel('LOS velocity (mm/yr)')
	plt.grid('on');
	plt.savefig(out_dir+'/velocity_hist_log.png');
	plt.close();

	plt.figure();
	plt.gca().set_yscale('linear');
	plt.title('Pixels by Velocity: mean=%.2fmm/yr, sdev=%.2fmm/yr' % (np.mean(zdata2), np.std(zdata2)) )
	plt.hist(zdata2,bins=80);
	plt.ylabel('Number of Pixels');
	plt.xlabel('LOS velocity (mm/yr)')
	plt.grid('on');
	plt.savefig(out_dir+'/velocity_hist_lin.png');
	plt.close();


	plt.figure(figsize=(8,10));
	plt.imshow(vel,aspect=0.5,cmap='jet',vmin=-30, vmax=30);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	plt.title("Velocity");
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label("mm/yr", size=16);
	plt.savefig(out_dir+"/vel_cutoff.png");
	plt.close();

	plt.figure(figsize=(8,10));
	plt.imshow(vel,aspect=0.5,cmap='jet',vmin=-150, vmax=150);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	plt.title("Velocity");
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label("mm/yr", size=16);
	plt.savefig(out_dir+"/vel.png");
	plt.close();

	return;
