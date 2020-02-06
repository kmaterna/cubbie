# This is in Python

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import glob, sys, math
import datetime as dt 
from subprocess import call
import sentinel_utilities
import readmytupledata as rmd
import netcdf_read_write as rwr 

def drive_nsbas(swath, intfs, nsbas_min_intfs, sbas_smoothing, wavelength, outdir):
    signal_spread_data=rwr.read_grd("F"+swath+"/"+outdir+"/signalspread.nc");
    intf_tuple = rmd.reader(intfs); 
    velocities = compute_nsbas(intf_tuple, nsbas_min_intfs, sbas_smoothing, wavelength, signal_spread_data); 
    rwr.produce_output_netcdf(x, y, velocities, 'mm/yr', 'F'+swath+'/'+outdir+'/velo_nsbas.grd')
    rwr.produce_output_plot('F'+swath+'/'+outdir+'/velo_nsbas.grd', 'Velocity Profile',
        'F'+swath+'/'+outdir+'/velo_nsbas.png', 'velocity (mm/yr)');
    return;


# ------------ COMPUTE ------------ #
def compute_nsbas(intf_tuple, nsbas_good_perc, smoothing, wavelength, signal_spread_data):
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an NSBAS matrix describing each image that's a real number (not nan). 	
	print("Performing NSBAS on %d files" % (len(intf_tuple.zvalues)) );
	vel = np.zeros([len(intf_tuple.yvalues), len(intf_tuple.xvalues)]);
	c = 0;

	# I can select a pixel and extract its time series. 
	# Should create the function to send this time series out to a file. 
	# Reference Pixel = 748, 321
	# a = do_nsbas_pixel(intf_tuple.zvalues[:,1500,200], intf_tuple.dates_correct, smoothing, wavelength);

	it = np.nditer(intf_tuple.zvalues[0,:,:], flags=['multi_index'], order='F');  # iterate through the 3D array of data
	while not it.finished:
		i=it.multi_index[0];
		j=it.multi_index[1];
		signal_spread = signal_spread_data[i,j];
		pixel_value = intf_tuple.zvalues[:,i,j];
		if signal_spread > nsbas_good_perc: # if we want a calculation for that day... 
			vel[i][j] = do_nsbas_pixel(pixel_value, intf_tuple.dates_correct, smoothing, wavelength); 
		else:
			vel[i,j] = np.nan;
		c=c+1;
		if np.mod(c,10000)==0:
			print('Done with ' + str(c) + ' out of ' + str(len(intf_tuple.xvalues)*len(intf_tuple.yvalues)) + ' pixels')        
		if i==748 and j==321:
			break;
		it.iternext();
	return vel;



def do_nsbas_pixel(pixel_value, date_pairs, smoothing, wavelength):
	# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
	# dates: if we have 35 images, this is the date of each image, in format 
	# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image, in format 2015157_2018177 (real julian day, 1-offset corrected)
	# This solves Gm = d for the movement of the pixel with smoothing. 

	d = np.array([]);
	dates_total=[];

	for i in range(len(date_pairs)):
		dates_total.append(date_pairs[i][0:7])
		dates_total.append(date_pairs[i][8:15])
	dates_total = set(dates_total);
	dates=sorted(dates_total);
	
	date_pairs_used=[];
	for i in range(len(pixel_value)):
		if not math.isnan(pixel_value[i]):
			d = np.append(d, pixel_value[i]);  # removes the nans from the computation. 
			date_pairs_used.append(date_pairs[i]);  # might be a slightly shorter array of which interferograms actually got used. 
	model_num=len(dates)-1;

	G = np.zeros([len(date_pairs_used)+model_num-1, model_num]);  # in one case, 91x35
	# print(np.shape(G));
	
	for i in range(len(d)):  # building G matrix line by line. 
		ith_intf = date_pairs_used[i];
		first_image=ith_intf.split('_')[0]; # in format '2017082'
		second_image=ith_intf.split('_')[1]; # in format '2017094'
		first_index=dates.index(first_image);
		second_index=dates.index(second_image);
		for j in range(second_index-first_index):
			G[i][first_index+j]=1;

	# Building the smoothing matrix with 1, -1 pairs
	for i in range(len(date_pairs_used),len(date_pairs_used)+model_num-1):
		position=i-len(date_pairs_used);
		G[i][position]=1*smoothing;
		G[i][position+1]=-1*smoothing;
		d = np.append(d,0);

	# solving the SBAS linear least squares equation for displacement between each epoch. 
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
	for i in range(len(m)):
		m_cumulative.append(np.sum(m[0:i]));  # The cumulative phase from start to finish! 


	# Solving for linear velocity
	x_axis_datetimes=[dt.datetime.strptime(x,"%Y%j") for x in dates];
	x_axis_days=[(x - x_axis_datetimes[0]).days for x in x_axis_datetimes];  # number of days since first acquisition. 

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

	# plt.figure();
	# plt.plot(x_axis_days[0:-1],m,'b.');
	# plt.plot(x_axis_days,m_cumulative,'g.');
	# plt.plot(x_axis_days, model_line,'--g');
	# plt.xlabel("days");
	# plt.ylabel("cumulative phase");
	# plt.text(0,0,str(vel)+"mm/yr slope");
	# plt.savefig('m_model.eps');

	return vel;



# ------------ OUTPUT ------------ #
def outputs(xdata, ydata, number_of_datas, zdim, vel, out_dir):
	
	netcdf_read_write.produce_output_netcdf(xdata,ydata, vel, 'mm/yr', out_dir+'/vel.grd');
	geocode(out_dir+'/vel.grd',out_dir);

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


def geocode(ifile, directory):
	# geocode: needs vel.grd, vel_ll.grd, vel_ll, and directory 
	stem = ifile.split('/')[-1]  # format: vel.grd
	stem = stem.split('.')[0]   # format: vel
	call(['geocode_mod.csh',stem+'.grd',stem+'_ll.grd',stem+"_ll",directory],shell=False);
	return;


