# This is in Python

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import glob, sys, math
import datetime as dt 
from subprocess import call
import netcdf_read_write


# ------------- CONFIGURE ------------ # 
def configure(config_params):

	# Time Series parameters
	smoothing = 3;  
	out_dir='nsbas';


	# Setting up the input and output directories. 
	dir_of_interest='referenced_unwrap.grd'
	file_dir="intf_all/"+dir_of_interest;
	file_names=glob.glob(file_dir+"/*_*_unwrap.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern within "+dir_of_interest); sys.exit(1);
	call(['mkdir','-p',out_dir],shell=False);
	return [file_names, config_params.nsbas_min_intfs, smoothing, config_params.wavelength, out_dir];




# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = netcdf_read_write.read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:  # this happens to be in date order on my mac
		data = netcdf_read_write.read_grd(ifile);
		data_all.append(data);
	date_pairs=[];
	dates=[];
	for name in file_names:
		pairname=name.split('/')[-1][0:15];
		date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
		splitname=pairname.split('_');
		dates.append(splitname[0])
		dates.append(splitname[1])
	dates=list(set(dates));
	dates=sorted(dates);
	print(dates);
	return [xdata, ydata, data_all, dates, date_pairs];




# ------------ COMPUTE ------------ #
def compute(xdata, ydata, zdata_all, nsbas_good_num, dates, date_pairs, smoothing, wavelength):
	[zdim, xdim, ydim] = np.shape(zdata_all);
	number_of_datas=np.zeros([xdim, ydim]);
	vel=np.zeros([xdim,ydim]);
	[number_of_datas,zdim] = analyze_coherent_number(zdata_all);
	vel = analyze_velocity_nsbas(zdata_all, number_of_datas, nsbas_good_num, dates, date_pairs, smoothing, wavelength);
	return [vel,number_of_datas,zdim];


def analyze_velocity_nsbas(zdata, number_of_datas, nsbas_good_num, dates, date_pairs, smoothing, wavelength):
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an SBAS matrix describing each image that's a real number (not nan). 
	print("Analyzing the nsbas timeseries per pixel.")
	[zdim, xdim, ydim] = np.shape(zdata)
	vel = np.zeros([xdim, ydim]);
	
	for i in range(xdim):  # A loop through each pixel. 
		for j in range(ydim):
			pixel_value = [zdata[k][i][j] for k in range(zdim)];  # slicing the values of phase for a pixel across the various interferograms
			
			if number_of_datas[i][j] >= nsbas_good_num:  # If we have a pixel that will be analyzed: Do SBAS
				print("%d %d " % (i, j) )

				vel[i][j] = do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength); 
				# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
				# dates: if we have 35 images, this is the date of each image
				# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image
				# This solves Gm = d for the movement of the pixel with smoothing. 

			else:
				vel[i][j]=np.nan;

	return vel;



def do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength):
	# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
	# dates: if we have 35 images, this is the date of each image
	# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image	
	# for x in range(len(dates)-1):

	d = np.array([]);
	dates=sorted(dates);
	date_pairs_used=[];
	for i in range(len(pixel_value)):
		if not math.isnan(pixel_value[i]):
			d = np.append(d, pixel_value[i]);  # removes the nans from the computation. 
			date_pairs_used.append(date_pairs[i]);  # might be a slightly shorter array of which interferograms actually got used. 
	model_num=len(dates)-1;

	G = np.zeros([len(date_pairs_used)+model_num-1, model_num]);  # in one case, 91x35
	print(np.shape(G));
	
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
	# plt.close();

	# if vel>1000:
	# 	vel=1000;
	# if vel<-1000:
	# 	vel=-1000;

	return vel;


def analyze_coherent_number(zdata):
	# Analyze the number of coherent acquisitions for each pixel
	print("Analyzing the number of coherent interferograms per pixel.")
	[zdim, xdim, ydim] = np.shape(zdata)
	number_of_datas=np.zeros([xdim,ydim]);
	for k in range(zdim):
		for i in range(xdim):
			for j in range(ydim):
				if not math.isnan(zdata[k][i][j]):
					number_of_datas[i][j]=number_of_datas[i][j]+1;

	return [number_of_datas, zdim];





# ------------ COMPUTE ------------ #
def outputs(xdata, ydata, number_of_datas, zdim, vel, out_dir):
	
	netcdf_read_write.produce_output_netcdf(xdata, ydata, number_of_datas, 'coherent_intfs', out_dir+'/number_of_datas.grd');
	netcdf_read_write.flip_if_necessary(out_dir+'/number_of_datas.grd');
	netcdf_read_write.produce_output_plot(out_dir+'/number_of_datas.grd', "Number of Coherent Intfs (Total = "+str(zdim)+")", out_dir+'/number_of_coherent_intfs.eps', 'intfs');
	geocode(out_dir+'/number_of_datas.grd',out_dir);

	netcdf_read_write.produce_output_netcdf(xdata,ydata, vel, 'mm/yr', out_dir+'/vel.grd');
	netcdf_read_write.flip_if_necessary(out_dir+'/vel.grd');
	netcdf_read_write.produce_output_plot(out_dir+'/vel.grd','NSBAS LOS Velocity',out_dir+'/vel.eps', 'mm/yr');
	geocode(out_dir+'/vel.grd',out_dir);
	return;




def geocode(ifile, directory):
	# geocode: needs vel.grd, vel_ll.grd, vel_ll, and directory 
	stem = ifile.split('/')[-1]  # format: vel.grd
	stem = stem.split('.')[0]   # format: vel
	call(['geocode_mod.csh',stem+'.grd',stem+'_ll.grd',stem+"_ll",directory],shell=False);
	return;



# ------------ THE MAIN PROGRAM------------ # 

def do_nsbas(config_params):
	[file_names, nsbas_good_num, smoothing, wavelength, out_dir] = configure(config_params);
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names);
	[vel, number_of_datas, zdim] = compute(xdata, ydata, data_all, nsbas_good_num, dates, date_pairs, smoothing, wavelength);
	outputs(xdata, ydata, number_of_datas, zdim, vel, out_dir);
