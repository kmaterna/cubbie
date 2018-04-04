# This is in Python

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.io.netcdf as netcdf
import collections
import glob, sys, math
import datetime as dt 


# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="unwrap_ra";
	
	# Regular parameters
	nsbas_num_toss=15;  # out of 80, how many interferograms can be incoherent? 
	smoothing = 3;  # 
	wavelength = 56;  # mm 

	file_names=glob.glob(file_dir+"/*.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
	return [file_names, nsbas_num_toss, smoothing, wavelength];




# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:  # this happens to be in date order on my mac
		data = read_grd(ifile);
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

def read_grd(filename):
	data0 = netcdf.netcdf_file(filename,'r').variables['z'][::-1];
	data=data0.copy();
	return data;
def read_grd_xy(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['x'][::-1];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['y'][::-1];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	return [xdata, ydata]; 





# ------------ COMPUTE ------------ #
def compute(xdata, ydata, zdata_all, nsbas_num_toss, dates, date_pairs, smoothing, wavelength):
	[zdim, xdim, ydim] = np.shape(zdata_all)
	vel=np.zeros([xdim,ydim]);
	[number_of_datas,zdim] = analyze_coherent_number(zdata_all);
	vel = analyze_velocity_nsbas(zdata_all, number_of_datas, nsbas_num_toss, dates, date_pairs, smoothing, wavelength);
	return [vel,number_of_datas,zdim];


def analyze_velocity_nsbas(zdata, number_of_datas, nsbas_num_toss, dates, date_pairs, smoothing, wavelength):
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an SBAS matrix describing each image that's a real number (not nan). 
	[zdim, xdim, ydim] = np.shape(zdata)
	vel = np.zeros([xdim, ydim]);
	
	for i in range(xdim):  # A loop through each pixel. 
		for j in range(ydim):
			pixel_value = [zdata[k][i][j] for k in range(zdim)];  # slicing the values of phase for a pixel across the various interferograms
			# num_reals=0;
			# for k in range(zdim):
			# 	if not math.isnan(pixel_value[k]):
			# 		num_reals=num_reals+1;
			
			if number_of_datas[i][j] > zdim - nsbas_num_toss:  # If we have a pixel that will be analyzed: Do SBAS

				vel[i][j] = do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength); 
				# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
				# dates: if we have 35 images, this is the date of each image
				# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image
				# This solves Gm = d for the movement of the pixel with smoothing. 


			else:
				vel[i][j]=np.nan;
				continue;

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
	[zdim, xdim, ydim] = np.shape(zdata)
	number_of_datas=np.zeros([xdim,ydim]);
	for k in range(zdim):
		for i in range(xdim):
			for j in range(ydim):
				if not math.isnan(zdata[k][i][j]):
					number_of_datas[i][j]=number_of_datas[i][j]+1;

	return [number_of_datas, zdim];





# ------------ COMPUTE ------------ #
def outputs(xdata, ydata, number_of_datas, zdim, vel):
	
	plot_title="Number of Coherent Intfs (Total = "+str(zdim)+")"
	produce_output_netcdf(xdata, ydata, number_of_datas, 'coherent_intfs', 'number_of_datas.nc');
	produce_output_plot('number_of_datas.nc', plot_title, 'number_of_coherent_intfs.eps', 'intfs');
	produce_output_netcdf(xdata,ydata, vel, 'mm/yr', 'vel.nc', 'NSBAS LOS Velocity', 'vel.eps');
	produce_output_plot('vel.nc','NSBAS LOS Velocity','vel.eps', 'mm/yr');
	return;


def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname):
	# # Write the netcdf velocity grid file.  This works, but doesn't get read from gmt grdinfo. Not sure why. 
	f=netcdf.netcdf_file(netcdfname,'w');
	f.history = 'Created for a test';
	f.createDimension('x',len(xdata));
	f.createDimension('y',len(ydata));
	print(np.shape(vel));
	x=f.createVariable('x',float,('x',))
	x[:]=xdata;
	x.units = 'range';
	y=f.createVariable('y',float,('y',))
	y[:]=ydata;
	y.units = 'azimuth';
	z=f.createVariable('z',float,('y','x',));
	z[:,:]=zdata;
	z.units = zunits;
	f.close();
	return;



def produce_output_plot(netcdfname, plottitle, filename, clabel):

	# Read in the dataset you just wrote. 
	fr = netcdf.netcdf_file(netcdfname,'r');
	xread=fr.variables['x'];
	yread=fr.variables['y'];
	zread=fr.variables['z'];
	zread_copy=zread[:][:].copy();

	# Make a plot
	fig = plt.figure(figsize=(12,10));
	ax1 = fig.add_axes([-0.5, 0.1, 1.5, 0.8]);
	plt.imshow(zread_copy,aspect=0.3);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	plt.title(plottitle);
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label(cblabel, size=16);
	plt.savefig(filename);
	plt.close();
	return;




if __name__=="__main__":
	[file_names, nsbas_num_toss, smoothing, wavelength] = configure();
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names);
	[vel, number_of_datas, zdim] = compute(xdata, ydata, data_all, nsbas_num_toss, dates, date_pairs, smoothing, wavelength);
	outputs(xdata, ydata, number_of_datas, zdim, vel);

