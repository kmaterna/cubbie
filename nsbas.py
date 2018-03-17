# This is in Python3

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.io import netcdf_file as netcdf
import collections
import glob, sys, math


# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="unwrap_ra";
	
	# Regular parameters
	nsbas_num_toss=6;  # out of 62, how many interferograms can be incoherent? 
	smoothing = 0;  # 
	wavelength = 56;  # mm 

	file_names=glob.glob(file_dir+"/*.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
	return [file_names, nsbas_num_toss, smoothing, wavelength];




# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:
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
	return [xdata, ydata, data_all, dates, date_pairs];

def read_grd(filename):
	data0 = netcdf(filename,'r').variables['z'][::-1];
	data=data0.copy();
	return data;
def read_grd_xy(filename):
	xdata0 = netcdf(filename,'r').variables['x'][::-1];
	ydata0 = netcdf(filename,'r').variables['y'][::-1];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	return [xdata, ydata]; 







# ------------ COMPUTE ------------ #
def compute(xdata, ydata, zdata_all, nsbas_num_toss, dates, date_pairs, smoothing, wavelength):
	[zdim, xdim, ydim] = np.shape(zdata_all)
	print(np.shape(zdata_all[::]))
	# number_of_datas = analyze_coherent_number(zdata_all);
	analyze_velocity_nsbas(zdata_all, nsbas_num_toss, dates, date_pairs, smoothing, wavelength);
	return ;


def analyze_velocity_nsbas(zdata, nsbas_num_toss, dates, date_pairs, smoothing, wavelength):
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an SBAS matrix describing each image that's a real number (not nan). 
	[zdim, xdim, ydim] = np.shape(zdata)
	vel = np.zeros([xdim, ydim]);
	
	for i in range(xdim):  # A loop through each pixel. 
		for j in range(ydim):
			pixel_value = [zdata[k][i][j] for k in range(zdim)];  # slicing the values of phase for a pixel across the various interferograms
			num_reals=0;
			for k in range(zdim):
				if not math.isnan(pixel_value[k]):
					num_reals=num_reals+1;
			
			if num_reals > zdim - nsbas_num_toss:  # If we have a pixel that will be analyzed: Do SBAS

				vel[i][j] = do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength); 
				# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
				# dates: if we have 35 images, this is the date of each image
				# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image
				# This solves Gm = d for the movement of the pixel with smoothing. 

				sys.exit(0);

			else:
				vel[i][j]=np.nan;
				continue;

	return;



def do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength):
	# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
	# dates: if we have 35 images, this is the date of each image
	# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image	
	# for x in range(len(dates)-1):

	return 0;




def analyze_coherent_number(zdata):
	# Analyze the number of coherent acquisitions for each pixel
	[zdim, xdim, ydim] = np.shape(zdata)
	number_of_datas=np.zeros([xdim,ydim]);
	for k in range(zdim):
		for i in range(xdim):
			for j in range(ydim):
				if not math.isnan(zdata[k][i][j]):
					number_of_datas[i][j]=number_of_datas[i][j]+1;
	plt.figure();
	plt.imshow(number_of_datas);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.title("Number of Coherent Intfs (Total = "+str(zdim)+")");
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	plt.colorbar();
	plt.savefig("number_of_coherent_intfs.eps");
	plt.close();
	return number_of_datas;






# ------------ COMPUTE ------------ #
def outputs():
	
	return;








if __name__=="__main__":
	[file_names, nsbas_num_toss, smoothing, wavelength] = configure();
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names);
	compute(xdata, ydata, data_all, nsbas_num_toss, dates, date_pairs, smoothing, wavelength);
	outputs();