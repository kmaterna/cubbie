# This is in Python3

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.io import netcdf_file as netcdf
import collections
import glob, sys, math


# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="Sbas_tests/unwrap_ra";
	file_names=glob.glob(file_dir+"/*.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
	return file_names;


# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:
		data = read_grd(ifile);
		data_all.append(data);
	return [xdata, ydata, data_all];

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
def compute(xdata, ydata, zdata_all):
	[zdim, xdim, ydim] = np.shape(zdata_all)
	print(np.shape(zdata_all[::]))
	analyze_coherent_number(zdata_all);

	return;



def analyze_coherent_number(zdata_all):
	# Analyze the number of coherent acquisitions for each pixel
	[zdim, xdim, ydim] = np.shape(zdata_all)
	number_of_datas=np.zeros([xdim,ydim]);
	for k in range(zdim):
		for i in range(xdim):
			for j in range(ydim):
				if not math.isnan(zdata_all[k][i][j]):
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
	return;

# ------------ COMPUTE ------------ #
def outputs():
	return;




if __name__=="__main__":
	file_names = configure();
	[xdata, ydata, data_all] = inputs(file_names);
	compute(xdata, ydata, data_all);
	outputs();