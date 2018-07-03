# Atmospheric Phase Screen method of Tymofyeyeva and Fialko, 2015


import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call 
import glob, sys
import netcdf_read_write


def main_function(staging_directory, out_dir):
	[file_names] = configure(staging_directory, out_dir);
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names);
	compute(xdata, ydata, data_all, dates, date_pairs);
	return;

# ------------- CONFIGURE ------------ # 
def configure(staging_directory, out_dir):
	# Setting up the input and output directories. 
	file_names=glob.glob(staging_directory+"/*_*_unwrap.grd");
	file_names=sorted(file_names);
	print("Performing APS on files in %s " % staging_directory);
	if len(file_names)==0:
		print("Error! No files matching search pattern within "+staging_directory); sys.exit(1);
	call(['mkdir','-p',out_dir],shell=False);
	return [file_names];


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

# ----------- COMPUTE ------------- #
def compute(xdata, ydata, data_all, dates, date_pairs):
	
	# Goal: construct an atmoshperic phase mask for each image. 
	# Try one example first. Day 2017298, or dates[44]. 

	print(dates[44]);
	print(date_pairs[62]);
	print(date_pairs[63]);
	print(date_pairs[65]);
	print(date_pairs[66]);

	zdim, rowdim, coldim = np.shape(data_all);

	intf0=62;
	intf1=63;
	intf2=65;
	intf3=66;

	vmin=-20;
	vmax=20;

	my_aps=np.zeros((rowdim,coldim));
	intf1_corrected=np.zeros((rowdim,coldim));
	intf2_corrected=np.zeros((rowdim,coldim));
	loopcheck=np.zeros((rowdim,coldim));
	for i in range(rowdim):
		for j in range(coldim):

			#my_aps[i][j]=(1.0/2) * (data_all[intf1][i][j] - data_all[intf2][i][j]);
			inside_loop = 1;
			outside_loop = 1;
			if np.isnan(data_all[intf1][i][j]) or np.isnan(data_all[intf2][i][j]):
				inside_loop = 0;
			if np.isnan(data_all[intf0][i][j]) or np.isnan(data_all[intf3][i][j]):
				outside_loop = 0;
			loopcheck[i][j]=inside_loop+outside_loop;

			if inside_loop == 0 and outside_loop == 0:
				my_aps[i][j]=np.nan
			elif inside_loop == 1 and outside_loop == 0:
				my_aps[i][j] = (1.0/2) * (data_all[intf1][i][j] - data_all[intf2][i][j]);
			elif inside_loop == 0 and outside_loop == 1:
				my_aps[i][j] = (1.0/2) * (data_all[intf0][i][j] - data_all[intf3][i][j]);
			elif inside_loop == 1 and outside_loop ==1 :
				my_aps[i][j] = (1.0/4) * (data_all[intf1][i][j] - data_all[intf2][i][j] + data_all[intf0][i][j] - data_all[intf3][i][j]);

			intf1_corrected[i][j]=data_all[intf1][i][j] - my_aps[i][j];
			intf2_corrected[i][j]=data_all[intf2][i][j] + my_aps[i][j];




	plt.figure();
	plt.subplot(2,3,1);
	plt.imshow(data_all[intf1],cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title(date_pairs[intf1]);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,2)
	plt.imshow(data_all[intf2],cmap='hsv', vmin=vmin, vmax=vmax);
	plt.title(date_pairs[intf2]);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,3);
	plt.imshow(my_aps, cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title('APS for ' +dates[44])
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,4);
	plt.imshow(intf1_corrected,cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title('intf1_corr');
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,5)
	plt.imshow(intf2_corrected,cmap='hsv', vmin=vmin, vmax=vmax);
	plt.title('intf2_corr');
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,6)
	plt.imshow(loopcheck,cmap='hsv',vmin=0, vmax=2);
	plt.title('loopcheck');
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])
	plt.savefig('test_image2.png');
	plt.close();
	return;


# ----------- OUTPUTS ------------- # 
def outputs():
	return;