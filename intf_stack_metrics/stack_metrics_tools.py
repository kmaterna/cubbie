# August 2020
# Here we have a number of tools for post-analysis a stack of interferograms
# Hopefully reducing the number of times I have to write these functions over and over again. 

import numpy as np 
import glob
import matplotlib.pyplot as plt 
import matplotlib
import matplotlib.cm as cm
import netcdf_read_write 


def make_residual_plot(file1, file2, plotname, histname, vmin=-20, vmax=5, 
	title1='', title2='', scalelabel='LOS Velocity', units='mm/yr', flip_sign=False):
	"""
	A basic function that takes two co-registered grids and subtracts them, showing residuals in the third panel
	and histogram of residuals in separate plot. 
	"""
	data1 = netcdf_read_write.read_grd(file1);
	data2 = netcdf_read_write.read_grd(file2);
	if flip_sign:
		data1 = -1 * data1;
		data2 = -1 * data2;
	residuals = np.subtract(data1, data2);
	residuals_vector = np.reshape(residuals, (np.shape(residuals)[0]*np.shape(residuals)[1],));
	

	fig,axarr = plt.subplots(1,3,sharey=True, figsize=(20, 8), dpi=300);
	axarr[0].imshow(data1, vmin=vmin, vmax=vmax, cmap='rainbow');
	axarr[0].tick_params(labelsize=16);
	axarr[0].set_title(title1,fontsize=20);
	axarr[0].invert_yaxis()

	axarr[1].imshow(data2, vmin=vmin, vmax=vmax, cmap='rainbow');
	axarr[1].tick_params(labelsize=16);
	axarr[1].set_title(title2,fontsize=20);
	axarr[1].invert_yaxis()	

	axarr[2].imshow(residuals, vmin=-10, vmax=10, cmap='rainbow');
	axarr[2].tick_params(labelsize=16);
	axarr[2].set_title('Residuals',fontsize=20);
	axarr[2].invert_yaxis()	

	# Fancy color bar #1
	cbarax = fig.add_axes([0.85, 0.08, 0.1, 0.9],visible=False);
	color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
	custom_cmap.set_array(np.arange(vmin, vmax, 0.1));
	cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
	cb.set_label(scalelabel+' ('+units+')', fontsize=18);
	cb.ax.tick_params(labelsize=16);
	
	# Fancy color bar for residuals
	cbarax = fig.add_axes([0.68, 0.05, 0.1, 0.9],visible=False);
	color_boundary_object = matplotlib.colors.Normalize(vmin=-10, vmax=10);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
	custom_cmap.set_array(np.arange(vmin, vmax, 0.1));
	cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='horizontal');
	cb.set_label('Residual ('+units+')', fontsize=16);
	cb.ax.tick_params(labelsize=14);	
	plt.savefig(plotname);
	plt.close();

	plt.figure(dpi=300, figsize=(8,6))
	plt.hist(residuals_vector, bins=50, color='orange');
	plt.yscale('log');
	plt.ylabel('Number of Pixels',fontsize=20)
	plt.xlabel('Residuals ('+units+')', fontsize=20)
	plt.gca().tick_params(axis='both',which='major',labelsize=16);
	plt.savefig(histname);
	plt.close();
	return;


def plot_two_nonregistered_grids(file1, file2, plotname, vmin=-20, vmax=5, flip_sign1=False, flip_sign2=False,
	title1='', title2='', scalelabel='Velocity (mm/yr)'):
	"""
	A little function that plots two grid files in subplots side by side
	(they don't have to have the same registration, so no need to compute residuals)
	"""
	data1 = netcdf_read_write.read_grd(file1);
	data2 = netcdf_read_write.read_grd(file2);
	if flip_sign1:
		data1 = -1*data1;
	if flip_sign2:
		data2 = -1*data2;
	
	fig,axarr = plt.subplots(1,2,sharey=False, figsize=(15, 8), dpi=300);
	axarr[0].imshow(data1, vmin=vmin, vmax=vmax, cmap='rainbow');
	axarr[0].tick_params(labelsize=16);
	axarr[0].set_title(title1,fontsize=20);
	axarr[0].invert_yaxis()

	axarr[1].imshow(data2, vmin=vmin, vmax=vmax, cmap='rainbow');
	axarr[1].tick_params(labelsize=16);
	axarr[1].set_title(title2,fontsize=20);
	axarr[1].invert_yaxis()	

	cbarax = fig.add_axes([0.85, 0.08, 0.1, 0.9],visible=False);
	color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
	custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
	custom_cmap.set_array(np.arange(vmin, vmax, 0.1));
	cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
	cb.set_label(scalelabel, fontsize=18);
	cb.ax.tick_params(labelsize=16);
	plt.savefig(plotname);
	plt.close();	
	return;


def histogram_of_grd_file_values(filename, varname=Deviation, plotname=histogram_values.png):
	"""
	simple plot to make a histogram of a grid file
	"""
	z = netcdf_read_write.read_grd(filename);
	z_vector = np.reshape(z, (np.shape(z)[0]*np.shape(z)[1],))
	plt.figure(dpi=250, figsize=(8,7));
	plt.hist(z_vector, bins=50, color='orange');
	plt.yscale('log');
	plt.ylabel('Number of Pixels',fontsize=20)
	plt.xlabel(varname, fontsize=20)
	plt.gca().tick_params(axis='both',which='major',labelsize=16);
	plt.savefig(plotname);
	plt.close();
	return;




