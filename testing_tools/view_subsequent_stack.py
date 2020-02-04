# Here I'm going to stack lots of unwrapped interferograms together. 
# I will see when the strange grid problem begins. 
# Is it related to a single image? 
# Run this from the F1, F2, and F3 directories. 


import numpy as np
import matplotlib.pyplot as plt
import subprocess
import datetime as dt 
import collections
import sys
import netcdf_read_write as rwr
import glob
import readmytupledata as rmd 
import Super_Simple_Stack as sss


def compute_and_plot(filelist, signal_spread_data):

	wavelength = 56; 

	plt.figure(figsize=(20,20));
	count=1;

	filepathslist = filelist[0:24];
	velocities, x, y = sss.velocity_simple_stack(filepathslist, wavelength, signal_spread_data, 25);  # signal threshold < 100%.  lower signal threshold allows for more data into the stack.  
	plt.imshow(velocities,cmap='jet',vmin=-4, vmax=26,aspect=0.5);
	plt.gca().invert_yaxis();
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	count=count+1;
	plt.savefig("stacking/testing_stack/showing_all_vels.eps");
	rwr.produce_output_netcdf(x, y, velocities, "mm/yr", "stacking/testing_stack/vel.nc");
	return;

def test_for_nans(filename):
	vels = rwr.read_grd(filename);
	plt.figure(figsize=(20,20));
	plt.imshow(np.isnan(vels));
	plt.gca().invert_yaxis();
	plt.savefig("stacking/testing_stack/vel_is_nans.eps");
	
	plt.imshow(vels,cmap='jet',vmin=-4, vmax=26,aspect=0.5);
	plt.gca().invert_yaxis();
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	plt.colorbar()
	plt.savefig("stacking/testing_stack/showing_all_vels.eps");	

	rwr.give_metrics_on_grd(filename);
	return;

if __name__=="__main__":
	filelist = glob.glob("stacking/ref_unwrapped/*_unwrap.grd");
	signal_spread_data=rwr.read_grd("stacking/simple_stack/signalspread.nc");	
	compute_and_plot(filelist, signal_spread_data);

	test_for_nans("stacking/testing_stack/vel.nc");
	# test_for_nans("stacking/simple_stack/velo_simple_stack.grd");
	# test_for_nans("stacking/simple_stack/signalspread.nc");





