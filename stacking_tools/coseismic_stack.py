# Code to take a set of interferograms that span a particular earthquake and generate an average. 
# The average should contain less noise than the original interferograms.

import numpy as np
import sys
import readmytupledata as rmd
import netcdf_read_write as rwr 


def drive_coseismic_stack_gmtsar(intf_files, wavelength, outdir):
	intf_tuple = rmd.reader(intfs); 
	average_coseismic = get_avg_coseismic(intf_tuple, wavelength);
	rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, 'mm', outdir+'/coseismic.grd');
	rwr.produce_output_plot(outdir+'/coseismic.grd', 'LOS Displacement', outdir+'/coseismic.png', 'displacement (mm)');
	return;

def drive_coseismic_stack_isce(intf_files, wavelength, outdir):
	intf_tuple = rmd.reader_isce(intf_files); 
	average_coseismic = get_avg_coseismic(intf_tuple, wavelength);
	rwr.produce_output_netcdf(intf_tuple.xvalues, intf_tuple.yvalues, average_coseismic, 'mm', outdir+'/coseismic.grd');
	rwr.produce_output_plot(outdir+'/coseismic.grd', 'LOS Displacement', outdir+'/coseismic.png', 'displacement (mm)', 
		aspect=1/5, invert_yaxis=False);
	return;

def get_avg_coseismic(intf_tuple, wavelength):
	# I could send this into the iterator_func in NSBAS if I wanted to. 
	disp = np.zeros(np.shape(intf_tuple.zvalues[0,:,:]));
	for i in range(len(intf_tuple.yvalues)):
		for j in range(len(intf_tuple.xvalues)):
			disp[i][j]=np.nanmean(intf_tuple.zvalues[:,i,j]) * wavelength/(4*np.pi);
	return disp;