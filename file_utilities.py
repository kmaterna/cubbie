# Utilities to convert file types

import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call
import glob
import isce_read_write
import netcdf_read_write


def gmtsar_nc_stack_2_isce_stack(ts_file, output_dir, bands=2):
	# Decompose a time series object into a series of slices
	# Write the slices into isce unwrapped format.
	# We choose to use a single band for isce unwrapped format. 
	call(["mkdir","-p",output_dir],shell=False);
	tdata, xdata, ydata, zdata = netcdf_read_write.read_3D_netcdf(ts_file);
	for i in range(np.shape(zdata)[0]):
		call(["mkdir","-p",output_dir+"/scene_"+str(i)]);
		temp=zdata[i,:,:];

		# Write data out in isce format
		ny, nx = np.shape(temp);
		name = "ts_slice_"+str(i);
		filename = output_dir+"/scene_"+str(i)+"/"+name+".unw";
		temp=np.float32(temp);
		isce_read_write.write_isce_unw(temp, temp, nx, ny, "FLOAT", filename);

		isce_read_write.plot_scalar_data(filename, band=bands,colormap='rainbow',datamin=-50, datamax=200,
			aspect=1/5,outname=output_dir+"/scene_"+str(i)+"/isce_unw_band.png");
	return;