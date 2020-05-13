# Utilities to convert file types

import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call
import glob
import isce_read_write
import netcdf_read_write


def gmtsar_nc_2_isce_stack(ts_file, output_dir):
	# Decompose a time series object into a series of slices
	# Write the slices into isce unwrapped format.
	call(["mkdir","-p",output_dir],shell=False);
	tdata, xdata, ydata, zdata = netcdf_read_write.read_3D_netcdf(ts_file);
	for i in range(np.shape(zdata)[0]):
		call(["mkdir","-p",output_dir+"/scene_"+str(i)]);
		temp=zdata[i,:,:];

		# Double check: lets plot everything, as it came in and as it went out. 
		plt.figure();
		plt.imshow(temp,aspect=1/5,vmin=-50, vmax=200,cmap='rainbow')
		plt.savefig(output_dir+"/scene_"+str(i)+"/ncdata.png");
		plt.close();

		# Write data out in isce format
		ny, nx = np.shape(temp);
		name = "ts_slice_"+str(i);
		filename = output_dir+"/scene_"+str(i)+"/"+name+".unw";
		isce_read_write.write_isce_unw(temp, nx, ny, "FLOAT", filename);

		# More inspecting
		test_data = isce_read_write.read_scalar_data(filename,band=2,flush_zeros=False);
		print("Nanmax of gmtsar:")
		print(np.nanmax(temp))
		print("Nanmax of isce:")
		print(np.nanmax(test_data))
		print("Nanmin of gmtsar:")
		print(np.nanmin(temp))
		print("Nanmin of isce:")
		print(np.nanmin(test_data))
		print("\n");
		isce_read_write.plot_scalar_data(filename, band=2,colormap='rainbow',datamin=-50, datamax=200,
			aspect=1/5,outname=output_dir+"/scene_"+str(i)+"/isce_band2.png");
	return;