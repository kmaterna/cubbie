# A series of loops that overwrite the files they started with. 
# Ex1: Multiply a bunch of grids by -1 to match someone else's velocity format
# Ex2: Multiply a bunch of grids by a nanmask


import numpy as np 
import matplotlib.pyplot as plt 
import netcdf_read_write
import subprocess


def multiply_stack_by_mask(filelist, maskfile):
	# Hasn't been tested
	for file in filelist:
		subprocess.call(['gmt','grdmath',file,maskfile,'MUL','-G'+file],shell=False);
	return;

def multiply_by_minus1_multiple(filelist):
	# Hasn't been tested
	for file in filelist:
		subprocess.call(['gmt','grdmath',file,'-1','MUL','-G'+file],shell=False);
	return;	

def multiply_file_by_minus1(filename, new_filename):
	print("multiplying %s by -1 " % filename);
	x, y, z = netcdf_read_write.read_netcdf4_variables(filename,'lon','lat','z');
	z = np.multiply(z, -1);
	filestem = filename.split('.grd')[0];
	netcdf_read_write.produce_output_netcdf(x, y, z, "mm/yr", new_filename,dtype=np.float32);
	return;
