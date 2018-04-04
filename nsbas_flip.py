import numpy as np 
import matplotlib.pyplot as plt 
import scipy.io.netcdf as netcdf

# ------------- CONFIGURE ------------ # 
def configure():
	file_name="vel.nc";
	return file_name;

# ------------- INPUTS ------------ # 
def inputs(file_name):
	[xdata,ydata] = read_grd_xy(file_name);
	data = read_grd(file_name);
	return [xdata, ydata, data];

def read_grd(filename):
	data0 = netcdf.netcdf_file(filename,'r').variables['z'][::-1];
	data=data0.copy();
	return data;
def read_grd_xy(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['x'][::-1];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['y'][::-1];
	xdata=xdata0.copy();

	# This is the key! Flip the x-axis when necessary.  
	#xdata=np.flip(xdata,0);  # This is sometimes necessary and sometimes not!  Not sure why. 

	return [xdata, ydata]; 


# ----- OUTPUTS ----- # 
def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname):
	# # Write the netcdf velocity grid file.  This works, but doesn't get read from gmt grdinfo. Not sure why. 
	f=netcdf.netcdf_file(netcdfname,'w');
	f.createDimension('x',len(xdata));
	f.createDimension('y',len(ydata));
	print(np.shape(zdata));
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


if __name__=="__main__":
	file_name=configure();
	[xdata, ydata, data] = inputs(file_name);
	produce_output_netcdf(xdata, ydata, data, 'mm/yr','vel.grd');