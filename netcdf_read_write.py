# Netcdf reading and writing functions
# Bring a netcdf3 file into python!

import numpy as np 
import scipy.io.netcdf as netcdf
import matplotlib.pyplot as plt 
import subprocess


# --------------- READING ------------------- # 

def read_grd(filename):
	data0 = netcdf.netcdf_file(filename,'r').variables['z'][::];
	data=data0.copy();
	return data;
def read_grd_xy(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['x'][:];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['y'][:];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	return [xdata, ydata]; 

def read_grd_xyz(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['x'][:];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['y'][:];
	zdata0 = netcdf.netcdf_file(filename,'r').variables['z'][::];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	zdata=zdata0.copy();
	return [xdata, ydata, zdata]; 

def read_grd_lonlatz(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['lon'][:];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['lat'][:];
	zdata0 = netcdf.netcdf_file(filename,'r').variables['z'][::];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	zdata=zdata0.copy();
	return [xdata, ydata, zdata]; 

def read_grd_variables(filename, var1, var2, var3):
	xdata0= netcdf.netcdf_file(filename,'r').variables[var1][:];
	xdata=xdata0.copy();
	ydata0= netcdf.netcdf_file(filename,'r').variables[var2][:];
	ydata=ydata0.copy();
	zdata0= netcdf.netcdf_file(filename, 'r').variables[var3][::];
	zdata=zdata0.copy();
	return [xdata, ydata, zdata];

def read_netcdf4_xy(filename):
	netcdf4file=filename;
	netcdf3file=filename+'nc3';
	subprocess.call('nccopy -k classic '+netcdf4file+' '+netcdf3file,shell=True); 
	[xdata, ydata] = read_grd_xy(netcdf3file);
	return [xdata, ydata];

def read_netcdf4(filename):
	netcdf4file=filename;
	netcdf3file=filename+'nc3';
	subprocess.call('nccopy -k classic '+netcdf4file+' '+netcdf3file,shell=True); 
	data = read_grd(netcdf3file);
	return data;

def read_netcdf4_xyz(filename):
	netcdf4file=filename;
	netcdf3file=filename+'nc3';
	subprocess.call('nccopy -k classic '+netcdf4file+' '+netcdf3file,shell=True); 
	zdata = read_grd(netcdf3file);
	[xdata, ydata] = read_grd_xy(netcdf3file);
	return [xdata, ydata, zdata];

def read_netcdf4_variables(filename, var1, var2, var3):
	netcdf4file=filename;
	netcdf3file=filename+'nc3';
	subprocess.call('nccopy -k classic '+netcdf4file+' '+netcdf3file,shell=True); 
	[xdata, ydata, zdata] = read_grd_variables(filename, var1, var2, var3);
	return [xdata, ydata, zdata];


def read_any_grd_xyz(filename):
	# Switch between netcdf4 and netcdf3 automatically. 
	try:
		[xdata, ydata, zdata] = read_grd_xyz(filename);
	except TypeError: 
		[xdata, ydata, zdata] = read_netcdf4_xyz(filename);
	return [xdata, ydata, zdata];


def read_any_grd_variables(filename, var1, var2, var3):
	# Switch between netcdf4 and netcdf3 automatically. 
	try:
		[xdata, ydata, zdata] = read_grd_variables(filename, var1, var2, var3);
	except TypeError: 
		[xdata, ydata, zdata] = read_netcdf4_variables(filename, var1, var2, var3);
	return [xdata, ydata, zdata];

def give_metrics_on_grd(filename):
	grid_data = read_grd(filename);
	nan_pixels = np.count_nonzero(np.isnan(grid_data));
	total_pixels=np.shape(grid_data)[0]*np.shape(grid_data)[1];
	print("Shape of %s is [%d, %d]" % (filename, np.shape(grid_data)[0], np.shape(grid_data)[1]) );
	print("Min data is %f " % (np.nanmin(grid_data)) );
	print("Max data is %f " % (np.nanmax(grid_data)) );
	print("Nans: %d of %d pixels are nans (%.3f percent)" % (nan_pixels, total_pixels, nan_pixels/total_pixels*100) );
	return; 


# --------------- WRITING ------------------- # 

def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname):
	# # Write the netcdf velocity grid file.  
	f=netcdf.netcdf_file(netcdfname,'w');
	f.history = 'Created for a test';
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
	flip_if_necessary(netcdfname);
	return;

def flip_if_necessary(filename):
	# IF WE NEED TO FLIP DATA:
	xinc = subprocess.check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $8}\'',shell=True);  # the x-increment
	yinc = subprocess.check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $9}\'',shell=True);  # the x-increment
	xinc=float(xinc.split()[0]);
	yinc=float(yinc.split()[0]);

	if xinc < 0:  # FLIP THE X-AXIS
		print("flipping the x-axis");
		[xdata,ydata] = read_grd_xy(filename);
		data = read_grd(filename);
		# This is the key! Flip the x-axis when necessary.  
		#xdata=np.flip(xdata,0);  # This is sometimes necessary and sometimes not!  Not sure why. 
		produce_output_netcdf(xdata, ydata, data, 'mm/yr',filename);
		xinc = subprocess.check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $8}\'',shell=True);  # the x-increment
		xinc = float(xinc.split()[0]);
		print("New xinc is: %f " % (xinc) );
	if yinc < 0:
		print("flipping the y-axis");
		[xdata,ydata] = read_grd_xy(filename);
		data = read_grd(filename);
		# Flip the y-axis when necessary.  
		# ydata=np.flip(ydata,0);  
		produce_output_netcdf(xdata, ydata, data, 'mm/yr',filename);
		yinc = subprocess.check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $9}\'',shell=True);  # the x-increment
		yinc = float(yinc.split()[0]);
		print("New yinc is: %f" % (yinc) );
	return;


def produce_output_plot(netcdfname, plottitle, plotname, cblabel):
	# Read in the dataset
	[xread, yread, zread] = read_grd_xyz(netcdfname);

	# Make a plot
	fig = plt.figure(figsize=(7,10));
	ax1 = fig.add_axes([0.0, 0.1, 0.9, 0.8]);
	plt.imshow(zread,aspect=1.0);
	plt.gca().invert_yaxis()  # for imshow, rows get labeled in the downward direction
	# plt.gca().get_xaxis().set_ticks([]);
	# plt.gca().get_yaxis().set_ticks([]);
	plt.title(plottitle);
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label(cblabel, size=16);
	plt.savefig(plotname);
	plt.close();
	return;


def produce_output_contourf(netcdfname, plottitle, plotname, cblabel):
	# Read in the dataset
	[xread, yread, zread] = read_grd_xyz(netcdfname);

	# Make a plot
	fig = plt.figure(figsize=(7,10));
	plt.contourf(xread, yread, zread)
	# plt.gca().get_xaxis().set_ticks([]);
	# plt.gca().get_yaxis().set_ticks([]);

	plt.title(plottitle);
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label(cblabel, size=16);
	plt.savefig(plotname);
	plt.close();
	return;	

