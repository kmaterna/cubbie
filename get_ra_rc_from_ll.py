# A set of python scripts
# What is the range/azimuth and row/col of a particular geographic coordinate? 
# This will use trans.dat, and some inputs of an example .grd file and geographic coords. 

import numpy as np 
import subprocess
import netcdf_read_write

def get_nearest_row_col(example_grd, ra, az):
	[xdata, ydata, zdata] = netcdf_read_write.read_grd_xyz(example_grd);	
	col_idx = (np.abs(xdata - ra)).argmin()  # xdata is columns
	row_idx = (np.abs(ydata - az)).argmin()  # ydata is rows
	return [row_idx, col_idx];


def get_ra_from_ll(trans_dat, lon, lat):
	# This works on a single point
	print("converting ll to ra")
	ofile=open("ll_temp.txt",'w');
	ofile.write("%f %f 0\n" % (lon, lat) );
	ofile.write("%f %f 0\n" % (lon+0.1, lat) );
	ofile.write("%f %f 0\n" % (lon+0.1, lat+0.1) );
	ofile.write("%f %f 0\n" % (lon, lat+0.1) );
	ofile.close();
	subprocess.call(['proj_ll2ra_ascii.csh',trans_dat,'ll_temp.txt','ra_temp.txt'],shell=False);
	subprocess.call(['rm','gmt.history'],shell=False);
	[ra, az, z] = np.loadtxt('ra_temp.txt',unpack=True)
	ra=ra[0]
	az=az[0]
	return [ra, az];


def get_ll_from_ra(trans_dat, ra, az):
	# Works on a single point
	print("converting ra to ll")
	ofile=open("ra_temp.txt",'w');
	ofile.write("%f %f 0\n" % (ra, az) );
	ofile.write("%f %f 0\n" % (ra+100, az+150) );
	ofile.write("%f %f 0\n" % (ra+150, az) );
	ofile.close();
	subprocess.call(['proj_ra2ll_ascii.csh',trans_dat,'ra_temp.txt','ll_temp.txt'],shell=False);
	subprocess.call(['rm','gmt.history'],shell=False);
	[lon, lat, z] = np.loadtxt('ll_temp.txt',unpack=True)
	lon=lon[0]
	lat=lat[0]
	return [lon, lat];


def get_ll_from_row_col(row, col, example_grd, trans_dat):
	[xdata, ydata, zdata] = netcdf_read_write.read_grd_xyz(example_grd)
	ra = xdata[col]; 
	az = ydata[row]; # check this
	[lon, lat] = get_ll_from_ra(trans_dat, ra, az);
	return [lon, lat];






if __name__=="__main__":
	trans_dat="topo/trans.dat";
	example_grd="stacking/unwrapped/2015157_2018177_unwrap.grd"
	lon,lat = -116.572, 35.321;  # P617
	ra = 10296.8986328
	az = 5984.77251953
	[ra, az] = get_ra_from_ll(trans_dat, lon, lat);
	[lon1, lat1] = get_ll_from_ra(trans_dat, ra, az);
	[row, col] = get_nearest_row_col(example_grd, ra, az);
	print(row, col)



