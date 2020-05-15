import numpy as np 
import sys, glob
from subprocess import call
import netcdf_read_write

# This is for when you've run a large SBAS in chunks of several million pixels each 
# Because it saves time to run in parallel. 

def get_input_dirs():
	input_files=[];
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/0_500000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/500000_1250000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/1250000_2000000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/2000000_3500000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/3500000_4500000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/4500000_5500000");
	input_files.append("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/5500000_7000000");
	return input_files;

def get_datestrs():
	files = glob.glob("/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/0_500000/*.grd");
	datestrs = [i.split('/')[-1][0:8] for i in files];
	print(datestrs);
	return datestrs;

def combine_all_files(datestr, input_dirs, output_dir):
	print("\nCombining files for date %s" % datestr);

	xdata, ydata, zdata0 = netcdf_read_write.read_grd_xyz(input_dirs[0]+"/"+datestr+".grd");
	xdata, ydata, zdata1 = netcdf_read_write.read_grd_xyz(input_dirs[1]+"/"+datestr+".grd");
	xdata, ydata, zdata2 = netcdf_read_write.read_grd_xyz(input_dirs[2]+"/"+datestr+".grd");
	xdata, ydata, zdata3 = netcdf_read_write.read_grd_xyz(input_dirs[3]+"/"+datestr+".grd");
	xdata, ydata, zdata4 = netcdf_read_write.read_grd_xyz(input_dirs[4]+"/"+datestr+".grd");
	xdata, ydata, zdata5 = netcdf_read_write.read_grd_xyz(input_dirs[5]+"/"+datestr+".grd");
	xdata, ydata, zdata6 = netcdf_read_write.read_grd_xyz(input_dirs[6]+"/"+datestr+".grd");
	zdata_total = np.zeros(np.shape(zdata0));
	
	for j in range(len(ydata)):
		if np.mod(j,200)==0:
			print(j)
		for k in range(len(xdata)):
			vector = [ zdata0[j][k], zdata1[j][k], zdata2[j][k], zdata3[j][k], zdata4[j][k], zdata5[j][k], zdata6[j][k] ];
			zdata_total[j][k]=np.sum(vector);
	output_file = output_dir+"/"+datestr+".grd";
	output_plot = output_dir+"/"+datestr+".png";
	netcdf_read_write.produce_output_netcdf(xdata, ydata, zdata_total, "mm", output_file);
	netcdf_read_write.produce_output_plot(output_file, datestr, output_plot, "mm", aspect=1.0, invert_yaxis=True, vmin=-50, vmax=100);
	return;

if __name__=="__main__":
	output_dir = "/Volumes/Ironwolf/Track_71/stacking/nsbas_apr20/combined/"
	call(["mkdir","-p",output_dir],shell=False);
	input_dirs=get_input_dirs();
	datestrs=get_datestrs();
	# for i in range(len(datestrs)):
		# combine_all_files(datestrs[i], input_dirs, output_dir);
	# In the morning, comment out the combine_all_files bit and just do the geocoding. 


	# In the morning, quickly geocode all the time series files. 
	for i in range(len(datestrs)):
		call(["quick_geocode.csh","stacking/nsbas_apr20/combined","merged",datestrs[i]+".grd",datestrs[i]+"_ll"],shell=False);



