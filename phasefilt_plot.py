import matplotlib.pyplot as plt 
import numpy as np
import glob as glob
import sys
import datetime as dt 
import netcdf_read_write

# TOP LEVEL DRIVER
def top_level_driver():
	[file_names]=configure();
	[xdata,ydata,data_all,date_pairs]=inputs(file_names);
	make_plots(xdata,ydata,data_all,date_pairs);
	return;


# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="intf_all";
	file_type="phasefilt.grd";
	
	file_names=glob.glob(file_dir+"/*/"+file_type);
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
	print("Reading "+str(len(file_names))+" files.");
	return [file_names];


# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = netcdf_read_write.read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:  # this happens to be in date order on my mac
		data = netcdf_read_write.read_grd(ifile);
		data_all.append(data);
	date_pairs=[];
	for name in file_names:
		pairname=name.split('/')[-2][0:15];
		date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
		print(pairname)
	return [xdata, ydata, data_all, date_pairs];


def make_plots(xdata,ydata,data_all,date_pairs):
	num_plots_x=4;
	num_plots_y=3;

	for i in range(len(data_all)):
		if np.mod(i,num_plots_y*num_plots_x)==0:
			count=i;

			fignum=i/(num_plots_y*num_plots_x); # counting figures up 0 to 1 to 2....

			# Looping forward and plotting the next 12 plots... 
			f,axarr = plt.subplots(num_plots_y, num_plots_x,figsize=(10,10));
			for k in range(num_plots_y):
				for m in range(num_plots_x):
					if count==len(data_all):
						break;

					# How many days separate this interferogram? 
					day1=date_pairs[count].split('_')[0];
					day2=date_pairs[count].split('_')[1];
					dt1=dt.datetime.strptime(day1,'%Y%j');
					dt2=dt.datetime.strptime(day2,'%Y%j');
					deltat=dt2-dt1;
					daysdiff=deltat.days;

					# The actual plotting
					axarr[k][m].imshow(data_all[count],cmap='jet',aspect=0.5);
					axarr[k][m].invert_yaxis();
					axarr[k][m].invert_xaxis();
					axarr[k][m].get_xaxis().set_ticks([]);
					axarr[k][m].get_yaxis().set_ticks([]);
					axarr[k][m].set_title(str(date_pairs[count])+'   '+str(daysdiff)+' days',fontsize=8);
					count=count+1;
			plt.savefig("selected_data_"+str(fignum)+".eps");
			plt.close();
	return;


if __name__=="__main__":
	top_level_driver();

