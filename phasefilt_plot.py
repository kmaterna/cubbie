import matplotlib.pyplot as plt 
import numpy as np
import glob as glob
import netcdf_read_write


# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="phasefilt";
	
	file_names=glob.glob(file_dir+"/*.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
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
		pairname=name.split('/')[-1][0:15];
		date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
	return [xdata, ydata, data_all, date_pairs];

def make_plots(xdata,ydata,data_all,date_pairs):
	num_plots_x=4;
	num_plots_y=3;
	count=40;
	plt.figure();
	f,axarr = plt.subplots(num_plots_y, num_plots_x);
	for i in range(num_plots_y):
		for j in range(num_plots_x):
			axarr[i][j].imshow(data_all[count]);
			axarr[i][j].invert_yaxis();
			axarr[i][j].invert_xaxis();
			axarr[i][j].get_xaxis().set_ticks([]);
			axarr[i][j].get_yaxis().set_ticks([]);
			axarr[i][j].set_title(str(date_pairs[count]),fontsize=10);
			count=count+1;
	plt.savefig("selected_data.eps");
	plt.close();
	return;


if __name__=="__main__":
	[file_names]=configure();
	[xdata,ydata,data_all,date_pairs]=inputs(file_names);
	make_plots(xdata,ydata,data_all,date_pairs);

