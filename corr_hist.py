import matplotlib.pyplot as plt 
import numpy as np
import glob as glob
import sys
import netcdf_read_write


def top_level_driver():
	[file_names]=configure();
	[xdata,ydata,corr_all,date_pairs]=inputs(file_names);
	make_plots(xdata,ydata,corr_all,date_pairs);
	return;

# ------------- CONFIGURE ------------ # 
def configure():
	file_dir="intf_all";
	file_type="corr.grd";
	
	file_names=glob.glob(file_dir+"/*/"+file_type);
	if len(file_names)==0:
		print("Error! No files matching search pattern."); sys.exit(1);
	print("Reading "+str(len(file_names))+" files.");
	return [file_names[0:-2]];


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

# -------- OUTPUT --------------- # 
def make_plots(xdata,ydata,data_all,date_pairs):
	num_plots_x=4;
	num_plots_y=3;

	for i in range(len(data_all)):
		if np.mod(i,num_plots_y*num_plots_x)==0:
			count=i;
			fignum=i/(num_plots_y*num_plots_x);
			f,axarr = plt.subplots(num_plots_y, num_plots_x,figsize=(10,10));
			for i in range(num_plots_y):
				for j in range(num_plots_x):
					if count==len(data_all)-1:
						break;
					numelements=np.shape(data_all[count]);
					mycorrs=np.reshape(data_all[count],(numelements[0]*numelements[1],1));
					nonans=mycorrs[~np.isnan(mycorrs)];
					axarr[i][j].hist(nonans);
					axarr[i][j].set_title(str(date_pairs[count]),fontsize=8);
					axarr[i][j].set_yscale('log');
					axarr[i][j].set_xlim([0,1]);
					axarr[i][j].set_ylim([1,1000*1000]);
					axarr[i][j].plot([0.1,0.1],[0,1000*1000],'--r');

					count=count+1;
			plt.savefig("corr_hist_"+str(fignum)+".eps");
			plt.close();
	return;


if __name__=="__main__":
	top_level_driver();