import matplotlib.pyplot as plt 
import numpy as np
import glob as glob
import sys
import datetime as dt
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

# -------- OUTPUT --------------- # 
def make_plots(xdata,ydata,data_all,date_pairs):
	num_plots_x=4;
	num_plots_y=3;

	for i in range(len(data_all)):
		if np.mod(i,num_plots_y*num_plots_x)==0:
			count=i;
			fignum=i/(num_plots_y*num_plots_x);
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

					numelements=np.shape(data_all[count]);
					mycorrs=np.reshape(data_all[count],(numelements[0]*numelements[1],1));
					nonans=mycorrs[~np.isnan(mycorrs)];
					above_threshold=mycorrs[np.where(mycorrs>0.2)];
					above_threshold=int(len(above_threshold)/1000);

					axarr[k][m].hist(nonans);

					axarr[k][m].set_title(str(date_pairs[count])+'   '+str(daysdiff)+' days',fontsize=8);
					axarr[k][m].set_yscale('log');
					axarr[k][m].set_xlim([0,1]);
					axarr[k][m].set_ylim([1,1000*1000]);
					axarr[k][m].plot([0.1,0.1],[0,1000*1000],'--r');
					axarr[k][m].text(0.75,200*1000,str(above_threshold)+'K',fontsize=8);

					count=count+1;
			plt.savefig("corr_hist_"+str(fignum)+".eps");
			plt.close();
	return;


if __name__=="__main__":
	top_level_driver();