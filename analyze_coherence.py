# The purpose of this script is to run some diagnostics on the mean coherence. 
# How did the coherence vary with other things, like perpendicular baseline and season? 
# What do histograms of coherence look like? 
# A little messy, could use a little TLC, but passable for now. 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import sys
import glob as glob
from subprocess import call, check_output
import sentinel_utilities
import netcdf_read_write

# ------------- DRIVERS ------------------ # 

def analyze_coherence_function():
    corr_file='corr_results.txt';
    baseline_table = 'raw/baseline_table.dat';

    # Correlation vs. Other Things. 
    write_corr_results(corr_file);
    [stem1, stem2, mean_corr] = sentinel_utilities.read_corr_results(corr_file);
    [stems_blt,tbaseline,xbaseline,mission_days]=sentinel_utilities.read_baseline_table(baseline_table);
    make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline);

    # Histograms
	[file_names]=configure_histograms();
	[xdata,ydata,corr_all,date_pairs]=inputs(file_names);
	make_histograms(xdata,ydata,corr_all,date_pairs); 
    return;

# ------------ CORR vs STUFF -------------- # 

def write_corr_results(filename):
	ifile=open(filename,'w');

	dir_list=glob.glob('intf_all/2*'); # get a list of all the directories
	for item in dir_list:
		directname = item.split('/')[-1];  # format: 2015153_2015177
		call("gmt grdmath "+item+"/corr.grd MEAN = "+item+"/out.grd",shell=True);
		corr = check_output("gmt grdinfo "+item+"/out.grd | grep z | awk \'{print $3}\'", shell=True);
		corr = corr.split()[0];
		SLCs=check_output("ls "+item+"/*.SLC",shell=True); # THIS WILL BREAK IF USING PYTHON3.  USE PYTHON2. PROBLEM WITH BYTES vs STRING. 
		# IN CHECK OUTPUT, THE RETURN VALUE IS A BYTES OBJECT IN PYTHON3. OOPS. 
                # CAN I FIX THIS NOW? 
		slc1=SLCs.split('\n')[0];
		slc1=slc1.split('/')[-1];
		slc2=SLCs.split()[1];
		slc2=slc2.split('/')[-1];
		call("rm "+item+"/out.grd",shell=True);
		ifile.write('%s %s %s %s\n' % (directname, slc1, slc2, corr) )
		# Format: 2018005_2018017 S1A20180106_ALL_F1.SLC S1A20180118_ALL_F1.SLC 0.152738928795
	ifile.close();
	return;


def make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline):

	b_perp_baseline   = [];
	temporal_baseline = [];
	season = [];
	for i in range(len(mean_corr)):  # define a bunch of things for each interferogram.
		# How long in time? 
		temp1=stem1[i].split('_')[0];
		date1=dt.datetime.strptime(temp1[3:12],'%Y%m%d');
		temp2=stem2[i].split('_')[0];
		date2=dt.datetime.strptime(temp2[3:12],'%Y%m%d');
		temporal_baseline.append(abs((date1-date2).days));

		# Perpendicular Baseline?
		myindex1 = int(np.where(stems_blt == stem1[i])[0])
		myindex2 = int(np.where(stems_blt == stem2[i])[0])
		bl1=xbaseline[myindex1];
		bl2=xbaseline[myindex2];
		b_perp_baseline.append(abs(bl1-bl2));

		# Which season?
		season.append(int(date1.strftime("%j")));

	summer0=160;
	summer1=300;

	f, axarr = plt.subplots(3, figsize=(15, 15));
	axarr[0].set_title('Mean Coherence of Sentinel-1A Interferograms in Mendocino',fontsize=24)
	axarr[0].plot(temporal_baseline, mean_corr,'.',markersize=13);
	axarr[0].set_xlabel('Time (days)',fontsize=20)
	axarr[0].set_ylabel('Mean Coherence',fontsize=20)
	axarr[0].tick_params(axis='both',labelsize=20)

	for i in range(len(season)):
		if season[i]>summer0 and season[i]<summer1:
			line1, = axarr[1].plot(b_perp_baseline[i], mean_corr[i],'.',color='red',markersize=13, label='summer');
		else:
			line2, = axarr[1].plot(b_perp_baseline[i], mean_corr[i],'.',color='blue',markersize=13, label='not summer');
	axarr[1].set_xlabel('Baseline (m)',fontsize=20)
	axarr[1].set_ylabel('Mean Coherence',fontsize=20)
	axarr[1].tick_params(axis='both',labelsize=20)
	axarr[1].legend(handles=[line1, line2], loc=1, fontsize=18);


	for i in range(len(season)):
		if temporal_baseline[i]==12:
			line1, =axarr[2].plot(season[i], mean_corr[i],'.',color='magenta',markersize=13, label='12 days');
		else:
			line2, =axarr[2].plot(season[i], mean_corr[i],'.',color='green',markersize=13, label='>12 days');
	axarr[2].plot([summer0, summer0],[0.05, 0.35],'--k');
	axarr[2].plot([summer1, summer1],[0.05, 0.35],'--k');
	axarr[2].set_xlabel('Day of Year',fontsize=20)
	axarr[2].set_ylabel('Mean Coherence',fontsize=20)
	axarr[2].legend(handles=[line1, line2],loc=2, fontsize=18)	
	axarr[2].tick_params(axis='both',labelsize=20)
	plt.savefig('coherence_stats.eps')

	return;



# ------------- HISTOGRAMS ------------ # 

def configure_histograms():
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
def make_histograms(xdata,ydata,data_all,date_pairs):
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
	analyze_coherence_function();

