# The purpose of this script is to see how unwrapping has proceeded. 
# 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import glob as glob
from subprocess import call, check_output
import sentinel_utilities


def analyze_coherence_function():
    corr_file='corr_results.txt';
    baseline_table = 'raw/baseline_table.dat';
    write_corr_results(corr_file);
    [stem1, stem2, mean_corr] = sentinel_utilities.read_corr_results(corr_file);
    [stems_blt,tbaseline,xbaseline,mission_days]=sentinel_utilities.read_baseline_table(baseline_table);
    make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline);
    return;

def write_corr_results(filename):
	ifile=open(filename,'w');

	dir_list=glob.glob('intf_all/2*'); # get a list of all the directories
	for item in dir_list:
		directname = item.split('/')[-1];  # format: 2015153_2015177
		call("gmt grdmath "+item+"/corr.grd MEAN = "+item+"/out.grd",shell=True);
		corr = check_output("gmt grdinfo "+item+"/out.grd | grep z | awk \'{print $3}\'", shell=True);
		corr = corr.split()[0];
		SLCs=check_output("ls "+item+"/*.SLC",shell=True);
		slc1=SLCs.split()[0];
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



if __name__=="__main__":
	analyze_coherence_function();

