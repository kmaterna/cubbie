# The purpose of this script is to see how unwrapping has proceeded. 
# 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import sentinel_utilities


def analyze_unwrapping_function():
    [stem1, stem2, mean_corr] = sentinel_utilities.read_corr_results('corr_results.txt');
    [stems_blt,tbaseline,xbaseline,mission_days]=sentinel_utilities.read_baseline_table('raw/baseline_table.dat');
    make_coh_vs_others_plots(stem1, stem2, mean_corr, stems_blt, xbaseline);
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
	axarr[0].set_title('Mean Coherence of Sentinel-1A Interferograms',fontsize=24)
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
	analyze_unwrapping_function();

