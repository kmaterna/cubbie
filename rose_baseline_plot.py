# The purpose of this plot is to read baseline tables
# And give a visual tool to determine which are the best 
# one-year interferograms to make. 

import matplotlib.pyplot as plt 
import numpy as np 
import sentinel_utilities



def top_level_driver():
	baselinefile='baseline_table.dat'
	[stems, times, baselines, missiondays] = sentinel_utilities.read_baseline_table(baselinefile);
	rose_plot(times, baselines);
	return;

def rose_plot(times, baselines):

	t2015=[]; r2015=[]; d2015=[]; t2016=[]; r2016=[]; d2016=[]; t2017=[]; r2017=[]; d2017=[]; t2018=[]; r2018=[]; d2018=[];
	plt.figure();

	for i in range(len(times)):
		year=str(times[i])[0:4];
		day=str(times[i])[4:7];

		theta=2*np.pi*float(day)/365.25;
		radius=baselines[i]-min(baselines);		

		if year=='2015':
			r2015.append(radius);
			t2015.append(theta);
			d2015.append(day);
		elif year=='2016':
			r2016.append(radius);
			t2016.append(theta);
			d2016.append(day);
		elif year=='2017':
			r2017.append(radius);
			t2017.append(theta);
			d2017.append(day);
		else:
			r2018.append(radius);
			t2018.append(theta);
			d2018.append(day);
	
	dot1=plt.polar(t2015,r2015,'.',color='b',label='2015');
	dot2=plt.polar(t2016,r2016,'.',color='k',label='2016');
	dot3=plt.polar(t2017,r2017,'.',color='r',label='2017');	
	dot4=plt.polar(t2018,r2018,'.',color='g',label='2018');

	for i in range(len(d2017)):
		plt.annotate(d2017[i],xy=(t2017[i],r2017[i]),fontsize=6,color='r');
	for i in range(len(d2016)):
		plt.annotate(d2016[i],xy=(t2016[i],r2016[i]),fontsize=6,color='k');

	plt.legend(loc=4);
	plt.savefig('roseplot.eps');
	return;



