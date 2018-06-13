# The purpose of this script is to choose a reference pixel
# Criteria:
# 1. High coherence in pretty much all interferograms
# 2. Located in the un-deforming region of the field (if that exists)
# 3. Likely to unwrap correctly

# The code will find all unwrap.grd files located in intf_all/unwrap.grd

import numpy as np 
import matplotlib.pyplot as plt 
import glob, sys, math
import netcdf_read_write



def configure():
	# Setting up the input and output directories. 
	file_of_interest='unwrap.grd'
	file_dir="intf_all/"+file_of_interest;
	file_names=glob.glob(file_dir+"/*_*_"+file_of_interest);
	if len(file_names)==0:
		print("Error! No files matching search pattern with "+file_of_interest); sys.exit(1);
	range_bounds = [160, 180]
	azimuth_bounds = [245, 260]
	return [file_names, range_bounds, azimuth_bounds];


# ------------- INPUTS ------------ # 
def inputs(file_names):
	[xdata,ydata] = netcdf_read_write.read_grd_xy(file_names[0]);
	data_all=[];
	for ifile in file_names:  # this happens to be in date order on my mac
		data = netcdf_read_write.read_grd(ifile);
		data_all.append(data);
	date_pairs=[];
	dates=[];
	for name in file_names:
		pairname=name.split('/')[-1][0:15];
		date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
		splitname=pairname.split('_');
		dates.append(splitname[0])
		dates.append(splitname[1])
	dates=list(set(dates));
	dates=sorted(dates);
	print(dates);
	print("%d interferograms read from %d acquisitions." % (len(date_pairs), len(dates)) );
	return [xdata, ydata, data_all, dates, date_pairs];


# ----------- COMPUTE ------------ # 
def analyze_coherent_number(zdata):
	# Analyze the number of coherent acquisitions for each pixel
	print("Analyzing the number of coherent interferograms per pixel.")
	[zdim, xdim, ydim] = np.shape(zdata)
	number_of_datas=np.zeros([xdim,ydim]);
	for k in range(zdim):
		for i in range(xdim):
			for j in range(ydim):
				if not math.isnan(zdata[k][i][j]):
					number_of_datas[i][j]=number_of_datas[i][j]+1;

	print("Maximum Possible Interferograms = %d. " % (len(date_pairs)) );
	best_pixels=0;
	for i in range(xdim):
		for j in range(ydim):
			if number_of_datas[i][j]==len(date_pairs):
				best_pixels=best_pixels+1;
	print("Number of pixels equal to max_possible: %d out of %d ( %f %%)" % ( best_pixels, xdim*ydim, 100.0*best_pixels/(xdim*ydim) ) );

	return [number_of_datas, zdim];


def pick_reference_pixel(range_bounds, azimuth_bounds, number_of_datas, maxnum):
	xref=0;
	yref=0;
	counting_maxnum=0;
	for i in range(range_bounds[0], range_bounds[1]):
		for j in range(azimuth_bounds[0], azimuth_bounds[1]):
			counting_maxnum=max(counting_maxnum, number_of_datas[i][j]);
			if number_of_datas[i][j]==maxnum:
				xref=i;
				yref=j;
				return [xref, yref];
	if xref==0 and yref==0:
		for i in range(range_bounds[0], range_bounds[1]):
			for j in range(azimuth_bounds[0], azimuth_bounds[1]):
				if number_of_datas[i][j]==counting_maxnum:
					xref=i;
					yref=j;
					return [xref, yref];
	return [xref, yref];



# -------------- OUTPUTS -------------- # 
def outputs(xref, yref, range_bounds, azimuth_bounds, xdata, ydata, number_of_datas):

	[xdim, ydim]=np.shape(number_of_datas);

	zdata=np.reshape(number_of_datas, [xdim*ydim, 1])
	plt.figure();
	plt.hist(zdata,bins=81);
	plt.ylabel('Number of Pixels');
	plt.xlabel('Number of Coherent Interferograms')
	plt.savefig('Number_of_pixels.eps');
	plt.close();


	plt.figure();
	plt.imshow(number_of_datas);
	plt.plot([range_bounds[0], range_bounds[1] ],[azimuth_bounds[0], azimuth_bounds[1]], color='k')
	plt.plot(xref, yref, '.k');
	plt.savefig('number_of_datas.eps');
	plt.close();

	print("Reference pixel is: (%d, %d)" % (xref, yref) );
 
	return;




if __name__=="__main__":
	print("Starting to choose a reference pixel")
	
	[filenames, range_bounds, azimuth_bounds]=configure();
	[xdata, ydata, data_all, dates, date_pairs] = inputs(filenames);
	[number_of_datas, zdim] = analyze_coherent_number(data_all);
	[xref, yref] = pick_reference_pixel(range_bounds, azimuth_bounds, number_of_datas, zdim);
	outputs(xref, yref, range_bounds, azimuth_bounds, xdata, ydata, number_of_datas);






