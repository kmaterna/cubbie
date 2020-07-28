# The purpose of this script is to choose a reference pixel
# Criteria:
# 1. High coherence in pretty much all interferograms- will give warning if you don't have a pixel with perfect coherence. 
# 2. Located in the un-deforming region of the field (if that exists). Target region defined manually. 
# 3. Likely to unwrap correctly

# The code will find all unwrap.grd files located in file_dir

import numpy as np 
import matplotlib.pyplot as plt 
import glob, sys, math
import netcdf_read_write


def main_function(file_dir):
	print("Starting to choose a reference pixel")
	[filenames, range_bounds, azimuth_bounds]=configure(file_dir);
	[xdata, ydata, data_all, dates, date_pairs] = inputs(filenames);
	[number_of_datas, zdim] = analyze_coherent_number(data_all);
	[rowref, colref] = pick_reference_pixel(range_bounds, azimuth_bounds, number_of_datas, zdim);
	outputs(rowref, colref, range_bounds, azimuth_bounds, xdata, ydata, number_of_datas);
	return [rowref, colref];


def configure(file_dir):
	# Setting up the input and output directories. 
	# Range: columns
	# Azimuth: rows
	# file_dir is usually something like 'intf_all/unwrap.grd'
	file_of_interest='unwrap.grd'
	file_names=glob.glob(file_dir+"/*_*_"+file_of_interest);
	if len(file_names)==0:
		print("Error! No files matching search pattern with "+file_of_interest); sys.exit(1);
	range_bounds = [150, 190]
	azimuth_bounds = [235, 270]
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
	[zdim, rowdim, coldim] = np.shape(zdata)
	# these 2D arrays are indexed by row-column convention. 
	number_of_datas=np.zeros([rowdim,coldim]);
	for k in range(zdim):
		for i in range(rowdim):
			for j in range(coldim):
				if not math.isnan(zdata[k][i][j]):
					number_of_datas[i][j]=number_of_datas[i][j]+1;

	print("Maximum Possible Interferograms = %d. " % (zdim) );
	best_pixels=0;
	for i in range(rowdim):
		for j in range(coldim):
			if number_of_datas[i][j]==zdim:
				best_pixels=best_pixels+1;
	print("Number of pixels equal to max_possible: %d out of %d ( %f %%)" % ( best_pixels, rowdim*coldim, 100.0*best_pixels/(rowdim*coldim) ) );

	return [number_of_datas, zdim];


def pick_reference_pixel(range_bounds, azimuth_bounds, number_of_datas, maxnum):
	rowref=0;
	colref=0;
	counting_maxnum=0;
	for i in range(azimuth_bounds[0], azimuth_bounds[1]):
		for j in range(range_bounds[0], range_bounds[1]):
			counting_maxnum=max(counting_maxnum, number_of_datas[i][j]);
			if number_of_datas[i][j]==maxnum:
				rowref=i;
				colref=j;
				return [rowref, colref];
	if counting_maxnum<maxnum:
		print("WARNING! There is no pixel within range with %d coherent images!" % maxnum);
		print("The maximum number of coherent pixels is: %d " % counting_maxnum);
		print("Proceeding to select reference pixel with coherence in %d images." % counting_maxnum);
	# if there is no pixel that's perfectly coherent in the box: 
	if rowref==0 and colref==0:
		for i in range(azimuth_bounds[0], azimuth_bounds[1]):
			for j in range(range_bounds[0], range_bounds[1]):
				if number_of_datas[i][j]==counting_maxnum:
					rowref=i;
					colref=j;
					return [rowref, colref];
	return [rowref, colref];



# -------------- OUTPUTS -------------- # 
def outputs(rowref, colref, range_bounds, azimuth_bounds, xdata, ydata, number_of_datas):

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
	plt.plot([range_bounds[0], range_bounds[1],range_bounds[1], range_bounds[0], range_bounds[0] ],[azimuth_bounds[0], azimuth_bounds[0], azimuth_bounds[1], azimuth_bounds[1], azimuth_bounds[0]], color='k')
	plt.plot(colref, rowref, '.k');
	plt.savefig('number_of_datas.eps');
	plt.close();

	print("Reference pixel is: (%d, %d) (rowref, colref)" % (rowref, colref) );
	print("Number of coherent images for reference pixel: %d" % (number_of_datas[rowref][colref]) ); 
	return;



if __name__=="__main__":
	[rowref, colref] = main_function('intf_all/unwrap.grd');






