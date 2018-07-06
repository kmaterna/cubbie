
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import datetime as dt 
import collections
import sys
import netcdf_read_write
import aps



def how_many_nans(filename):
	[xdata, ydata, zdata]=netcdf_read_write.read_grd_xyz(filename);
	nan_pixels = np.count_nonzero(np.isnan(zdata));
	total_pixels=np.shape(zdata)[0]*np.shape(zdata)[1];
	print("For file %s: %d pixels of %d are NaNs (%f percent)." % (filename, nan_pixels, total_pixels, 100*float(nan_pixels/float(total_pixels))) )
	return [nan_pixels, total_pixels];

def number_below_value(filename, value):
	[xdata, ydata, zdata]=netcdf_read_write.read_grd_xyz(filename);
	count=0;
	for i in range(len(ydata)):
		for j in range(len(xdata)):
			if zdata[i][j] < value:
				count=count+1;
	total_pixels=np.shape(zdata)[0]*np.shape(zdata)[1];
	print("For file %s: %d pixels of %d are below %f (%f percent)." % (filename, count, total_pixels, value, 100*float(count/float(total_pixels))) )
	return;

def plot_grid_file(filename, figname):
	[xdata, ydata, zdata]=netcdf_read_write.read_grd_xyz(filename);
	plt.figure();
	plt.imshow(zdata);
	plt.colorbar();
	plt.savefig(figname+'.eps');
	plt.close();
	return;

def plot_grid_data(griddata, figname):
	plt.figure();
	plt.imshow(griddata);
	plt.colorbar();
	plt.savefig(figname+'.eps');
	plt.close();
	return;


def process_by_boxes(xc, yc, zc, xp, yp, zp, xdec, ydec, value):
	# c means correlation, p means phase
	xdec=int(xdec);
	ydec=int(ydec);
	number_of_xboxes=int(np.ceil(len(xc)/xdec));
	number_of_yboxes=int(np.ceil(len(yc)/ydec));
	newx    =np.zeros([number_of_xboxes]);
	newy    =np.zeros([number_of_yboxes]);
	numabove=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array
	numbelow=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array

	for i in range(number_of_xboxes):
		for j in range(number_of_yboxes):

			# Go through each box of (xdec, ydec) pixels and find the max correlated pixel. 
			min_x=i*xdec;  # start counting at 0
			max_x=min(i*xdec+xdec,len(xc));  # end counting at 0+xdec, or the end of the array
			min_y=j*ydec;
			max_y=min(j*ydec+ydec,len(yc));
			newx[i]=np.median(xc[min_x:max_x]);
			newy[j]=np.median(yc[min_y:max_y]);			

			myzc = zc[min_y:max_y,min_x:max_x];  # correlation
			#myzp = zp[min_y:max_y,min_x:max_x];  # phase.     
			# Numpy note: zp[0:2][0:2] is not the same thing as zp[0:2,0:2];
			

			# The math itself: 
			# Get the number of boxes above and below the cutoff value inside of the box. 

			number_of_nans = np.count_nonzero(np.isnan(myzc));
			if (xdec*ydec==number_of_nans):
				number_above=np.nan;
				number_below=np.nan;
			else:
				number_above=0;
				number_below=0;
				for n in range(xdec):
					for m in range(ydec):
						if myzc[m][n]>=value:
							number_above=number_above+1;
						else:
							number_below=number_below+1;  # below value, or nan
			numabove[i][j] = number_above;
			numbelow[i][j] = number_below;

	numabove=numabove.T;
	numbelow=numbelow.T;

	return [newx, newy, numabove, numbelow];



def remove_nans_array(myarray):
	numarray=[];
	for i in range(len(myarray)):
		if ~np.isnan(myarray[i]):
			numarray.append(myarray[i][0]);
	return numarray;



def develop_mean_phase(phase_array):
	# This function takes an array of phase values, and determines the mean phase value (sensitive to cycle slips)
	# It uses the "mean of circular quantities" technique. 
	xarray= [np.cos(i) for i in phase_array];
	yarray= [np.sin(i) for i in phase_array];
	xmean=np.mean(xarray);
	ymean=np.mean(yarray);
	tolerance=0.00001;
	if abs(xmean)<tolerance and abs(ymean)<tolerance:
		print("Error! The mean phase is undefined!");
		sys.exit(0);
	else:
		meanphase=np.arctan2(ymean,xmean);
		# Mean phase already wrapped into the -pi to pi range. 
	return meanphase;


def develop_median_phase(phase_array):
	# This function takes an array of phase values, and determines the median phase value (sensitive to cycle slips)
	# It uses the "median of circular quantities" technique. 
	xarray= [np.cos(i) for i in phase_array];
	yarray= [np.sin(i) for i in phase_array];
	xmedian=np.median(xarray);
	ymedian=np.median(yarray);
	tolerance=0.00001;
	if abs(xmedian)<tolerance and abs(ymedian)<tolerance:
		print("Error! The median phase is undefined!");
		sys.exit(0);
	else:
		medianphase=np.arctan2(ymedian,xmedian);
	return medianphase;






# Step 1: Figure out how many nan's are in a netcdf file. 
# And how many correlation values are below a certain value. 
folder='2018017_2018029'
filename='intf_all/'+folder+'/phasefilt.grd'
# [nanpixels, totalpixels] = how_many_nans(filename);
# plot_grid_file(filename,'phasefilt');



# filename='intf_all/'+folder+'/unwrap.grd'
# [nanpixels, totalpixels] = how_many_nans(filename);
# plot_grid_file(filename,'unwrap');
# filename='intf_all/'+folder+'/phasefilt_full.grd'
# [nanpixels, totalpixels] = how_many_nans(filename);
# plot_grid_file(filename,'phasefilt_full');
# filename='intf_all/'+folder+'/unwrap.grd'
# [nanpixels, totalpixels] = how_many_nans(filename);

# Investigate coherence
# filename='intf_all/'+folder+'/corr.grd';
# number_below_value(filename, 0.1);
# filename='intf_all/'+folder+'/corr_full.grd';
# number_below_value(filename, 0.1);
#[nanpixels, totalpixels] = how_many_nans(filename);


# Step 2: Figure out how many pixels in each box are above 0.1
# Figure out how many pixels in each box are below 0.1. 
# filename='intf_all/'+folder+'/corr_full.grd';
# [xc, yc, zc]=dec_corrbased.read_grd(filename);
# filename='intf_all/'+folder+'/phasefilt_full.grd';
# [xp, yp, zp]=dec_corrbased.read_grd(filename);




# [newx, newy, numabove, numbelow] = process_by_boxes(xc, yc, zc, xp, yp, zp, xdec, ydec, 0.1);
# plot_grid_data(numabove, 'numabove');
# plot_grid_data(numbelow, 'numbelow');

# Just investigate a couple of grid files to look at what's inside. 
# plot_grid_file('intf_all/'+folder+'/mask_full.grd','mask_full');
# plot_grid_file('intf_all/'+folder+'/mask.grd','mask_decimated');
# [nanpixels, totalpixels] = how_many_nans('intf_all/'+folder+'/mask_full.grd');
# [nanpixels, totalpixels] = how_many_nans('intf_all/'+folder+'/mask.grd');

#How many pixels are decorrelated? 
# counter=0;
# tolerance=0;
# for i in range(len(newx)):
# 	for j in range(len(newy)):
# 		if numbelow[j][i]>=(xdec*ydec)-tolerance:
# 			counter=counter+1;

# print("%d out of %d blocks have %d coherent pixels" % (counter, len(newx)*len(newy), tolerance) );


# APS TESTING
staging_directory='intf_all/referenced_unwrap.grd'
out_dir='intf_all/aps_unwrap.grd'
aps.main_function(staging_directory, out_dir);


