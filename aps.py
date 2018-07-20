# Atmospheric Phase Screen method of Tymofyeyeva and Fialko, 2015


import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call 
import glob, sys, os
from copy import deepcopy
import datetime as dt 
import netcdf_read_write


def main_function(staging_directory, out_dir, rowref, colref, starttime, endtime, run_type):
	[file_names, width_of_stencil, n_iter] = configure(staging_directory, out_dir, run_type);
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names, starttime, endtime, run_type);
	[aps_array, corrected_intfs] = compute(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, n_iter, rowref, colref, run_type);
	outputs(xdata, ydata, data_all, corrected_intfs, aps_array, date_pairs, dates, out_dir);
	return;



# ------------- CONFIGURE ------------ # 
def configure(staging_directory, out_dir, run_type):
	width_of_stencil = 30;   # in days.  
	n_iter = 3;  # The number of times we're running the APS solver. 

	# Setting up the input and output directories. 
	if run_type=='test':
		file_names=glob.glob(staging_directory+"/????????_????????/unwrap*.grd");
	else:
		file_names=glob.glob(staging_directory+"/*_*_unwrap.grd");
	
	file_names=sorted(file_names);
	print("Performing APS on %d files in %s " % (len(file_names), staging_directory));
	if len(file_names)==0:
		print("Error! No files matching search pattern within "+staging_directory); sys.exit(1);
	call(['mkdir','-p',out_dir],shell=False);
	return [file_names, width_of_stencil, n_iter];


# ------------- INPUTS ------------ # 
def inputs(file_names, start_time, end_time, run_type):

	# Read the input grd files. Support for netcdf3 and netcdf4. 
	try:
		[xdata,ydata] = netcdf_read_write.read_grd_xy(file_names[0]);
	except TypeError:
		[xdata,ydata] = netcdf_read_write.read_netcdf4_xy(file_names[0]);

	data_all=[];
	date_pairs=[];
	dates=[];
	start_dt = dt.datetime.strptime(str(start_time),"%Y%m%d");
	end_dt = dt.datetime.strptime(str(end_time),"%Y%m%d");

	# Get the dates of the acquisitions from the file names. 
	for ifile in file_names:  # this happens to be in date order on my mac
		if run_type=='test':  # testing with Kang's format. "20171117_20171123/unwrap_new.grd"
			pairname=ifile.split('/')[-2];
			image1=pairname.split('_')[0];
			image2=pairname.split('_')[1];
			image1_dt = dt.datetime.strptime(image1,"%Y%m%d");
			image2_dt = dt.datetime.strptime(image2,"%Y%m%d");
		else:  #  the usual GMTSAR format
			pairname=ifile.split('/')[-1][0:15];
			image1=pairname.split('_')[0];
			image2=pairname.split('_')[1];
			image1_dt = dt.datetime.strptime(image1,"%Y%j");
			image2_dt = dt.datetime.strptime(image2,"%Y%j");
		
		if image1_dt>=start_dt and image1_dt<= end_dt:
			if image2_dt>=start_dt and image2_dt <= end_dt:
				try:
					data = netcdf_read_write.read_grd(ifile);
				except TypeError:
					data = netcdf_read_write.read_netcdf4(ifile);

				data_all.append(data);
				pairname=dt.datetime.strftime(image1_dt,"%Y%j")+'_'+dt.datetime.strftime(image2_dt,"%Y%j");
				date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
				dates.append(dt.datetime.strftime(image1_dt,"%Y%j"));
				dates.append(dt.datetime.strftime(image2_dt,"%Y%j"));

	dates=list(set(dates));
	dates=sorted(dates);
	print(date_pairs);
	print("Reading %d interferograms from %d acquisitions. " % (len(date_pairs), len(dates) ) );
	return [xdata, ydata, data_all, dates, date_pairs];



# ----------- COMPUTE ------------- #

def form_APS_pairs(date_pairs, mydate, width_of_stencil):
	# Date pairs: a global list of available interferograms (in order).
	# containing: a small list of interferograms that contain the image of interest
	# width_of_stencil: the number of days we tolerate in stencils

	print("Calculating APS stencil on "+mydate+":")

	pairlist = [];
	containing = [item for item in date_pairs if mydate in item];  # which interferograms contain the image of interest? 

	for item in containing: # for each interferogram containing the main image...
		image1=item.split('_')[0];
		image2=item.split('_')[1];
		if image1==mydate:
			continue;
		dt1 = dt.datetime.strptime(image1,"%Y%j");
		dt2 = dt.datetime.strptime(image2,"%Y%j");
		imwidth = (dt2 - dt1).days;
		if abs(imwidth)<width_of_stencil:  # for each valid interferogram that lands on the master image
			# See if there is a pair in the stack. Does it exist? 

			item_potential_match = dt.datetime.strftime(dt2+dt.timedelta(days=abs(imwidth)),"%Y%j");
			item_potential_pair = image2+'_'+item_potential_match;

			if item_potential_pair in containing:
				pair1=item;
				pair2=item_potential_pair;
				#print("%s found matching %s" % (item_potential_pair, item));
				index1=date_pairs.index(pair1);
				index2=date_pairs.index(pair2);
				pairlist.append([index1, index2])

	# Debugging statements
	daydiff = "12"
	for i in range(len(pairlist)):
		index1=pairlist[i][0];
		intf1=date_pairs[index1];
		d1=dt.datetime.strptime(intf1.split('_')[0],'%Y%j');
		d2=dt.datetime.strptime(intf1.split('_')[1],'%Y%j');
		daydiff=abs(d1-d2).days;
		print("---->"+intf1+"      [%d days]" % daydiff );
	for i in range(len(pairlist)):
		index2=pairlist[i][1];
		intf2=date_pairs[index2];
		d1=dt.datetime.strptime(intf2.split('_')[0],'%Y%j');
		d2=dt.datetime.strptime(intf2.split('_')[1],'%Y%j');
		daydiff=abs(d1-d2).days;
		print("     "+intf2+"<---- [%d days]" % daydiff );

	# print(pairlist);
	# Example of pairlist: [[63,65],[62,66]], where the numbers are the indeces of paired interferograms used for APS construction

	return pairlist; 


def compute_aps_image(dates, date_pairs, data_all, width_of_stencil, mydate, run_type):
	# For a given date and stack of interferograms, estimate an atmospheric phase screen through Common Scene Stacking. 
	# Be aware: the NANs will propagate through your network if you iterate too many times.  

	print("Computing APS for image %s " % mydate);
	zdim, rowdim, coldim = np.shape(data_all);
	my_aps=np.zeros((rowdim,coldim));

	# A list of valid interferogram pairs by index in the stack, which we use for making APS
	pairlist = form_APS_pairs(date_pairs, mydate, width_of_stencil);
	if len(pairlist)==0:
		print("ERROR! No valid APS stencils were detected in the stack for image %s" % mydate)
		return my_aps;

	for i in range(rowdim):
		for j in range(coldim):  # for each pixel...

			aps_temp_sum=0;
			pair_counter=0;

			for k in range(len(pairlist)):  # for each pair of interferograms that we can use for APS formation... 
				intf1_index=pairlist[k][0];
				intf2_index=pairlist[k][1];
				if ~np.isnan(data_all[intf1_index][i][j]) and ~np.isnan(data_all[intf2_index][i][j]):  
				# if both interferograms are numbers, then add to APS estimate. 
					pair_counter=pair_counter+1;
 
					aps_temp_sum = aps_temp_sum + (data_all[intf1_index][i][j] - data_all[intf2_index][i][j]);
			
			if pair_counter>0:
				my_aps[i][j]=(1.0/(2*pair_counter))*aps_temp_sum;
			else:
				my_aps[i][j]=np.nan;

	return my_aps;



def remove_aps(data_all, APS, date_pairs, dates):
	# In this function, you start with a set of interferograms (3D array), and correct it for a stack of APS. 
	# data_all: 3D array [n_intfs, xpix, ypix]
	# APS : 3D array     [n_images, xpix, ypix]
	# Return the updated interferograms. 

	print("Correcting intf stack for APS");
	zdim, rowdim, coldim = np.shape(data_all);
	corrected_intfs = deepcopy(data_all);

	# For each interferogram.
	for k in range(len(date_pairs)):
		dates_in_intf=date_pairs[k].split('_');
		ind1 = dates.index(dates_in_intf[0]);
		ind2 = dates.index(dates_in_intf[1]);
		dphi = APS[ind2][:][:]-APS[ind1][:][:];

		corrected_intfs[k][:][:]= np.subtract(data_all[k][:][:],dphi);

	return corrected_intfs;



def compute_ANC(APS_array):
	# APS_array = 3D array of numeric values
	# A scaled RMS of the atmospheric phase screen.
	print("Computing ANC from APS_array.");
	ANC_array=[];
	nsar = np.shape(APS_array)[0];
	for i in range(nsar):
		myaps1d = [];
		myaps = APS_array[i][:][:];
		for j in range(np.shape(myaps)[0]):
			for k in range(np.shape(myaps)[1]):
				if ~np.isnan(myaps[j][k]):
					myaps1d.append(myaps[j][k]);
		abar = np.mean(myaps1d);
		res_sq_array = [(a-abar)*(a-abar) for a in myaps1d];
		M=len(myaps1d);
		ANC = np.sqrt((1.0/M) * np.sum(res_sq_array) );
		ANC_array.append(ANC);
	maxANC= np.max(ANC_array);
	ANC_array=[10*(1.0/maxANC)*i for i in ANC_array];  # a normalization factor. 
	return ANC_array;


def get_initial_ANC(dates, date_pairs, data_all, width_of_stencil):
	# For iteration 0, starting in the middle. The exact numbers don't really matter, just the shape. 
	zdim = len(dates);
	ANC_array=np.zeros(zdim);
	intervals = 10/(np.ceil(zdim/2.0));
	middlei=np.floor(zdim/2);
	for i in range(zdim):
		ANC_array[i]=10-intervals*abs(i-middlei);
	return ANC_array; 



def compute(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, n_iter, rowref, colref, run_type):

	# Automated selection of dates and interferograms for APS construction. Pseudocode:
	# Produce ANCs from the original data, without altering the interferograms. 
	# Iterate several complete times throughout the stack of APS
	# Each time, iterate closer to a nice estimate of each image's APS. 
	# Produce an array of corrected interferogram data.
	# Run NSBAS. 

	nrows = len(ydata);
	ncols = len(xdata);
	APS_array = np.zeros((len(dates),nrows,ncols));
	ANC_array = get_initial_ANC(dates, date_pairs, data_all, width_of_stencil);  # Start from the middle. 
	make_APS_plot(APS_array,dates,'pythonIT0');

	# Recording the initial ANC values. 
	ofile=open("test_ANC0.txt",'w');
	for i in range(len(ANC_array)):
		ofile.write("%s %f\n" % (dates[i], ANC_array[i]) );
	ofile.close();

	# The main loop to iterate through common scene stacking. 
	for i in range(1,n_iter+1):

		print("######## Beginning iteration %d #########" % i);

		ordered_ANCs = [x for x,_ in sorted(zip(ANC_array, dates),reverse=True)]
		ordered_dates = [x for _,x in sorted(zip(ANC_array, dates), reverse=True)];

		# EACH DAY IN ORDERED ANC.  This loop is "Calculate_aps_linear.m"
		for j in range(len(ordered_dates)): 

			given_date=ordered_dates[j];  # THIS IS THE KEY
			print("Solving APS for day %s " % given_date);
			given_date_index = dates.index(given_date);
			APS_array[given_date_index][:][:] = np.zeros((nrows,ncols));  # Set the current day's APS to zero, because we're going to solve for it. 

			# Fix interferograms
			los_out_new = remove_aps(data_all, APS_array, date_pairs, dates);

			# Compute a single APS using corrected interferograms. 
			my_aps = compute_aps_image(dates, date_pairs, los_out_new, width_of_stencil, given_date, run_type);

			# Update the APS array. 
			APS_array[given_date_index][:][:]=my_aps;

		# At this point we have an updated stack of APS. 
		ANC_array = compute_ANC(APS_array);  # compute ANCs for the next iteration. 
		los_out_new = remove_aps(data_all, APS_array, date_pairs, dates);  # update the original intfs with the new APS. 

		make_APS_plot(APS_array,dates,'pythonIT'+str(i));

		# Recording the ANCs. 
		ofile=open("test_ANC"+str(i)+".txt",'w');
		for j in range(len(ANC_array)):
			ofile.write("%s %f\n" % (dates[j], ANC_array[j]) );
		ofile.close();

	return [APS_array, los_out_new];


# ----------- OUTPUTS ------------- # 


def display_many_intfs(data_all, out_dir, filename):

	vmin=-20;
	vmax=20;

	# Start by making a nice group of plots. 
	number_of_plots = np.shape(data_all)[0];
	sidelength = int(np.floor(np.sqrt(number_of_plots)));

	# display interferograms
	fig, ax = plt.subplots(sidelength,sidelength,sharey=True,sharex=True);
	for i in range(sidelength):
		for j in range(sidelength):
			ax[i,j].imshow(data_all[i*sidelength+j][:][:],aspect=0.5,cmap='hsv',vmin=vmin,vmax=vmax);
			plt.gca().invert_yaxis()
			plt.gca().invert_xaxis()
			plt.gca().get_xaxis().set_ticks([]);
			plt.gca().get_yaxis().set_ticks([]);

	plt.savefig(out_dir+'/'+filename+'.eps');
	plt.close();

	return;



def make_APS_plot(myaps_array, dates, filename):
	number_of_plots = np.shape(myaps_array)[0];
	sidelength = int(np.floor(np.sqrt(number_of_plots)));

	# display APS
	fig, ax = plt.subplots(sidelength,sidelength,sharey=True,sharex=True);
	for i in range(sidelength):
		for j in range(sidelength):
			ax[i,j].imshow(myaps_array[i*sidelength+j][:][:],aspect=0.5,cmap='YlGnBu',vmin=-4,vmax=4);
			plt.gca().get_xaxis().set_ticks([]);
			plt.gca().get_yaxis().set_ticks([]);
			ax[i,j].set_title(dates[i*sidelength+j],fontsize=9);			
	plt.savefig('../RESULTS_python/'+filename+'.eps');
	plt.close();
	return;


def display_myaps_array(myaps_array, dates, out_dir, filename):
	number_of_plots = np.shape(myaps_array)[0];
	sidelength = int(np.floor(np.sqrt(number_of_plots)));

	# display APS
	fig, ax = plt.subplots(sidelength,sidelength,sharey=True,sharex=True);
	for i in range(sidelength):
		for j in range(sidelength):
			ax[i,j].imshow(myaps_array[i*sidelength+j][:][:],aspect=0.5,cmap='hsv',vmin=-20,vmax=20);
			plt.gca().invert_yaxis()
			plt.gca().invert_xaxis()
			plt.gca().get_xaxis().set_ticks([]);
			plt.gca().get_yaxis().set_ticks([]);
			ax[i,j].set_title(dates[i*sidelength+j],fontsize=9);

	plt.savefig(out_dir+'/'+filename+'.eps');
	plt.close();

	return;

def intf_stack_to_grds(xdata, ydata, corrected_intfs, date_pairs, out_dir):
	print("Writing corrected interferograms to file. ")
	zdim, rowdim, coldim = np.shape(corrected_intfs);
	for i in range(zdim):
		item_name=date_pairs[i]
		netcdfname=out_dir+'/'+item_name+'_unwrap.grd';
		netcdf_read_write.produce_output_netcdf(xdata, ydata, corrected_intfs[i][:][:], 'unwrapped_phase', netcdfname);
		netcdf_read_write.flip_if_necessary(netcdfname);
	return;


def outputs(xdata, ydata, data_all, corrected_intfs, myaps_array, date_pairs, dates, out_dir):

	print("Starting the output process");

	# Display raw intfs
	display_many_intfs(data_all, out_dir,'raw_intfs');
	display_many_intfs(corrected_intfs, out_dir,'aps_intfs');
	display_myaps_array(myaps_array, dates, out_dir, 'APS_arrays');

	intf_stack_to_grds(xdata, ydata, corrected_intfs, date_pairs, out_dir);

	return;
