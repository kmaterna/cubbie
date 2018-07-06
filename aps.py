# Atmospheric Phase Screen method of Tymofyeyeva and Fialko, 2015


import numpy as np 
import matplotlib.pyplot as plt 
from subprocess import call 
import glob, sys, os
from copy import deepcopy
import datetime as dt 
import netcdf_read_write


def main_function(staging_directory, out_dir):
	[file_names, width_of_stencil, ANC_threshold] = configure(staging_directory, out_dir);
	[xdata, ydata, data_all, dates, date_pairs] = inputs(file_names);
	[aps_array, corrected_intfs] = compute(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, ANC_threshold, out_dir);
	outputs(data_all, corrected_intfs, aps_array, dates, out_dir);
	return;

# ------------- CONFIGURE ------------ # 
def configure(staging_directory, out_dir):
	width_of_stencil = 25;   # in days.  
	ANC_threshold = 3;

	# Setting up the input and output directories. 
	file_names=glob.glob(staging_directory+"/*_*_unwrap.grd");
	file_names=sorted(file_names);
	print("Performing APS on files in %s " % staging_directory);
	if len(file_names)==0:
		print("Error! No files matching search pattern within "+staging_directory); sys.exit(1);
	call(['mkdir','-p',out_dir],shell=False);
	return [file_names, width_of_stencil, ANC_threshold];


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
	return [xdata, ydata, data_all, dates, date_pairs];



# ----------- COMPUTE ------------- #

def form_APS_pairs(date_pairs, mydate, containing, width_of_stencil):
	# Date pairs: a global list of available interferograms (in order).
	# containing: a small list of interferograms that contain the image of interest
	# width_of_stencil: the number of days we tolerate in stencils

	pairlist = [];

	# print("images that contain "+mydate+":")
	# print(containing);
	# Debugging print statements. 

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

	# print(pairlist);
	# Example of pairlist: [[63,65],[62,66]], where the numbers are the indeces of paired interferograms used for APS construction

	return pairlist; 


def compute_aps_image(dates, date_pairs, data_all, width_of_stencil, mydate):

	containing = [item for item in date_pairs if mydate in item];  # which interferograms contain the image of interest? 
	zdim, rowdim, coldim = np.shape(data_all);

	my_aps=np.zeros((rowdim,coldim));
	myaps_1darray=[];

	# A list of valid interferogram pairs by index in the stack, which we use for making APS
	pairlist = form_APS_pairs(date_pairs, mydate, containing, width_of_stencil);
	if len(pairlist)==0:
		print("ERROR! No valid APS stencils were detected in the stack for image %s" % mydate)
		return [my_aps,0];

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
				myaps_1darray.append(my_aps[i][j]);
			else:
				my_aps[i][j]=np.nan;

	print(mydate);
	ANC=compute_ANC(myaps_1darray);
	print("ANC: %.2f" % ANC);

	return [my_aps, ANC];


def compute_ANC(myaps):
	# myaps = 1d array of numeric values
	# A scaled RMS of the atmospheric phase screen.
	normalizer=10/5.0;
	M=len(myaps);
	abar = np.mean(myaps);
	res_sq_array = [(a-abar)*(a-abar) for a in myaps];
	ANC = np.sqrt((1.0/M) * np.sum(res_sq_array) );
	return ANC;


def correct_intf_stack(data_all, my_aps, given_date, date_pairs):
	# In this function, you start with a set of interferograms (3D array), and correct it for a given APS. 
	# You return the updated interferograms. 

	print("Correcting intf stack for APS on day %s" % given_date);
	zdim, rowdim, coldim = np.shape(data_all);
	corrected_intfs = data_all;

	for i in range(zdim):
		# For each interferogram, you might want to correct the image for the new APS. 
		if given_date in date_pairs[i]:
			image_split = date_pairs[i].split('_');
			if given_date==image_split[1]:  # if the given date is the second part of the image. 
				# These are the interferograms made where the APS is for the later date. The APS is subtracted. 
				for j in range(rowdim):
					for k in range(coldim): # each pixel takes a correction from the APS
						corrected_intfs[i][j][k] = (data_all[i][j][k]) - (my_aps[j][k]); 
			else:  # If the given date is the first part of the image
				for j in range(rowdim):
					for k in range(coldim): # each pixel takes a correction from the APS
						corrected_intfs[i][j][k] = (data_all[i][j][k]) + (my_aps[j][k]); 

	return corrected_intfs;



def get_initial_ANC_ranking(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, out_dir):
	# If you've done initial ranking before, you can read the results from a file (saves time)
	# If you're starting fresh, then you want to compute the ANCs for each APS as they are first calculated. 
	# There is no over-writing of data. 
	ANC_array=[];
	aps_array_initial = np.zeros((len(dates),len(ydata),len(xdata)));
	ANCfile = out_dir+'/initial_ANC.txt';
	if os.path.isfile(ANCfile):  # read the ANCs if the file exists. 
		dates=[];
		ifile=open(ANCfile,'r');
		for line in ifile:
			dates.append(line.split(': ')[0]);
			ANC_array.append(float(line.split(': ')[1]));
		ifile.close();
	else:  # Otherwise we need to do a computation. 
		ofile = open(ANCfile,'w');
		for i in range(len(dates)):
			[my_aps, ANC] = compute_aps_image(dates, date_pairs, data_all, aps_array_initial, width_of_stencil, dates[i]);
			ANC_array.append(ANC);
			ANCfile.write('%s: %.2f \n' % (dates[i], ANC) );
		ANCfile.close();
	return ANC_array; 



def compute(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, ANC_threshold, out_dir):

	# Automated selection of dates and interferograms for APS construction. This will: 
	# Rank days based on ANC
	# Make APS iteratively, in order of ANC
	# Produce an array of corrected interferogram data

	aps_array_iterated = np.zeros((len(dates),len(ydata),len(xdata)));

	# ITERATION 1: EITHER COMPUTE ANCs INITIALLY OR READ THEM FROM A FILE. 
	ANC_array = get_initial_ANC_ranking(xdata, ydata, data_all, dates, date_pairs, width_of_stencil, out_dir);

	print("Before iterations, max ANC is: %.2f " % np.max(ANC_array));
	max_index=ANC_array.index(np.max(ANC_array));
	print("Occurs on day: %s" % dates[max_index]);


	# ITERATION 2: SOLVE FOR ATMOSPHERIC SCREEN IN ORDER OF ANC
	# re-order the dates and ANCs
	corrected_intfs = deepcopy(data_all);  # Start iteratively adjusting the interferograms. 

	for i in range(len(dates)):  # this is how we control the number of iterations. 

		ordered_ANCs = [x for x,_ in sorted(zip(ANC_array, dates),reverse=True)]
		ordered_dates = [x for _,x in sorted(zip(ANC_array, dates), reverse=True)];
		# print(ordered_ANCs);
		# print(ordered_dates);
		if ordered_ANCs[0]<ANC_threshold:
			print("Exiting loop. %d images iteratively corrected. ")
			break;

		given_date=ordered_dates[0];  # THIS IS THE KEY
		# given_date='2017298';
		[my_aps,ANC]=compute_aps_image(dates, date_pairs, corrected_intfs, width_of_stencil, given_date);

		# Fix interferograms
		corrected_intfs = correct_intf_stack(corrected_intfs, my_aps, given_date, date_pairs);

		# At this stage, the ANC for your chosen image is zero by construction. 
		# At the same time, the ANC for the neighboring images have been changed. 
		# I will compute them, record them, and properly iterate later. 

		# Compute my_aps,ANC again based on fixed interferograms 
		[my_aps,new_ANC]=compute_aps_image(dates, date_pairs, corrected_intfs, width_of_stencil, given_date);

		# Replace the ANC measure with the updated one. 
		print("Iterating over %s: ANC %.2f -> ANC %.2f " % (given_date, ANC, new_ANC));
		i = dates.index(given_date);
		ANC_array[i]=new_ANC;

	return [aps_array_iterated, corrected_intfs];


# ----------- OUTPUTS ------------- # 
def view_one_image(data_all, date_pairs, mydate, my_aps, loopcheck, intf1_index, intf1_corrected, intf2_index, intf2_corrected, ANC, out_dir):

	vmin=-20;
	vmax=20;

	plt.figure();
	plt.subplot(2,3,1);
	plt.imshow(data_all[intf1_index],cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title(date_pairs[intf1_index]);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,2)
	plt.imshow(data_all[intf2_index],cmap='hsv', vmin=vmin, vmax=vmax);
	plt.title(date_pairs[intf2_index]);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,3);
	plt.imshow(my_aps, cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title('ANC: %.2f' % ANC)
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,4);
	plt.imshow(intf1_corrected,cmap='hsv',vmin=vmin, vmax=vmax);
	plt.title('intf1_corr');
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,5)
	plt.imshow(intf2_corrected,cmap='hsv', vmin=vmin, vmax=vmax);
	plt.title('intf2_corr');
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])

	plt.subplot(2,3,6)
	plt.imshow(loopcheck,cmap='jet',vmin=-0.1, vmax=2.1);
	plt.title('loops for '+mydate);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([])
	plt.savefig(out_dir+'/initialmask_'+mydate+'.png');
	plt.close();
	return;


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


def outputs(data_all, corrected_intfs, myaps_array, dates, out_dir):

	# Display raw intfs
	display_many_intfs(data_all, out_dir,'raw_intfs');
	display_many_intfs(corrected_intfs, out_dir,'aps_intfs');
	display_myaps_array(myaps_array, dates, out_dir, 'APS_arrays');

	# view_one_image(data_all, date_pairs, mydate, my_aps, loopcheck, intf1_index, intf1_corrected, intf2_index, intf2_corrected, ANC, out_dir);

	return;
