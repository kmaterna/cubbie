# Reviewing the NSBAS core functions
# July 2019

import numpy as np 
import matplotlib.pyplot as plt 
import glob
import scipy.io.netcdf as netcdf
import readmytupledata as rmd
import netcdf_read_write as rwr 


def configure():
	myfiles="glob.glob(intf_all/*/unwrap_ref.grd)";
	manual_exclude="Metadata/manual_exclude.txt";
	wavelength=56; # mm 
	nsbas_good_num=65; # % of images above coherence threshold
	smoothing=1.0;
	outfile="velocity_NSBAS_reasonable.grd"
	signal_spread_file="signalspread.nc"
	return myfiles, manual_exclude, signal_spread_file, wavelength, nsbas_good_num, smoothing, outfile;


# ------------ INPUTS ------------ #
def inputs(myfiles, signal_spread_file, manual_exclude):
	signal_spread_data=somehow_read(signal_spread_file);  # HELP! Somehow read this into an array. 
	myfiles=somehow_exclude_manual_remove_here(myfiles);  # HELP! SOMEHOW EXCLUDE SOME FILES. NOT SURE YOUR INPUT FORMAT. 

	datatuple=rmd.reader(myfiles);
	dates, date_pairs = read_dates_and_datepairs(myfiles);
	print("Reading %d interferograms from %d acquisitions. " % (len(date_pairs), len(dates) ) );
	return datatuple, signal_spread_data, dates, date_pairs;


def read_dates_and_datepairs(myfiles):
	date_pairs=[];
	dates=[];
	for ifile in myfiles:  
		pairname=ifile.split('/')[-2][0:15];
		image1=pairname.split('_')[0];
		image2=pairname.split('_')[1];
		date_pairs.append(pairname);  # returning something like '2016292_2016316' for each intf
		dates.append(image1);
		dates.append(image2);
	dates=list(set(dates)); # find the unique days. 
	dates=sorted(dates);
	print(dates);
	return dates, datepairs;


# ------------ COMPUTE ------------ #
def compute(datatuple, nsbas_good_num, signal_spread_data, dates, date_pairs, smoothing, wavelength, outfile):
	zdata=datatuple.zvalues;
	# The point here is to loop through each pixel, determine if there's enough data to use, and then 
	# make an SBAS matrix describing each image that's a real number (not nan). 
	print("Analyzing the nsbas timeseries per pixel.")
	ofile=open(outfile,'w');
	[zdim, xdim, ydim] = np.shape(zdata)
	vel = np.zeros([xdim, ydim]);
	
	for i in range(xdim):  # A loop through each pixel. 
		for j in range(ydim):

			pixel_value = [zdata[k][i][j] for k in range(zdim)];  # slicing the values of phase for a pixel across the various interferograms
			
			if signal_spread_data[i][j] >= nsbas_good_num:  # If we have a pixel that will be analyzed: Do SBAS
				print("%d %d " % (i, j) )
				vel[i][j] = do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength); 
				ofile.write("%d %d %f\n" % (i, j, vel[i][j]) );
				# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
				# dates: if we have 35 images, this is the date of each image
				# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image
				# This solves Gm = d for the movement of the pixel with smoothing. 
			else:
				vel[i][j]=np.nan;

	ofile.close();

	return vel;



def do_nsbas_pixel(pixel_value, dates, date_pairs, smoothing, wavelength):
	# pixel_value: if we have 62 intf, this is a (62,) array of the phase values in each interferogram. 
	# dates: if we have 35 images, this is the date of each image
	# date_pairs: if we have 62 intf, this is a (62) list with the image pairs used in each image	
	# for x in range(len(dates)-1):

	d = np.array([]);
	dates=sorted(dates);
	date_pairs_used=[];
	for i in range(len(pixel_value)):
		if not math.isnan(pixel_value[i]):
			d = np.append(d, pixel_value[i]);  # removes the nans from the computation. 
			date_pairs_used.append(date_pairs[i]);  # might be a slightly shorter array of which interferograms actually got used. 
	model_num=len(dates)-1;

	G = np.zeros([len(date_pairs_used)+model_num-1, model_num]);  # in one case, 91x35
	# print(np.shape(G));
	
	for i in range(len(d)):  # building G matrix line by line. 
		ith_intf = date_pairs_used[i];
		first_image=ith_intf.split('_')[0]; # in format '2017082'
		second_image=ith_intf.split('_')[1]; # in format '2017094'
		first_index=dates.index(first_image);
		second_index=dates.index(second_image);
		for j in range(second_index-first_index):
			G[i][first_index+j]=1;

	# Building the smoothing matrix with 1, -1 pairs
	for i in range(len(date_pairs_used),len(date_pairs_used)+model_num-1):
		position=i-len(date_pairs_used);
		G[i][position]=1*smoothing;
		G[i][position+1]=-1*smoothing;
		d = np.append(d,0);

	# solving the SBAS linear least squares equation for displacement between each epoch. 
	m = np.linalg.lstsq(G,d)[0];  

	# modeled_data=np.dot(G,m);
	# plt.figure();
	# plt.plot(d,'.b');
	# plt.plot(modeled_data,'.--g');
	# plt.savefig('d_vs_m.eps')
	# plt.close();

	# Adding up all the displacement. 
	m_cumulative=[];
	m_cumulative.append(0);
	for i in range(len(m)):
		m_cumulative.append(np.sum(m[0:i]));  # The cumulative phase from start to finish! 


	# Solving for linear velocity
	x_axis_datetimes=[dt.datetime.strptime(x,"%Y%j") for x in dates];
	x_axis_days=[(x - x_axis_datetimes[0]).days for x in x_axis_datetimes];  # number of days since first acquisition. 

	x=np.zeros([len(x_axis_days),2]);
	y=np.array([]);
	for i in range(len(x_axis_days)):
		x[i][0]=x_axis_days[i];
		x[i][1]=1;  
		y=np.append(y,[m_cumulative[i]]);
	model_slopes = np.linalg.lstsq(x,y)[0];  # units: phase per day. 
	model_line = [model_slopes[1]+ x*model_slopes[0] for x in x_axis_days];

	# Velocity conversion: units in mm / year
	vel=model_slopes[0];  # in radians per day
	vel=vel*wavelength*365.24/2.0/(2*np.pi);

	return vel;





if __name__=="__main__":
	myfiles, manual_exclude, signal_spread_file, wavelength, nsbas_good_num, smoothing, outfile = configure();
	datatuple, signal_spread_data, dates, date_pairs = inputs(myfiles, signal_spread_file, manual_exclude);
	vel = compute(datatuple, nsbas_good_num, signal_spread_data, dates, date_pairs, smoothing, wavelength, outfile_writer);
	rwr.produce_output_netcdf(datatuple.xvales, datatuple.yvales, vel, 'velocity','vel.grd');
	rwr.flip_if_necessary('vel.grd');
