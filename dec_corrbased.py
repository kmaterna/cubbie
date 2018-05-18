import numpy as np 
import matplotlib.pyplot as plt 
import scipy.io.netcdf as netcdf
import glob, sys
import os.path
from subprocess import call, check_output
import nsbas

# The purpose of this function is to decimate a full-resolution image. 
# We decimate to reduce file size, speed up unwrapping, and maximize coherence. 
# In this script, we take the one pixel with highest coherence in each box. 
# This script does some renaming: 
# phase.grd --> phase_full.grd
# phasefilt.grd --> phasefilt_full.grd
# corr.grd --> corr.grd
# mask.grd --> mask_full.grd
# decimated output --> phasefilt.grd
# decimated corr --> corr.grd
# decimated mask --> mask.grd

def decimate_main_function(xdec, ydec):
	print("Decimating phasefilt.grd files based on maximum correlation pixels. ")
	[corrfiles] = configure_manyfiles();
	for i in range(len(corrfiles)):
		print(corrfiles[i]);
		perform_decimation_one_file(xdec, ydec, corrfiles[i]);
	return;


# This is the main function for a single file (preserving the functionality of working with a single file, for later use)
def perform_decimation_one_file(xdec, ydec, corrfile):
	print("Decimating phasefilt.grd based on "+corrfile);
	[phase_infile, outfile, corr_full, corroutfile, maskfile_full, maskoutfile, computeflag]=configure_singlefile(corrfile);
	if computeflag==0:
		return;
	else:
		[xc, yc, zc, xp, yp, zp] = inputs(corr_full, phase_infile);
		[xstar, ystar, phasestar, corrstar] = process_by_correlation(xc, yc, zc, xp, yp, zp, xdec, ydec);
		[xm, ym, zm] = read_grd(maskfile_full);
		#[xstar, ystar, maskstar]=compute_maskfile(xm, ym, zm, xdec, ydec);
		#produce_output_netcdf(xstar, ystar, maskstar, 'mask',maskoutfile);  # commented for specific experiment
		produce_output_netcdf(xstar, ystar, phasestar,'radians', outfile);
		produce_output_netcdf(xstar, ystar, corrstar, 'correlation', corroutfile);
	return;

# This is assuming we call from the processing directory
def configure_manyfiles():
	corrs=glob.glob("intf_all/20*/corr.grd");  # corr.grd files with their associated directories
	return [corrs];


def configure_singlefile(corrfile):
	folder_components=corrfile.split('/');  # a list of strings
	folder=''
	for i in range(len(folder_components)-1):
		folder=folder+folder_components[i];
		folder=folder+'/'  # something like "intf_all/2015177_2016010/"
	phasefilt=folder + 'phasefilt.grd';  # something like 'phasefilt'
	phasenofilt=folder + 'phase.grd';  # something like 'phase'
	phasefilt_full=folder+'phasefilt_full.grd';
	phasenofilt_full=folder+'phase_full.grd';
	corr_full=folder+'corr_full.grd';
	maskfile=folder+'mask.grd';
	maskfile_full=folder+'mask_full.grd';
	control_file=folder+'amp.grd';  # This is a full-size file that NEVER CHANGES. 


	# Check if we should copy the phase over into phase_full.grd based on file size. 
	phasesize = check_output("ls -lh "+phasenofilt+" | awk \'{print $5}\'", shell=True);
	corrsize  = check_output("ls -lh "+corrfile+" | awk \'{print $5}\'", shell=True);
	masksize  = check_output("ls -lh "+maskfile+" | awk \'{print $5}\'", shell=True);
	ampsize  = check_output("ls -lh "+control_file+" | awk \'{print $5}\'", shell=True);
	if phasesize==ampsize:  
	# the phase.grd is not decimated or the decimated file has been deleted. We should proceed to decimate. 
		# Move phasefilt.grd --> phasefilt_full.grd; phase.grd --> phase_full.grd. 
		# This will not be satisfied for mid-stream starts, since phase.grd has already been moved to phase_full.grd
		print("Moving phase.grd into phase_full.grd before decimating");
		if phasesize==ampsize:
			call("mv "+phasefilt+" "+phasefilt_full, shell=True);  # move phasefilt.grd into phasefilt_full.grd
			call("mv "+phasenofilt+" "+phasenofilt_full, shell=True);  # move phase.grd into phase_full.grd
		if corrsize == ampsize:
			call("mv "+corrfile+" "+corr_full, shell=True);  # move phase.grd into phase_full.grd
		if masksize == ampsize:
			call("mv "+maskfile+" "+maskfile_full, shell=True);  # move phase.grd into phase_full.grd
		computeflag=1;  # do the user-defined decimation based on correlation. 
	else:
		print("phase.grd and amp.grd are not the same size. No need to decimate. Skipping. ")
		computeflag=1;  # do not do anything; the decimating has already been done. 
		# Set to 1 for specific experiment. 

	dec_phasefile= folder+'phasefilt.grd';  # the output will be phasefilt.grd (so that SNAPHU will unwrap it)
	corroutfile=folder+'corr.grd';
	maskoutfile=folder+'mask.grd';
	phase_infile=folder+'phasefilt_full.grd';  # what file are we actually reading and decimating?  
	return [phasenofilt_full, dec_phasefile, corr_full, corroutfile, maskfile_full, maskoutfile, computeflag];



# --------------- INPUTS --------------- # 
def inputs(corrfile, phasefile):
	[xc,yc,zc] = read_grd(corrfile);
	[xp,yp,zp] = read_grd(phasefile);
	return [xc, yc, zc, xp, yp, zp];

def read_grd(filename):
	xdata0 = netcdf.netcdf_file(filename,'r').variables['x'][::-1];
	ydata0 = netcdf.netcdf_file(filename,'r').variables['y'][::-1];
	zdata0 = netcdf.netcdf_file(filename,'r').variables['z'][::-1];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	zdata=zdata0.copy();
	return [xdata, ydata, zdata]; 


# --------------- COMPUTE --------------- # 
def process_by_correlation(xc, yc, zc, xp, yp, zp, xdec, ydec):
	xdec=int(xdec);
	ydec=int(ydec);
	number_of_xboxes=int(np.ceil(len(xc)/xdec));
	number_of_yboxes=int(np.ceil(len(yc)/ydec));
	newx    =np.zeros([number_of_xboxes]);
	newy    =np.zeros([number_of_yboxes]);
	newphase=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array
	newcorr=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array
	#print(np.shape(newphase));

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
			myzp = zp[min_y:max_y,min_x:max_x];  # phase.     
			# Numpy note: zp[0:2][0:2] is not the same thing as zp[0:2,0:2];
			
			# The math itself: find the maximum correlation and take the corresponding phase.
			try:
				mymax=np.nanargmax(myzc);
			except ValueError:
				newphase[i][j]=np.nan;
				newcorr[i][j]=0.0;
			else:
				indices=np.unravel_index(mymax,myzc.shape)
				newphase[i][j]=myzp[indices[0]][indices[1]];  # the indices of the maximum correlation.
				newcorr[i][j] =myzc[indices[0]][indices[1]];
				# newphase[i][j]=myzp[0][0]; In case we want to see random decimation

	newphase=newphase.T;
	newcorr=newcorr.T;

	return [newx, newy, newphase, newcorr];


def compute_maskfile(xm, ym, zm, xdec, ydec):
	xdec=int(xdec);
	ydec=int(ydec);
	number_of_xboxes=int(np.ceil(len(xm)/xdec));
	number_of_yboxes=int(np.ceil(len(ym)/ydec));
	newx    =np.zeros([number_of_xboxes]);
	newy    =np.zeros([number_of_yboxes]);
	newmask=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array

	for i in range(number_of_xboxes):
		for j in range(number_of_yboxes):

			# Go through each box of (xdec, ydec) pixels and find the max correlated pixel. 
			min_x=i*xdec;  # start counting at 0
			max_x=min(i*xdec+xdec,len(xm));  # end counting at 0+xdec, or the end of the array
			min_y=j*ydec;
			max_y=min(j*ydec+ydec,len(ym));
			newx[i]=np.median(xm[min_x:max_x]);
			newy[j]=np.median(ym[min_y:max_y]);			

			myzm = zm[min_y:max_y,min_x:max_x];  # correlation
			
			# The math itself: find the maximum correlation and take the corresponding phase.
			try:
				mymax=np.nanargmax(myzm);
			except ValueError:
				newmask[i][j]=np.nan;
			else:
				indices=np.unravel_index(mymax,myzm.shape)
				newmask[i][j] =myzm[indices[0]][indices[1]];
	newmask=newmask.T;

	return [newx, newy, newmask];


# -------------- OUTPUTS ------------- # 
def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname):
	# # Write the netcdf velocity grid file.  
	f=netcdf.netcdf_file(netcdfname,'w');
	f.history = 'Decimated by correlation method';
	f.createDimension('x',len(xdata));
	f.createDimension('y',len(ydata));
	print(np.shape(zdata));
	x=f.createVariable('x',float,('x',))
	x[:]=xdata;
	x.units = 'range';
	y=f.createVariable('y',float,('y',))
	y[:]=ydata;
	y.units = 'azimuth';
	z=f.createVariable('z',float,('y','x',));
	z[:,:]=zdata;
	z.units = zunits;
	f.close();

	nsbas.flip_if_necessary(netcdfname);
	return;


if __name__=="__main__":
	decimate_main_function(xdec, ydec);
