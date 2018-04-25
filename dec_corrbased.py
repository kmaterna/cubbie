import numpy as np 
import matplotlib.pyplot as plt 
import scipy.io.netcdf as netcdf
import glob
from subprocess import call, check_output
import nsbas

# The purpose of this function is to decimate a full-resolution image. 
# We decimate to reduce file size, speed up unwrapping, and maximize coherence. 
# In this script, we take the one pixel with highest coherence in each box. 
# This script does some crazy renaming: 
# phase.grd --> phase_full.grd
# pphasefilt.grd --> phasefilt_full.grd
# decimated output --> phase.grd

def decimate_main_function(xdec, ydec):
	print("Decimating phasefilt.grd files based on maximum correlation pixels. ")
	[corrfiles, phasefiles, phasenofilts] = configure_manyfiles();
	for i in range(len(corrfiles)):
		perform_decimation_one_file(xdec, ydec, corrfiles[i], phasefiles[i], phasenofilts[i]);
	return;


# This is the main function for a single file (preserving the functionality of working with a single file, for later use)
def perform_decimation_one_file(xdec, ydec, corrfile, phasefile, phasenofilt):
	print("Decimating for "+phasefile);
	[phase_infile, outfile]=configure_singlefile(phasefile, phasenofilt);
	[xc, yc, zc, xp, yp, zp] = inputs(corrfile, phase_infile);
	[xstar, ystar, phasestar] = process_by_correlation(xc, yc, zc, xp, yp, zp, xdec, ydec);
	produce_output_netcdf(xstar, ystar, phasestar,'radians',outfile);
	return;

# This is assuming we call from the processing directory
def configure_manyfiles():
	phasefilts=glob.glob("intf_all/20*/phasefilt.grd");
	phasenofilts=glob.glob("intf_all/20*/phase.grd");
	corrs=glob.glob("intf_all/20*/corr.grd");
	return [corrs, phasefilts, phasenofilts];

def configure_singlefile(phasefile, phasenofilt):
	phasefilt_components=phasefile.split('/');  # a list of strings
	local_phasefilt_nogrd=phasefilt_components[-1].split('.')[0];  # something like 'phasefilt'
	local_phasenofilt_nogrd=phasenofilt.split('/')[-1].split('.')[0];  # something like 'phase'
	folder=''
	for i in range(len(phasefilt_components)-1):
		folder=folder+phasefilt_components[i];
		folder=folder+'/'  # something like "intf_all/2015177_2016010/"

	# Check if we should copy the phase over into phase_full.grd based on file size. 
	phase_infile=folder+local_phasefilt_nogrd+"_full.grd";
	phase_nofilt_full=folder+local_phasenofilt_nogrd+"_full.grd"
	phasesize = check_output("ls -lh "+folder+"phase.grd | awk \'{print $5}\'", shell=True);
	corrsize  = check_output("ls -lh "+folder+"corr.grd | awk \'{print $5}\'", shell=True);
	if phasesize==corrsize:  # the phase.grd is not decimated. We should proceed to decimate. 
		# Move phasefilt.grd --> phasefilt_full.grd; phase.grd --> phase_full.grd. 
		print("Moving phase.grd into phase_full.grd before decimating");
		call("mv "+phasefile+" "+phase_infile, shell=True);  # move phasefilt.grd into phasefilt_full.grd
		call("mv "+phasenofilt+" "+phase_nofilt_full, shell=True);  # move phase.grd into phase_full.grd

	dec_phasefile= folder+"phase.grd";  # the output will be phase.grd (so that SNAPHU will unwrap it)
	return [phase_infile, dec_phasefile];

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

def process_by_correlation(xc, yc, zc, xp, yp, zp, xdec, ydec):
	xdec=int(xdec);
	ydec=int(ydec);
	number_of_xboxes=int(np.ceil(len(xc)/xdec));
	number_of_yboxes=int(np.ceil(len(yc)/ydec));
	newx    =np.zeros([number_of_xboxes]);
	newy    =np.zeros([number_of_yboxes]);
	newphase=np.zeros([number_of_xboxes, number_of_yboxes]);  # the new array
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
			else:
				indices=np.unravel_index(mymax,myzc.shape)
				newphase[i][j]=myzp[indices[0]][indices[1]];  # the indices of the maximum correlation. 
				# newphase[i][j]=myzp[0][0]; In case we want to see random decimation

	newphase=newphase.T;

	return [newx, newy, newphase];


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
