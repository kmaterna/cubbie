# A script to take a stack of images plus a DEM
# Solve for best-fitting linear trend
# (Globally or locally)
# Remove trend and save the adjusted stack in out_dir
# In order to work correctly, this script needs a specified reference pixel. 

import numpy as np
import matplotlib.pyplot as plt 
import glob, sys, os
import subprocess
import netcdf_read_write


def main_function(staging_directory, outdir, rowref, colref):
	[filenames,demfile]=configure(staging_directory, outdir);
	demdata = subsample_read_dem(filenames[0], demfile);

	for item in filenames:
		[xdata, ydata, zdata]=netcdf_read_write.read_grd_xyz(item);
		[corrected_zdata, zarray, corrarray, demarray]= global_compute_item(zdata,demdata, rowref, colref);
		output_item(xdata, ydata, corrected_zdata, zarray, corrarray, demarray, item, outdir);
	return;




# ------ CONFIGURE THE MAJOR LOOP ------------ # 
def configure(staging_directory, outdir):
	print("Detrending atm/topo for all files in %s and storing the result in %s " % (staging_directory, outdir) )
	file_names=glob.glob(staging_directory+"/*_*_unwrap.grd");
	if len(file_names)==0:
		print("Error! No files matching search pattern within "+staging_directory); sys.exit(1);
	subprocess.call(['mkdir','-p',outdir],shell=False);
	demfile='topo/topo_ra.grd';
	return [file_names, demfile];


# -------- PREPROCESS I/O FUNCTIONS ---------- # 
def subsample_read_dem(samplefile, demfile):
	# Take an example interferogram and subsample the topo_ra to exactly match this file size. 
	# You may have to play with -T/-r options in GMT grdsample to force the same gridcell/gridline registration
	# You also may have to force the netcdf4 file into netcdf3 for later reading in python. 

	if not os.path.isfile(demfile):
		print("ERROR! %s does not exist- exiting " % demfile);
		sys.exit(1);

	subsampled_file='topo/topo_ra_subsampled.grd';
	intervals=subprocess.check_output(['gmt','grdinfo','-I',samplefile],shell=False);
	intervals=intervals.split('\n')[0];
	ranges=subprocess.check_output(['gmt','grdinfo','-I-',samplefile],shell=False);
	ranges=ranges.split('\n')[0];
	command='gmt grdsample '+demfile+' -Gtopo/temp.grd -T '+intervals+' '+ranges;
	print(command);
	subprocess.call(command, shell=True);
	subprocess.call('nccopy -k classic topo/temp.grd '+subsampled_file,shell=True); 
	subprocess.call(['rm','topo/temp.grd'],shell=False);

	[x,y,z]=netcdf_read_write.read_grd_xyz(subsampled_file);
	return z;




# ---------- COMPUTE FUNCTIONS ------------ # 

def global_compute_item(zdata, demdata, rowref, colref):

	zarray=[];  # the 1D array with original phase values
	demarray=[];  # the 1D array with topography
	corrarray=[];  # the 1D array with corrected phase

	# Collect valid phase values
	rowdim,coldim=np.shape(zdata);
	for i in range(rowdim):
		for j in range(coldim):
			if ~np.isnan(zdata[i][j]):
				zarray.append(zdata[i][j]);
				demarray.append(demdata[i][j]);

	# Now generate a best-fitting slope between phase and topography
	coef = np.polyfit(demarray, zarray, 1);

	# Now remove the slope (pinning everything so that the reference pixel has phase 0)
	corrected_zdata=np.zeros(np.shape(zdata));
	rowdim,coldim=np.shape(zdata);
	reference_pixel_offset = zdata[rowref][colref]-coef[0]*demdata[rowref][colref];
	for i in range(rowdim):
		for j in range(coldim):
			if ~np.isnan(zdata[i][j]):
				corrected_zdata[i][j]=zdata[i][j]-coef[0]*demdata[i][j] - reference_pixel_offset;
				corrarray.append(zdata[i][j]-coef[0]*demdata[i][j] - reference_pixel_offset);
			else:
				corrected_zdata[i][j]=np.nan;

	#print(corrected_zdata[rowref][colref]);

	return [corrected_zdata, zarray, corrarray, demarray];



# --------------- OUTPUTS ---------------- # 
def output_item(xdata, ydata, corrected_zdata, zarray, corrarray, demarray, item, outdir):
	item_name=item.split('/')[-1];
	item_name_short=item_name.split('unwrap.grd')[0];
	netcdfname=outdir+'/'+item_name;
	netcdf_read_write.produce_output_netcdf(xdata, ydata, corrected_zdata, 'unwrapped_phase', netcdfname);
	netcdf_read_write.flip_if_necessary(netcdfname);


	plt.figure();
	plt.plot(demarray,zarray,'.');
	plt.plot(demarray,corrarray,'.r',alpha=0.15)
	plt.ylabel('phase');
	plt.xlabel('topo');
	plt.title('Initial and Corrected Phase vs. Topography')
	plt.savefig('intf_all/atm_topo_corrected.grd'+'/phase_topo_'+item_name_short+'.png');
	plt.close();

	return;


if __name__=="__main__":
	main_function('intf_all/unwrap.grd', 'intf_all/atm_topo_corrected.grd', 241, 175);
