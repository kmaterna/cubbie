# A script to take a stack of images plus a DEM
# Solve for best-fitting linear trend
# (Globally or locally)
# Remove trend and save the adjusted stack in out_dir

import glob, sys, os
import subprocess
import netcdf_read_write


def main_function(staging_directory, outdir):
	[filenames,demfile]=configure(staging_directory, outdir);
	demdata = subsample_read_dem(filenames[0], demfile);

	for item in filenames:
		[xdata, ydata, zdata]=netcdf_read_write.read_grd_xyz(item);
		corrected_zdata= global_compute_item(zdata,demdata);
		#write_item(xdata, ydata, corrected_zdata, item, outdir);

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


# -------- PREPROCESS and I/O FUNCTIONS ------ # 
def subsample_read_dem(samplefile, demfile):
	# Take an example interferogram and subsample the topo_ra to exactly match this file size. 
	# You may have to play with -T/-r options in GMT grdsample to force the same gridcell/gridline registration
	# You also may have to force the netcdf4 file into netcdf3 for later reading in python. 

	if not os.path.isfile(demfile):
		print("ERROR! %s does not exist- exiting " % demfile);
		sys.exit(1);

	subsampled_file='topo/topo_ra_subsampled.grd';
	intervals=subprocess.check_output('gmt grdinfo -I '+samplefile,shell=True);
	intervals=intervals.split('\n')[0];
	ranges=subprocess.check_output('gmt grdinfo -I- '+samplefile,shell=True);
	ranges=ranges.split('\n')[0];
	command='gmt grdsample '+demfile+' -Gtopo/temp.grd -T '+intervals+' '+ranges;
	print(command);
	subprocess.call(command, shell=True);
	subprocess.call('nccopy -k classic topo/temp.grd '+subsampled_file,shell=True); 
	subprocess.call(['rm','topo/temp.grd'],shell=False);

	[x,y,z]=netcdf_read_write.read_grd_xyz(subsampled_file);
	return z;



# --------- BASIC I/O FUNCTIONS ------------ # 
def read_dem(demfile):



	return 0;

def write_item(xdata, ydata, corrected_zdata, item, outdir):
	item_name=item.split('/')[-1];
	netcdfname=outdir+'/'+item_name;
	netcdf_read_write.produce_output_netcdf(xdata, ydata, corrected_zdata, 'unwrapped_phase', netcdfname);
	netcdf_read_write.flip_if_necessary(netcdfname);
	return;



# ---------- COMPUTE FUNCTIONS ------------ # 

def global_compute_item(zdata, demdata):
	return 0;






if __name__=="__main__":
	main_function('intf_all/unwrap.grd', 'intf_all/atm_topo_corrected.grd');
