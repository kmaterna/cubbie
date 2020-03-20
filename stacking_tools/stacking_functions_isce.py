import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
import matplotlib.pyplot as plt 
from subprocess import call, check_output
import stacking_utilities
import readmytupledata
import stack_corr
import netcdf_read_write as rwr
import isce_read_write


# --------------- UTILITIES ------------ # 
def get_100p_pixels_get_ref(filenameslist, ref_idx):
	# If we have a reference idx in the file, just return that. 
	# If not, then you might have to run through this function a number of times 
	# To select your boxes and your eventual reference pixel. 
	# I pick one in a stable area, outside of the deformation, ideally in a desert. 
	if ref_idx != "":
		print("Returning reference row,col from the config file: %s " % ref_idx);
		return int(ref_idx.split(',')[0]), int(ref_idx.split(',')[1])
	else:
		# We don't have a reference index chosen, so we go looking. 
		# This step needs some manual adjustment. 

		print("Finding the pixels that are 100 percent coherent");
		# Finding pixels that are completely non-nan
		mydata = readmytupledata.reader_isce(filenameslist,band=1); 
		# if it's straight from SNAPHU, band = 2
		# if it's manually processed first, band = 1
		total_pixels = np.shape(mydata.zvalues)[1]*np.shape(mydata.zvalues)[2];
		total_images = np.shape(mydata.zvalues)[0];
		xvalues = range(np.shape(mydata.zvalues)[2]);
		yvalues = range(np.shape(mydata.zvalues)[1]);
		count = 0;
		ypixels_good=[];
		xpixels_good=[];
		ypixels_options=[];
		xpixels_options=[];

		for i in range(np.shape(mydata.zvalues)[1]):
			for j in range(np.shape(mydata.zvalues)[2]):
				oneslice = mydata.zvalues[:,i,j];
				if np.sum(~np.isnan(oneslice))==total_images:  # if we have perfect coherence
					count=count+1;
					xpixels_good.append(j);
					ypixels_good.append(i);
					# Here we will adjust parameters until we find a reference pixel that we like. 
					if j>248 and j<262 and i>1970 and i<2000:
						xpixels_options.append(j);
						ypixels_options.append(i);

		idx_lucky = 120;
		xref = xpixels_options[idx_lucky];
		yref = ypixels_options[idx_lucky];

		# # Make stack-corr
		# # The first time through, you might want to run stack-corr. 
		# a=stack_corr.stack_corr(mydata, np.nan);
		# rwr.produce_output_netcdf(xvalues, yvalues, a, 'Percentage', 'signalspread_cut.nc');
		# rwr.produce_output_plot('signalspread_cut.nc', 'Signal Spread', 'signalspread_cut.png', 'Percentage of coherence', aspect=1/4, invert_yaxis=False )

		print("%d of %d (%f percent) are totally coherent. " % (count, total_pixels, 100*(count/total_pixels)) );
		print(np.shape(mydata.zvalues));
		print("%d pixels are good options for the reference pixel. " % (len(xpixels_options)) );

		# Make a plot that shows where those pixels are
		a=rwr.read_grd('signalspread_cut.nc');
		fig = plt.figure();
		plt.imshow(a,aspect=1/4, cmap='rainbow');
		plt.plot(xpixels_good, ypixels_good, '.', color='k');
		plt.plot(xpixels_options,ypixels_options, '.', color='g');
		plt.plot(xref, yref, '.', color='r');
		plt.savefig('best_pixels.png');
		plt.close();

		print("Selecting reference pixel at row/col %d, %d " % (yref, xref) );

		return yref, xref;


def make_referenced_unwrapped_isce(intf_files, rowref, colref, output_dir):
	# This is simpler than the Sentinel case, since you don't have to cross 
	# swath borders or solve for n pi. 
	# We just return the values minus a single reference pixel, 
	# and write them to a directory. 
	output_dir="F/"+output_dir;
	print("Imposing reference pixel on %d files; saving output in %s" % (len(intf_files), output_dir) );
	for item in intf_files:
		datestr = item.split('/')[2];
		outname=output_dir+"/"+datestr+".refunwrapped";
		print("Making %s " % outname);
		zdata = isce_read_write.read_scalar_data(item);
		refvalue = zdata[rowref, colref];
		xdata = range(0,np.shape(zdata)[1]);
		ydata = range(0,np.shape(zdata)[0]);
		referenced_zdata = stacking_utilities.apply_reference_value(xdata, ydata, zdata, refvalue);
		referenced_zdata = np.float32(referenced_zdata);
		isce_read_write.write_isce_data(referenced_zdata, len(xdata), len(ydata), dtype="FLOAT",filename=outname);
	print("Done making reference unwrapped");
	total_intf_list=glob.glob(output_dir+"/*.refunwrapped");
	print("%s contains %d files " % (output_dir, len(total_intf_list)) );
	return;



# --------------- STEP 2: Make ref_unwrap.grd ------------ # 

def collect_unwrap_ref(config_params):
    if config_params.startstage>2:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<2:   # if we're ending at intf, we don't do this. 
        return;

    print("Stage 2 - Collecting referenced unwrapped.");

    # Very general, takes all files and doesn't discriminate. 
    intf_files=stacking_utilities.get_list_of_intf_all(config_params);

    rowref, colref = get_100p_pixels_get_ref(intf_files, config_params.ref_idx);

    # Now we coalesce the files and reference them to the right value/pixel
    make_referenced_unwrapped_isce(intf_files, rowref, colref, config_params.ref_dir);

    return;

# --------------- STEP 3: Velocities and Time Series! ------------ # 
def vels_and_ts(config_params):
	if config_params.startstage>3:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<3:   # if we're ending at intf, we don't do this. 
		return;

	# This is where the hand-picking takes place. 
	# Ex: manual excludes, manual selects, long intfs only, ramp-removed, atm-removed, etc. 
	intfs = stacking_utilities.make_selection_of_intfs(config_params);
	
	if config_params.ts_type=="STACK":
		print("Running velocities by simple stack.")
		# sss.drive_velocity_simple_stack(config_params.swath, intfs, config_params.wavelength, config_params.ts_output_dir);
	if config_params.ts_type=="SBAS":
		print("Running velocities and time series by SBAS");
		# sbas.drive_velocity_sbas(config_params.swath, intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir);
		# sbas.drive_ts_sbas(config_params);
	if config_params.ts_type=="NSBAS":
		print("Running velocities and time series by NSBAS");
		# nsbas.drive_velocity_nsbas(config_params.swath, intfs, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir);
		# nsbas.drive_ts_nsbas(config_params);

	return; 



# --------------- STEP 4: Geocoding Velocities ------------ # 
def geocode_vels(config_params):
	if config_params.startstage>4:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<4:   # if we're ending at intf, we don't do this. 
		return; 
