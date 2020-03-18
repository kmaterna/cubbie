import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
import matplotlib.pyplot as plt 
from subprocess import call, check_output
import stacking_utilities
import readmytupledata
import stack_corr
import netcdf_read_write as rwr


# --------------- UTILITIES ------------ # 
def get_100p_pixels(filenameslist):
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

	for i in range(np.shape(mydata.zvalues)[1]):
		for j in range(np.shape(mydata.zvalues)[2]):
			oneslice = mydata.zvalues[:,i,j];
			if np.sum(~np.isnan(oneslice))==total_images:  # if we have perfect coherence
				count=count+1;

	a=stack_corr.stack_corr(mydata, np.nan);
	rwr.produce_output_netcdf(xvalues, yvalues, a, 'Percentage', 'signalspread_cut.nc');
	rwr.produce_output_plot('signalspread_cut.nc', 'Signal Spread', 'signalspread_cut.png', 'Percentage of coherence', aspect=1/4, invert_yaxis=False )

	print("%d of %d (%f percent) are totally coherent. " % (count, total_pixels, 100*(count/total_pixels)) );
	print(np.shape(mydata.zvalues));
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

    get_100p_pixels(intf_files);
    sys.exit(0);




    # Here we need to get ref_idx if we don't have it already
    # rowref, colref = stacking_utilities.get_ref_index(config_params.ref_swath, config_params.swath, config_params.ref_loc, config_params.ref_idx, intf_files);

    # Now we coalesce the files and reference them to the right value/pixel
    # stacking_utilities.make_referenced_unwrapped(intf_files, config_params.swath, config_params.ref_swath, rowref, colref, config_params.ref_dir);

    return;