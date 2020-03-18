import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
import matplotlib.pyplot as plt 
from subprocess import call, check_output
import stacking_utilities
import readmytupledata


# --------------- UTILITIES ------------ # 
def get_100p_pixels(filenameslist):
	print("Finding the pixels that are 100 percent coherent");
	# Finding pixels that are completely non-nan
	mydata = readmytupledata.reader_isce(filenameslist,band=2);
	total_pixels = np.shape(mydata.zvalues)[1]*np.shape(mydata.zvalues)[2];
	total_images = np.shape(mydata.zvalues)[0];
	count = 0;

	# Looks like we wrote the interpolated values, not masked interpolated values.
	# Will have to re-run the code tomorrow and see what actually happens.  
	fig = plt.figure();
	plt.imshow(mydata.zvalues[10,:,:],aspect=1/5);
	plt.savefig('example_unwrap_manually_masked.png');


	for i in range(np.shape(mydata.zvalues)[1]):
		for j in range(np.shape(mydata.zvalues)[2]):
			oneslice = mydata.zvalues[:,i,j];
			# print(oneslice);
			# print(np.sum(np.isnan(oneslice)));
			if np.sum(~np.isnan(oneslice))==total_images:  # if we have perfect coherence
				count=count+1;
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