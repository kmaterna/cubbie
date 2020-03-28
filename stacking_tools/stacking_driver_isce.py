
"""
	specific driver for stacking, velocities, and time series
	Only call this if you have isce installed on your machine. 
"""

import stacking_configparser
import stacking_main_functions
import stacking_functions_isce

if __name__=="__main__":
	config_params = stacking_configparser.read_config_isce();

	# Step 0: set up output directories
	stacking_main_functions.set_up_output_directories(config_params); 

	# Step 1: make atmospheric corrections, etc. prior to TS
	stacking_functions_isce.make_corrections_isce(config_params); 

	# Step 2: make reference unwrapped
	stacking_functions_isce.collect_unwrap_ref(config_params); 

	# Step 3: make velocity field
	stacking_functions_isce.vels_and_ts(config_params); 

	# Step 4: geocoding
	stacking_functions_isce.geocode_vels(config_params);