
"""
	specific driver for stacking, velocities, and time series
"""

import stacking_main_functions

if __name__=="__main__":
	config_params = stacking_main_functions.read_config();

	# Step 0: set up output directories
	stacking_main_functions.set_up_output_directories(config_params); 

	# Step 1: make reference unwrapped
	stacking_main_functions.collect_unwrap_ref(config_params); 

	# Step 2: make velocity field
	stacking_main_functions.velocities(config_params); 

	# Step 3: geocoding
	stacking_main_functions.geocode_vels(config_params);