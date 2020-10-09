#!/usr/bin/env python
"""
	specific driver for stacking, velocities, and time series
"""

import stacking_configparser
import stacking_functions_gmstar

if __name__ == "__main__":
	config_params = stacking_configparser.read_config();

	# Step 0: set up output directories
	stacking_functions_gmstar.set_up_output_directories(config_params);

	# Step 1: make atmospheric corrections, etc.
	stacking_functions_gmstar.make_corrections(config_params);

	# Step 2: Get reference information
	stacking_functions_gmstar.get_ref(config_params);

	# Step 3: make velocity field
	stacking_functions_gmstar.vels_and_ts(config_params);

	# Step 4: geocoding
	stacking_functions_gmstar.geocode_vels(config_params);
