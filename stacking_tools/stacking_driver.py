#!/usr/bin/env python
"""
	specific driver for stacking, velocities, and time series
"""

import stacking_configparser
import stacking_functions

if __name__ == "__main__":
	conf, config_params = stacking_configparser.read_config_general();

	# Step 0: set up output directories
	stacking_functions.set_up_output_directories(config_params);

	# Step 1: make atmospheric corrections, etc.
	stacking_functions.make_corrections(config_params);

	# Step 2: Get reference information
	stacking_functions.get_ref(config_params);

	# Step 3: make velocity field
	stacking_functions.vels_and_ts(config_params);

	# Step 4: geocoding
	stacking_functions.geocode_vels(config_params);
