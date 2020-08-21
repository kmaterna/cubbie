# August 2020 
# Calculate the misfit between a geocoded InSAR velocity field and a GPS velocity field 
# that has been projected into the Line of Sight

import numpy as np
import netcdf_read_write
import los_projection_tools


def top_level_driver(gps_los_file, geocoded_insar_file):
	[gps_los_velfield, xarray, yarray, LOS_array] = inputs(gps_los_file, geocoded_insar_file);
	compute(gps_los_velfield, xarray, yarray, LOS_array);
	return;

def inputs(gps_los_file, geocoded_insar_file):
	print("Reading files %s and %s for calculating misfit." % (gps_los_file, geocoded_insar_file) );
	[gps_los_velfield] = los_projection_tools.input_gps_as_los(gps_los_file);
	[xarray, yarray, LOS_array] = netcdf_read_write.read_netcdf4_xyz(geocoded_insar_file);
	if np.nanmean(xarray)>180:
		xarray = np.subtract(xarray,360);  # some files come in with 244 instead of -115.  Fixing that. 
	return [gps_los_velfield, xarray, yarray, LOS_array];

def compute(gps_los_velfield, xarray, yarray, LOS_array):
	misfit_array = los_projection_tools.misfit_gps_geocoded_insar(gps_los_velfield, xarray, yarray, LOS_array, window_pixels=15);
	print("Returning misfit values for %d desired GPS stations" % len(misfit_array) );
	rms_misfit = np.sqrt(np.mean(misfit_array**2));
	print("Results:");
	print("RMS Misfit Between these two fields is %f mm/yr\n" % rms_misfit);
	return;
