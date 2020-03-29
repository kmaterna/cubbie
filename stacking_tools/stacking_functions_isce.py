import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
import matplotlib.pyplot as plt 
from subprocess import call, check_output
import stacking_utilities
import readmytupledata
import stack_corr
import netcdf_read_write as rwr
import nsbas_isce
import isce_read_write
import mask_and_interpolate
import unwrapping_errors

# --------------- STEP 1: Make corrections and TS prep ------------ # 
def make_corrections_isce(config_params):
	if config_params.startstage>1:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<1:   # if we're ending before, we don't do this. 
		return;  
	print("Stage 1 - Doing optional corrections"); 

	# For ISCE, we might want to re-make all the interferograms and unwrap them in custom fashion. 
	# This operates on files in the Igram directory, no need to move directories yourself. 
	if config_params.solve_unwrap_errors:
		unwrapping_errors.main_function(config_params.rlks, config_params.alks, config_params.filt, 
			config_params.xbounds, config_params.ybounds, config_params.cor_cutoff_mask);

	# WE ALSO MAKE THE SIGNAL SPREAD FOR FULL IMAGES
	cor_value=0.5;
	filepathslist=glob.glob("../Igrams/????????_????????/filt*.cor");  # *** This may change
	cor_data = readmytupledata.reader_isce(filepathslist);
	a = stack_corr.stack_corr(cor_data, cor_value);
	rwr.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', config_params.ts_output_dir+'/signalspread_full.nc')
	rwr.produce_output_plot(config_params.ts_output_dir+'/signalspread_full.nc', 
		'Signal Spread above cor='+str(cor_value), config_params.ts_output_dir+'/signalspread_full.png', 
		'Percentage of coherence', aspect=1/4, invert_yaxis=False);

	return;


# --------------- UTILITIES ------------ # 
def get_100p_pixels_get_ref(filenameslist, ref_idx, outdir):
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
					if j>320 and j<337 and i>3500 and i<3700:
						xpixels_options.append(j);
						ypixels_options.append(i);

		idx_lucky = 800;
		xref = xpixels_options[idx_lucky];
		yref = ypixels_options[idx_lucky];

		print("%d of %d (%f percent) are totally coherent. " % (count, total_pixels, 100*(count/total_pixels)) );
		print(np.shape(mydata.zvalues));
		print("%d pixels are good options for the reference pixel. " % (len(xpixels_options)) );

		# Make a plot that shows where those pixels are
		# a=rwr.read_grd(outdir+'/signalspread_cut.nc');
		fig = plt.figure();
		# plt.imshow(a,aspect=1/4, cmap='rainbow');
		plt.plot(xpixels_good, ypixels_good, '.', color='k');
		plt.plot(xpixels_options,ypixels_options, '.', color='g');
		plt.plot(xref, yref, '.', color='r');
		plt.savefig(outdir+'/best_pixels.png');
		plt.close();

		print("Selecting reference pixel at row/col %d, %d " % (yref, xref) );

		return yref, xref;


def make_referenced_unwrapped_isce(intf_files, rowref, colref, output_dir):
	# This is simpler than the Sentinel case, since you don't have to cross 
	# swath borders or solve for n pi. 
	# We just return the values minus a single reference pixel, 
	# and write them to a directory. 
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
    call(['cp','stacking.config',config_params.ts_output_dir],shell=False);

    # Very general, takes all files and doesn't discriminate. 
    intf_files=stacking_utilities.get_list_of_intf_all(config_params);

    # This is an eaiser function, just one swath and one reference pixel. 
    rowref, colref = get_100p_pixels_get_ref(intf_files, config_params.ref_idx, config_params.ts_output_dir);

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
	intf_files = stacking_utilities.make_selection_of_intfs(config_params);
	call(['cp','stacking.config',config_params.ts_output_dir],shell=False);
	
	if config_params.ts_type=="STACK":
		print("Running velocities by simple stack.")
		# sss.drive_velocity_simple_stack(intfs, config_params.wavelength, config_params.ts_output_dir);
	if config_params.ts_type=="SBAS":
		print("Running velocities and time series by SBAS");
		# sbas.drive_velocity_sbas(config_params.swath, intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir);
		# sbas.drive_ts_sbas(config_params);
	if config_params.ts_type=="NSBAS":
		print("Running velocities and time series by NSBAS");
		nsbas_isce.drive_nsbas(config_params.swath, intf_files, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir);
	if config_params.ts_type=="WNSBAS":
		print("Running velocities and time series by WNSBAS");
		coh_files = stacking_utilities.make_selection_of_coh_files(config_params, intf_files);
		nsbas_isce.drive_nsbas(config_params.swath, intf_files, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir, coh_files);
	return; 



# --------------- STEP 4: Geocoding Velocities ------------ # 
def geocode_vels(config_params):
	if config_params.startstage>4:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<4:   # if we're ending before, we don't do this. 
		return; 
	# The goals here for UAVSAR:
	# Load lon/lat grids
	# Resample lon/lat grids
	# Cut them appropriately
	# Write them out in the output folder
	# Turn some images (vel, ts steps, etc.) into KMLs
	if config_params.SAT=="UAVSAR":
		geocode_UAVSAR_stack(config_params);
	return;


def geocode_UAVSAR_stack(config_params):
	# Collect initial information and set things up
	llh_array = np.fromfile(config_params.llh_file, dtype=np.float32);  # this is a vector. 
	lon=[]; lat=[]; hgt=[];
	lat=llh_array[np.arange(0,len(llh_array),3)];  # ordered array opened from the provided UAVSAR files
	lon=llh_array[np.arange(1,len(llh_array),3)];
	hgt=llh_array[np.arange(2,len(llh_array),3)];
	example_igram=glob.glob("../Igrams/????????_????????/*.int")[0];  
	phase_array = isce_read_write.read_phase_data(example_igram);
	print("Shape of the interferogram: ", np.shape(phase_array));

	# Determine the shape of the llh array 
	# assuming there's a giant gap somewhere in the lat array
	# that can tell us how many elements are in the gridded array
	typical_gap = abs(lat[1]-lat[0]);
	for i in range(1,len(lat)):
		if abs(lat[i]-lat[i-1]) > 100*typical_gap:
			print(lat[i]-lat[i-1]);
			print("There are %d columns in the lon/lat arrays" % i);
			llh_pixels_range=i;
			break;
	llh_pixels_azimuth = int(len(lon)/llh_pixels_range);
	print("llh_pixels_azimuth: ", llh_pixels_azimuth);
	print("llh_pixels_range: ",llh_pixels_range);

	# We turn the llh data into 2D arrays.
	lat_array = np.reshape(lat,(llh_pixels_azimuth,llh_pixels_range));
	lon_array = np.reshape(lon,(llh_pixels_azimuth,llh_pixels_range));

	# # write the data into a GDAL format. 
	isce_read_write.write_isce_data(lon_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", config_params.ts_output_dir+"/lon_total.gdal");
	isce_read_write.write_isce_data(lat_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", config_params.ts_output_dir+"/lat_total.gdal");
	temp = isce_read_write.read_scalar_data(config_params.ts_output_dir+"/lon_total.gdal");
	print("Shape of the lon file: ",np.shape(temp));

	# Resampling in GDAL
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',config_params.ts_output_dir+'/lon_total.gdal',
		config_params.ts_output_dir+'/lon_igram_res.tif'],shell=False);
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',config_params.ts_output_dir+'/lat_total.gdal',
		config_params.ts_output_dir+'/lat_igram_res.tif'],shell=False);


	# Cut the data, and quality check. 
	# Writing the cut lon/lat into new files. 
	temp_lon = isce_read_write.read_scalar_data(config_params.ts_output_dir+"/lon_igram_res.tif");
	print("Shape of the lon file: ",np.shape(temp_lon));
	xbounds=[float(config_params.xbounds.split(',')[0]),float(config_params.xbounds.split(',')[1])];
	ybounds=[float(config_params.ybounds.split(',')[0]),float(config_params.ybounds.split(',')[1])];
	cut_lon = mask_and_interpolate.cut_grid(temp_lon, xbounds, ybounds, fractional=True, buffer_rows=3);
	print("Shape of the cut lon file: ",np.shape(cut_lon));
	nx = np.shape(cut_lon)[1];
	ny = np.shape(cut_lon)[0];
	isce_read_write.write_isce_data(cut_lon, nx, ny, "FLOAT", config_params.ts_output_dir+'/cut_lon.gdal');
	rwr.produce_output_netcdf(np.array(range(0,nx)), np.array(range(0,ny)), cut_lon, 
		"degrees", config_params.ts_output_dir+'/cut_lon.nc');
	temp_lat = isce_read_write.read_scalar_data(config_params.ts_output_dir+"/lat_igram_res.tif");
	print("Shape of the lat file: ",np.shape(temp_lat));
	cut_lat = mask_and_interpolate.cut_grid(temp_lat, xbounds, ybounds, fractional=True, buffer_rows=3);
	print("Shape of the cut lat file: ",np.shape(cut_lat));
	isce_read_write.write_isce_data(cut_lat, nx, ny, "FLOAT", config_params.ts_output_dir+'/cut_lat.gdal');
	rwr.produce_output_netcdf(np.array(range(0,nx)), np.array(range(0,ny)), cut_lat, 
		"degrees", config_params.ts_output_dir+'/cut_lat.nc');

	# Double checking the shape of the interferogram data (should match!)
	signalspread = isce_read_write.read_scalar_data(config_params.ts_output_dir+'/signalspread_cut.nc');
	print("For comparison, shape of cut data is: ",np.shape(signalspread));

	isce_read_write.plot_scalar_data(config_params.ts_output_dir+'/cut_lat.gdal',
		colormap='rainbow',aspect=1/4,outname=config_params.ts_output_dir+'/cut_lat_geocoded.png');

	# NEXT: HERE WE WILL TURN THE GRID INTO A KML
	# cut_lon, cut_lat, ts, or whatever we feel like. 
	# We should have a function that turns individual images
	# Or time series into KML. 

	return;


