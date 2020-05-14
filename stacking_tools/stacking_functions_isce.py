import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
import matplotlib.pyplot as plt 
from subprocess import call, check_output
import xml.etree.ElementTree as et
import stacking_utilities
import readmytupledata
import stack_corr
import coseismic_stack
import netcdf_read_write as rwr
import nsbas_accessing
import isce_read_write
import file_utilities
import mask_and_interpolate
import unwrapping_errors
import haversine

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
	# This function helps you manually choose a reference pixel. 
	# You might have to run through this function a number of times 
	# To select your boxes and your eventual reference pixel. 
	# I pick one in a stable area, outside of the deformation, ideally in a desert. 

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
				if j>280 and j<300 and i>2500 and i<2700:
					xpixels_options.append(j);
					ypixels_options.append(i);

	idx_lucky = 710;
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

	print("Based on 100p pixels, selecting reference pixel at row/col %d, %d " % (yref, xref) );
	print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.");

	return yref, xref;


def from_lonlat_get_rowcol(config_params):
	# Given a ref loc, get the geocoded grids and find the nearest point. 
	# Return its row and column. 
	# Step 1: Geocode properly based on a sample interferogram grid (done)
	# Step 2: extract nearest pixel (code in Brawley repo)
	geocode_UAVSAR_stack(config_params);

	reflon = float(config_params.ref_loc.split(',')[0]);
	reflat = float(config_params.ref_loc.split(',')[1]);
	# Next we get the nearest pixel from the rasters
	raster_lon = isce_read_write.read_scalar_data(config_params.ts_output_dir+"/cut_lon.gdal");
	raster_lat = isce_read_write.read_scalar_data(config_params.ts_output_dir+"/cut_lat.gdal");
	i_found, j_found = get_nearest_pixel_in_raster(raster_lon, raster_lat, reflon, reflat);
	print("From lon/lat, found Row and Column at %d, %d " % (i_found, j_found) );
	print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.");
	return i_found,j_found;


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat):
	# Take a raster and find the grid location closest to the target location
	dist = np.zeros(np.shape(raster_lon));
	lon_shape = np.shape(raster_lon);
	for i in range(lon_shape[0]):
		for j in range(lon_shape[1]):
			mypt = [raster_lat[i][j], raster_lon[i][j]];
			dist[i][j]=haversine.distance((target_lat, target_lon), mypt);
	minimum_distance = np.nanmin(dist);	
	if minimum_distance<0.25:  # if we're inside the domain.
		idx = np.where(dist==np.nanmin(dist));
		i_found = idx[0][0];
		j_found = idx[1][0];
		print(raster_lon[i_found][j_found], raster_lat[i_found][j_found]);
	else:
		i_found = -1;
		j_found = -1;  # error codes
	return i_found, j_found;



def stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, ts_output_dir,label=""):
	# WE MAKE THE SIGNAL SPREAD FOR THE CUT IMAGES
	cor_files = [i.replace("fully_processed.unwrappedphase","cut.cor") for i in intf_files]; # get for isce
	netcdfname = ts_output_dir+'/signalspread_cut_ref'+label+'.nc'
	cor_value=np.nan;
	cor_data = readmytupledata.reader_isce(cor_files);
	a = stack_corr.stack_corr(cor_data, cor_value);
	rwr.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', netcdfname)
	rwr.produce_output_plot(netcdfname, 'Signal Spread above cor='+str(cor_value), ts_output_dir+'/signalspread_cut_ref'+label+'.png', 
		'Percentage of coherence', aspect=1/4, invert_yaxis=False, dot_points=[[colref], [rowref]]);
	signal_spread_ref = a[rowref, colref];
	print("Signal Spread of the reference pixel = %.2f " % (signal_spread_ref) );
	if signal_spread_ref<50:
		print("WARNING: Your reference pixel has very low coherence. Consider picking a different one.");
		print("STOPPING ON PURPOSE.");
		sys.exit(0);
	return;


# --------------- STEP 2: Get Reference Pixel ------------ # 

def get_ref(config_params):
    if config_params.startstage>2:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<2:   # if we're ending at intf, we don't do this. 
        return;

    print("Stage 2 - Collecting referenced unwrapped.");
    call(['cp','stacking.config',config_params.ts_output_dir],shell=False);

    # Very general, takes all files and doesn't discriminate. 
    intf_files=stacking_utilities.get_list_of_intf_all(config_params);

    # If we are starting manually, we find reference pixel by using 100% pixels...
    if config_params.ref_idx =="" and config_params.ref_loc=="":
        rowref, colref = get_100p_pixels_get_ref(intf_files, config_params.ref_idx, config_params.ts_output_dir);
        sys.exit(0);
    # .... or by using a lon/lat to get the reference pixel
    if config_params.ref_idx=="" and config_params.ref_loc!="":  
        rowref, colref = from_lonlat_get_rowcol(config_params);
        sys.exit(0);
    # .... or we already know the reference indices:
    if config_params.ref_idx!="":  
        rowref = int(config_params.ref_idx.split("/")[0]);
        colref = int(config_params.ref_idx.split("/")[1]);
        print("From the config file, the Rowref and Colref are %d, %d\n" % (rowref, colref) );

    stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, config_params.ts_output_dir)

    return;


# --------------- STEP 3: Velocities and Time Series! ------------ # 

def vels_and_ts(config_params):
	if config_params.startstage>3:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<3:   # if we're ending at intf, we don't do this. 
		return;

	# This is where the hand-picking takes place. 
	# Ex: manual excludes, manual selects, long intfs only, ramp-removed, etc. 
	intf_files = stacking_utilities.make_selection_of_intfs(config_params);
	rowref=int(config_params.ref_idx.split('/')[0]);
	colref=int(config_params.ref_idx.split('/')[1]);
	call(['cp','stacking.config',config_params.ts_output_dir],shell=False);

	# Make signal_spread here. Can be commented if you already have it. 
	# stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, config_params.ts_output_dir, label='_selected');

	if config_params.ts_type=="STACK":
		print("Running velocities by simple stack.")
		sss.drive_velocity_simple_stack(intf_files, config_params.wavelength, rowref, colref, config_params.ts_output_dir);
	if config_params.ts_type=="COSEISMIC":
		print("Making a simple coseismic stack");
		coseismic_stack.drive_coseismic_stack_isce(intf_files, config_params.wavelength, rowref, colref, config_params.ts_output_dir); 
	if config_params.ts_type=="SBAS":
		print("Running velocities and time series by SBAS: SBAS currently broken.");
	if config_params.ts_type=="NSBAS":
		print("Running velocities and time series by NSBAS");
		nsbas_accessing.drive_full_TS_isce(intf_files, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir);
	if config_params.ts_type=="WNSBAS":
		print("Running velocities and time series by WNSBAS");
		coh_files = stacking_utilities.make_selection_of_coh_files(config_params, intf_files);
		nsbas_accessing.drive_full_TS_isce(intf_files, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir, coh_files);
	return; 



# --------------- STEP 4: Geocoding Velocities ------------ # 
def geocode_vels(config_params):
	if config_params.startstage>4:  # if we're starting after, we don't do this. 
		return;
	if config_params.endstage<4:   # if we're ending before, we don't do this. 
		return; 
	if config_params.SAT=="UAVSAR":
		geocode_directory=config_params.ts_output_dir+"/isce_geocode";
		# Deleting the contents of this folder would be a good automatic step in the future. 
		file_utilities.gmtsar_nc_2_isce_stack(config_params.ts_output_dir+"/TS.nc",geocode_directory, bands=2);  # write the TS data into isce binaries
		W, E, S, N = geocode_UAVSAR_stack(config_params, geocode_directory);  # do this once or more than once
		create_isce_unw_geo(geocode_directory, W, E, S, N);
		create_isce_rdr_geo(geocode_directory, W, E, S, N);
		inspect_isce(geocode_directory);
	return;

def cut_resampled_grid(outdir, filename, variable, config_params):
	# This is for metadata like lon, lat, and lookvector
	# Given an isce file and a set of bounds to cut the file, 
	# Produce the isce data and gmtsar netcdf that match each pixel. 
	temp = isce_read_write.read_scalar_data(outdir+"/"+filename);
	print("Shape of the "+variable+" file: ",np.shape(temp));
	xbounds=[float(config_params.xbounds.split(',')[0]),float(config_params.xbounds.split(',')[1])];
	ybounds=[float(config_params.ybounds.split(',')[0]),float(config_params.ybounds.split(',')[1])];
	cut_grid = mask_and_interpolate.cut_grid(temp, xbounds, ybounds, fractional=True, buffer_rows=3);
	print("Shape of the cut lon file: ",np.shape(cut_grid));
	nx = np.shape(cut_grid)[1];
	ny = np.shape(cut_grid)[0];
	isce_read_write.write_isce_data(cut_grid, nx, ny, "FLOAT", outdir+'/cut_'+variable+'.gdal');
	rwr.produce_output_netcdf(np.array(range(0,nx)), np.array(range(0,ny)), cut_grid, 
		"degrees", outdir+'/cut_'+variable+'.nc');
	return;

def normalize_look_vector(lkve, lkvn, lkvu):
	east_sq = np.square(lkve)
	north_sq = np.square(lkvn)
	up_sq = np.square(lkvu)
	sumarray = np.add(east_sq, north_sq)
	sumarray = np.add(sumarray, up_sq);
	magnitude = np.sqrt(sumarray);	
	norm_lkve = np.divide(lkve, magnitude)
	norm_lkvn = np.divide(lkvn, magnitude)
	norm_lkvu = np.divide(lkvu, magnitude)
	return norm_lkve, norm_lkvn, norm_lkvu;

def calc_isce_azimuth_incidence(lkve, lkvn, lkvu):
	# lkve, lkvn, lkvu describe vector from plane to ground
	# ISCE convention: Azimuth angle measured from North in Anti-clockwise direction
	# ISCE convention: Incidence angle measured from vertical at target (always +ve)
	east_sq = np.square(lkve);
	north_sq = np.square(lkvn);
	sumarray = np.add(east_sq, north_sq);
	magnitude = np.sqrt(sumarray);
	azimuth_standard=np.arctan2(-lkvn, -lkve);
	azimuth_standard = np.rad2deg(azimuth_standard);
	azimuth = np.add(azimuth_standard,-90);

	incidence = np.arctan2(magnitude, -lkvu);
	incidence = np.rad2deg(incidence);
	return azimuth, incidence;

def geocode_UAVSAR_stack(config_params, geocoded_folder):
	# The goals here for UAVSAR:
	# Load lon/lat grids and look vector grids
	# Resample and cut the grids appropriately
	# Write pixel-wise metadata out in the output folder
	# All these grids have only single band. 
	call(["mkdir","-p",geocoded_folder],shell=False);
	llh_array = np.fromfile(config_params.llh_file, dtype=np.float32);  # this is a vector. 
	lkv_array = np.fromfile(config_params.lkv_file, dtype=np.float32);
	lon=[]; lat=[]; hgt=[]; lkv_e = []; lkv_n = []; lkv_u = [];
	lat=llh_array[np.arange(0,len(llh_array),3)];  # ordered array opened from the provided UAVSAR files
	lon=llh_array[np.arange(1,len(llh_array),3)];
	hgt=llh_array[np.arange(2,len(llh_array),3)];
	lkv_e=lkv_array[np.arange(0,len(lkv_array),3)]
	lkv_n=lkv_array[np.arange(1,len(lkv_array),3)]
	lkv_u=lkv_array[np.arange(2,len(lkv_array),3)]
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
	# The look vector is in meters from the aircraft to the ground. 
	lat_array = np.reshape(lat,(llh_pixels_azimuth,llh_pixels_range));
	lon_array = np.reshape(lon,(llh_pixels_azimuth,llh_pixels_range));
	lkve_array = np.reshape(lkv_e, (llh_pixels_azimuth, llh_pixels_range));
	lkvn_array = np.reshape(lkv_n, (llh_pixels_azimuth, llh_pixels_range));
	lkvu_array = np.reshape(lkv_u, (llh_pixels_azimuth, llh_pixels_range));
	lkve_array, lkvn_array, lkvu_array = normalize_look_vector(lkve_array, lkvn_array, lkvu_array);
	azimuth, incidence = calc_isce_azimuth_incidence(lkve_array, lkvn_array, lkvu_array);

	# # write the data into a GDAL format. 
	isce_read_write.write_isce_data(lon_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/lon_total.gdal");
	isce_read_write.write_isce_data(lat_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/lat_total.gdal");
	# isce_read_write.write_isce_data(lkve_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/lkve_total.gdal");
	# isce_read_write.write_isce_data(lkvn_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/lkvn_total.gdal");
	# isce_read_write.write_isce_data(lkvu_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/lkvu_total.gdal");
	isce_read_write.write_isce_data(azimuth, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/azimuth_total.gdal");
	isce_read_write.write_isce_data(incidence, llh_pixels_range, llh_pixels_azimuth, "FLOAT", geocoded_folder+"/incidence_total.gdal");

	# Resampling in GDAL to match the interferogram sampling
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/lon_total.gdal',
		geocoded_folder+'/lon_igram_res.tif'],shell=False);
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/lat_total.gdal',
		geocoded_folder+'/lat_igram_res.tif'],shell=False);
	# call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
	# 	'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
	# 	'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/lkve_total.gdal',
	# 	geocoded_folder+'/lkve_igram_res.tif'],shell=False);
	# call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
	# 	'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
	# 	'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/lkvn_total.gdal',
	# 	geocoded_folder+'/lkvn_igram_res.tif'],shell=False);
	# call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
	# 	'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
	# 	'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/lkvu_total.gdal',
	# 	geocoded_folder+'/lkvu_igram_res.tif'],shell=False);
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/incidence_total.gdal',
		geocoded_folder+'/incidence_igram_res.tif'],shell=False);
	call(['gdalwarp','-ts',str(np.shape(phase_array)[1]),str(np.shape(phase_array)[0]),
		'-r','bilinear','-to','SRC_METHOD=NO_GEOTRANSFORM',
		'-to','DST_METHOD=NO_GEOTRANSFORM',geocoded_folder+'/azimuth_total.gdal',
		geocoded_folder+'/azimuth_igram_res.tif'],shell=False);

	# Cut the data, and quality check. 
	# Writing the cut lon/lat into new files. 
	cut_resampled_grid(geocoded_folder, "lon_igram_res.tif", "lon", config_params);
	cut_resampled_grid(geocoded_folder, "lat_igram_res.tif", "lat", config_params);
	# cut_resampled_grid(geocoded_folder, "lkve_igram_res.tif", "lkve", config_params);
	# cut_resampled_grid(geocoded_folder, "lkvn_igram_res.tif", "lkvn", config_params);
	# cut_resampled_grid(geocoded_folder, "lkvu_igram_res.tif", "lkvu", config_params);
	cut_resampled_grid(geocoded_folder, "incidence_igram_res.tif", "incidence", config_params);
	cut_resampled_grid(geocoded_folder, "azimuth_igram_res.tif", "azimuth", config_params);

	isce_read_write.plot_scalar_data(geocoded_folder+'/cut_lat.gdal',
		colormap='rainbow',aspect=1/4,outname=geocoded_folder+'/cut_lat_geocoded.png');
	cut_lon = isce_read_write.read_scalar_data(geocoded_folder+'/cut_lon.gdal');
	cut_lat = isce_read_write.read_scalar_data(geocoded_folder+'/cut_lat.gdal');
	W,E = np.min(cut_lon), np.max(cut_lon);
	S,N = np.min(cut_lat), np.max(cut_lat);

	# This last thing may not work when finding the reference pixel, only when geocoding at the very last. 
	# Double checking the shape of the interferogram data (should match!)
	signalspread = isce_read_write.read_scalar_data(config_params.ts_output_dir+'/signalspread_cut.nc');
	print("For comparison, shape of cut data is: ",np.shape(signalspread));

	return W, E, S, N;

def create_isce_unw_geo(geocoded_dir, W, E, S, N):
	# With pixel-wise lat and lon and lookvector information, 
	# Can we make isce geocoded unwrapped products? 
	# Goal 1: .unw.geo / .unw.geo.xml (seems to work)
	# Goal 2: los.rdr.geo / los.rdr.geo.xml
	# geocodeGdal.py -l cut_lat.gdal -L cut_lon.gdal -f cut_something.gdal -b "S N W E"
	# After that, the BIL arrangement can be switched to BSQ, 
	# So I need to make an adjustment
	folders = glob.glob(geocoded_dir+"/scene*");
	i=0;
	for folder_i in folders:
		# Run the geocode command. 
		# This places the geocoded .unw.geo into each sub-directory. 
		datafile = glob.glob(folder_i+"/*.unw");
		datafile = datafile[0]
		command = "geocodeGdal.py -l "+geocoded_dir+"/cut_lat.gdal -L "+geocoded_dir+"/cut_lon.gdal "+"-f "+datafile+" -b \""+str(S)+" "+str(N)+" "+str(W)+" "+str(E)+"\" -x 0.00025 -y 0.00025"
		print(command);
		print("\n");
		call(command,shell=True);

		# Unfortunately, after geocodeGdal, the files end up BSQ instead of BIL.  This is necessary to reshape them. 
		filename = datafile+".geo"
		isce_read_write.plot_scalar_data(filename,colormap='rainbow',datamin=-50, datamax=200,outname='test_after_geocode.png',band=2);
		nlat=508
		nlon=1325  # obviously will have to change depending on the situation. 
		disp=np.memmap(filename,dtype='<f4').reshape(nlat*2, nlon)  # The right way to read this file.
		band1 = disp[0:nlat,:];
		band2 = disp[nlat:,:];

		# Places the two images side by side . This is necessary for reading by Kite
		properdata = np.hstack((band1, band2));
		f2 = plt.plot();
		plt.imshow(properdata)
		plt.savefig('after_rearranging.png')
		plt.close();
		properdata.tofile(filename);

		# There's a bit of a metadata problem when I do this, but hey, as long as it works in Kite, right? 
		# If you read this as a normal two-band ISCE file, it'll get messed up. 
		# But Kite is fine.  Oh well. 
		isce_read_write.plot_scalar_data(filename,colormap='rainbow',datamin=-50, datamax=200,outname='test_after_geocode_band2.png',band=2);
		i=i+1;
	return;

def create_isce_rdr_geo(geocoded_dir, W, E, S, N):
	# Create a geocoded azimuth and geocoded incidence file
	# Then concatenate them into a two-band-file (los.rdr.geo)
	# Then update the xml metadata. 
	print("Creating los.rdr.geo")
	datafile=geocoded_dir+"/cut_azimuth.gdal"
	command = "geocodeGdal.py -l "+geocoded_dir+"/cut_lat.gdal -L "+geocoded_dir+"/cut_lon.gdal "+"-f "+datafile+" -b \""+str(S)+" "+str(N)+" "+str(W)+" "+str(E)+"\" -x 0.00025 -y 0.00025"
	print(command);
	print("\n");
	call(command,shell=True);
	datafile=geocoded_dir+"/cut_incidence.gdal"
	command = "geocodeGdal.py -l "+geocoded_dir+"/cut_lat.gdal -L "+geocoded_dir+"/cut_lon.gdal "+"-f "+datafile+" -b \""+str(S)+" "+str(N)+" "+str(W)+" "+str(E)+"\" -x 0.00025 -y 0.00025"
	call(command,shell=True);
	grid_inc = isce_read_write.read_scalar_data(geocoded_dir+"/cut_incidence.gdal.geo", flush_zeros=False);
	grid_az = isce_read_write.read_scalar_data(geocoded_dir+"/cut_azimuth.gdal.geo",flush_zeros=False);
	ny, nx = np.shape(grid_inc);
	filename=geocoded_dir+"/los.rdr.geo"
	isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", filename);
	return;

def inspect_isce(geocoded_dir):
	# Plot things. 
	folders = glob.glob(geocoded_dir+"/scene*");
	for folder_i in folders:
		datafile = glob.glob(folder_i+"/*.unw.geo");
		datafile = datafile[0];
		grid = isce_read_write.read_scalar_data(datafile, flush_zeros=False);
		print("Statistics:")
		print("shape: ",np.shape(grid))
		print("max: ",np.nanmax(grid))
		print("min: ",np.nanmin(grid))
		isce_read_write.plot_scalar_data(datafile, colormap="rainbow",datamin=-50,datamax=200,outname=folder_i+"/geocoded_data.png");
	return;	
