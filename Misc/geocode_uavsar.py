# The purpose of this code is to use the provided lat/lon grids
# with the UAVSAR stack products 
# to geocode some data into a grid and a KML. 

import numpy as np
from osgeo import gdal            ## GDAL support for reading virtual files
from read_write_insar_utilities import isce_read_write


def new_geo_transform(lon_grid, lat_grid):
	# Finds the GDAL transformation from the grid coordinates to the lat/lon coordinates in the associated grids. 
	# The lat/lon grids just contain values for lat and lon at each point. 
	# Can handle non-north and not-straight-oriented grids. 
	# Remember: XGeo = GT(0)+xpixel*GT(1)+yline*GT(2)
	#           YGeo = GT(3)+xpixel*GT(4)+yline*GT(5)

	num_rows = np.shape(lon_grid)[0];  # example: 8332
	num_cols = np.shape(lon_grid)[1];  # example: 4869
	gt0 = lon_array[0][0];
	gt3 = lat_array[0][0];

	# Setting up the matrix equation: 
	# Using the 0,0 corner and two other corners for the inversion. 
	# in row, column convention
	corner1 = (0, num_cols-1);
	corner2 = (num_rows-1, num_cols-1); 
	corner3 = (num_rows-1, 0);

	# In a pixel,line x,y convention, different from the array values, which are row,column
	matrix = [[num_cols-1, 0, 0, 0],[0, 0, num_cols-1, 0],[0, num_rows-1, 0, 0],[0, 0, 0, num_rows-1]]; 
	data = [[lon_array[corner1]-gt0],[lat_array[corner1]-gt3],[lon_array[corner3]-gt0],[lat_array[corner3]-gt3]];
	inv_matrix = np.linalg.inv(matrix);
	transform_values = np.matmul(inv_matrix, data);

	gt1=transform_values[0][0];
	gt2=transform_values[1][0];
	gt4=transform_values[2][0];
	gt5=transform_values[3][0];

	new_transform = (gt0, gt1, gt2, gt3, gt4, gt5);
	print("New geotransform: ");
	print(new_transform);

	print("Checking sanity of this transform. ");

	test_point(lon_array, lat_array, (0,0), new_transform);
	test_point(lon_array, lat_array, (0, num_cols-1), new_transform);
	test_point(lon_array, lat_array, (num_rows-1, 0), new_transform);
	test_point(lon_array, lat_array, (num_rows-1, num_cols-1), new_transform);

	return new_transform;


def test_point(lon_array, lat_array, index, new_transform):
	# index is in row, col coordinates
	# xpixel and ypixel are in col, row coordinates. 
	print("Testing point: ");
	print("row/col index: ", index);
	print("input:", lon_array[index], lat_array[index]);
	print("return:", test_transform(new_transform, index[1], index[0]))
	print();
	return;


def test_transform(GT, xpixel, ypixel):
	xgeo = GT[0]+xpixel*GT[1]+ypixel*GT[2];
	ygeo = GT[3]+xpixel*GT[4]+ypixel*GT[5];
	return xgeo, ygeo;


if __name__=="__main__":
	date_string = "20100412_20130524"
	rlks=10;
	alks=10;  # for the interferogram. 
	rlks_llh = 2;  # for the downloaded llh file
	alks_llh = 8;

	llh_file = "/home/kmaterna/Documents/Brawley/Data/26509/SanAnd_26509_01_BC_s2_2x8.llh";

	llh_array = np.fromfile(llh_file, dtype=np.float32)
	phase_array = isce_read_write.read_phase_data(date_string + ".int");
	print("Shape of the interferogram: ", np.shape(phase_array));

	downlooked_pixels = np.shape(phase_array);
	total_pixels_azimuth = downlooked_pixels[0]*alks;  # this corresponds to the 1x1 case. 
	total_pixels_range = downlooked_pixels[1]*rlks;    # this corresponds to the 1x1 case. 
	llh_pixels_azimuth = int(np.ceil(total_pixels_azimuth / alks_llh)); # this is how many we have in the downsampled array. 
	llh_pixels_range = int(np.ceil(total_pixels_range / rlks_llh));

	lon=[]; lat=[]; hgt=[];
	num_pixels = len(llh_array);
	lat=llh_array[np.arange(0,num_pixels,3)];
	lon=llh_array[np.arange(1,num_pixels,3)];
	hgt=llh_array[np.arange(2,num_pixels,3)];

	# We turn the llh data into 2D arrays.
	lat_array = np.reshape(lat,(llh_pixels_azimuth,llh_pixels_range));
	lon_array = np.reshape(lon,(llh_pixels_azimuth,llh_pixels_range));
	hgt_array = np.reshape(hgt,(llh_pixels_azimuth,llh_pixels_range));

	# write the hgt into a GDAL format. 
	isce_read_write.write_isce_data(np.abs(hgt_array), llh_pixels_range, llh_pixels_azimuth, "FLOAT", "hgt.gdal");

	# Reading the gdal array. 
	ds = gdal.Open("hgt.gdal", gdal.GA_ReadOnly)
	data = ds.GetRasterBand(1).ReadAsArray()
	transform = ds.GetGeoTransform()
	print("Original geotransform: ")
	print(transform);
	# isce_read_write.plot_scalar_data("hgt.gdal", outname='hgt.png');

	# Get a new geotransform. This will turn pixel/line numbers into the proper coordinates on the ground. 
	new_transform = new_geo_transform(lon_array, lat_array);

	# NEXT: RESAMPLE IN THE SAME WAY AS THE INTERFEROGRAMS
	# NEXT: CUT THE GRID IN THE SAME WAY AS THE INTERFEROGRAMS

	# Write the original data back into a GDAL raster that has the right geotransform attached. 
	format="GTiff"
	driver = gdal.GetDriverByName( format );
	target_ds = driver.Create("test.tiff", np.shape(hgt_array)[1], np.shape(hgt_array)[0], 1, gdal.GDT_Byte);
	target_ds.SetGeoTransform(new_transform);
	target_ds.GetRasterBand(1).WriteArray(np.abs(hgt_array));
	target_ds = None

# A nice command, but it doesnt apply the translation except GT0 and GT3. 
# gdaldem color-relief -nearest_color_entry -co format=png -of KMLSUPEROVERLAY test.tiff color.txt raster-kml.kmz