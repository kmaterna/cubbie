# July 2020
# Process a wrapped JPL UAVSAR interferogram from the UAVSAR website
# Here we read the grid file

import numpy as np 
import matplotlib.pyplot as plt 
import sys
import struct
import subprocess

def read_igram_data(data_file, ann_file, dtype='f',igram_type='ground'):
	# Note: The data file is not a regular netcdf grd file
	# dtype float works for corr, unwrapped, etc.
	# dtype double works for interferogram
	# igram_type is ground or slant
	print("Reading %s-range file %s" % (igram_type, data_file) );
	num_rows, num_cols = get_rows_cols(ann_file, igram_type);
	f = open(data_file,'rb')
	final_shape=(num_rows,num_cols);
	num_data = final_shape[0]*final_shape[1]*2;  # 2 for real/complex
	rawnum = f.read();
	f.close();
	floats = np.array(struct.unpack(dtype*num_data, rawnum))
	real = floats[::2];
	imag = floats[1::2];
	phase = np.arctan2(imag,real);
	data = phase.reshape(final_shape);
	return data;

def read_corr_data(data_file, ann_file, dtype='f',igram_type='ground'):
	# Note: The data file is not a regular netcdf grd file
	# dtype float works for corr, unwrapped, etc.
	# dtype double works for interferogram
	# igram_type is ground or slant
	print("Reading %s-range file %s" % (igram_type, data_file) );
	num_rows, num_cols = get_rows_cols(ann_file, igram_type);
	f = open(data_file,'rb')
	final_shape=(num_rows,num_cols);
	num_data = final_shape[0]*final_shape[1];  
	rawnum = f.read();
	f.close();
	floats = np.array(struct.unpack(dtype*num_data, rawnum))
	data = floats.reshape(final_shape);
	return data;	

def get_rows_cols(ann_file, igram_type):
	for line in open(ann_file):
		if igram_type=='ground':
			if 'Ground Range Data Latitude Lines' in line:
				temp=line.split('=');
				num_rows = int(temp[1]);
			if 'Ground Range Data Longitude Samples' in line:
				temp=line.split('=');
				num_cols = int(temp[1]);
		elif igram_type=='slant':
			if 'Slant Range Data Azimuth Lines' in line:
				temp=line.split('=');
				num_rows = int(temp[1]);
			if 'Slant Range Data Range Samples' in line:
				temp=line.split('=');
				num_cols = int(temp[1]);
		else:
			print("Error! Igram type not recognized");

	return num_rows, num_cols;

