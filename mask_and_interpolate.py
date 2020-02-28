import numpy as np 
import scipy.interpolate



def cut_grid(data, xmin, xmax, ymin, ymax):
	return data[ymin:ymax, xmin:xmax];

def make_coherence_mask(cor, threshold):
	print("Making coherence mask.")
	mask = np.ones(np.shape(cor));
	for i in range(np.shape(cor)[0]):
		for j in range(np.shape(cor)[1]):
			if np.isnan(cor[i][j]) or cor[i][j]<threshold:
				mask[i][j]=np.nan;
	return mask;


def apply_coherence_mask(data, mask, is_complex=0):
	if is_complex==1:
		masked = np.complex64(np.multiply(data,mask));
	else:
		masked = np.float64(np.multiply(data,mask));
	return masked;

def interpolate_2d(data_array, is_complex=0):
	print("Performing 2d interpolation");
	if is_complex==1:
		data_array=np.angle(data_array);

	ymax, xmax = np.shape(data_array);
	yarray = range(ymax);
	xarray = range(xmax);
	interpolated_values = np.zeros(np.shape(data_array));

	x_interps=[]; y_interps=[]; z_interps=[]; xy_interps=[];
	xy_targets=[];

	# Time to get rid of the nan's. 
	for i in range(len(yarray)):
		for j in range(len(xarray)):
			# xy_targets.append([xarray[j], yarray[i]]);
			# Get the real values for use in interpolating
			if not np.isnan(data_array[i][j]):
				x_interps.append(xarray[j]);
				y_interps.append(yarray[i]);
				xy_interps.append([xarray[j], yarray[i]]);
				z_interps.append(data_array[i][j]);
			# Collect the points where we interpolate
			else:
				xy_targets.append([xarray[j], yarray[i]]);

	# f = scipy.interpolate.interp2d(y_interps, x_interps, z_interps);
	z_targets = scipy.interpolate.griddata(xy_interps, z_interps, xy_targets, method='linear', fill_value = 1);

	# Fill in the gaps using the interpolated value of z
	smoothdata=np.copy(data_array);
	for i in range(len(xy_targets)):
		idxx=xarray.index(xy_targets[i][0])
		idxy=yarray.index(xy_targets[i][1])
		if is_complex==1:
			r=1
			smoothdata[idxy][idxx] = np.cmath.rect(r,z_targets[i]);
		else:
			smoothdata[idxy][idxx] = z_targets[i];

	return smoothdata;