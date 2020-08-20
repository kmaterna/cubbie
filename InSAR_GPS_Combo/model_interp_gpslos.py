
""" 
This code models a cubic 2D interpolation on GPS velocities to 
generate a predicted field of InSAR velocities. 
It assumes one simple look angle and flight vector (doesn't change incidence angle over the InSAR range)
It uses Scipy interpolate: cubic or linear
matplotlib.path: can replicate the matlab inpolygon() function.
"""

import numpy as np 
import matplotlib.pyplot as plt 
import sys
import matplotlib.path as path
import collections
from scipy import interpolate
import gps_io_functions
import los_projection_math
import haversine

param_collection = collections.namedtuple("param_collection",['input_file','ca_file','or_file','type_of_interp','ascending_flight_angle','descending_flight_angle','incidence_angle','bounds','point1', 'point2', 'reference_point']);
Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']); 

def do_interpolation():
	my_param_collection = configure();
	[vel_tuple,ca_border,or_border] = inputs(my_param_collection.input_file, my_param_collection.bounds, my_param_collection.or_file, my_param_collection.ca_file);
	los_ascending_velfield, los_descending_velfield, eastnorthfield = compute(vel_tuple, my_param_collection, ca_border, or_border);
	outputs(vel_tuple,my_param_collection, ca_border, or_border, los_ascending_velfield, los_descending_velfield, eastnorthfield);


# -------------- CONFIG  --------------- # 
def configure():
	gps_input_file="../../../Mendocino_Geodesy/GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt"
	or_file = "../../../Misc/Mapping_Resources/oregon_bdr"
	ca_file = "../../../Misc/Mapping_Resources/california_bdr"
	type_of_interp="cubic"
	ascending_flight_angle=360-14;  # degrees from north (like strike)
	descending_flight_angle=180+14; # degrees from north (like strike)
	incidence_angle = 30;     # degrees from vertical (looking straight down is 0 degrees). 

	# # Oregon
	# bounds=[-125,-122, 41.5, 46.0];
	# gradient_point1=[-123.0,43.0];
	# gradient_point2=[-123.0,42.0];
	# reference_point=[];	

	# Mendocino
	bounds=[-125,-122, 39.0, 42.0];
	gradient_point1=[-124.0,40.0];
	gradient_point2=[-124.0,41.0];
	reference_point=[-123.4, 39.9];

	# # Bay Area
	# bounds=[-124,-121.3, 36.8, 39.0];
	# gradient_point1=[-122.3,37.2];
	# gradient_point2=[-121.5,37.6];
	# reference_point=[-121.5, 37.0];

	my_param_collection=param_collection(input_file=gps_input_file, ca_file=ca_file, or_file=or_file,
		type_of_interp=type_of_interp,
		ascending_flight_angle=ascending_flight_angle,descending_flight_angle=descending_flight_angle, incidence_angle=incidence_angle,
		bounds=bounds, point1=gradient_point1, point2=gradient_point2, reference_point=reference_point);
	return my_param_collection;


# -------------- INPUTS  --------------- # 
def inputs(input_file, bounds, or_file, ca_file):
	[vel_tuple] = gps_io_functions.read_pbo_vel_file(input_file);  # I actually have to interpolate over the whole field, shouldn't chunk it smaller with clean_velfield. 
	ca_border=np.loadtxt(ca_file);
	or_border=np.loadtxt(or_file);
	return [vel_tuple, ca_border, or_border];


# ---------- COMPUTE --------------- # 
# Steps: 
# Get interpolation points
# Interpolate based on function
# Project velocities into LOS
# Subtract reference point
# Report LOS gradient between points

def get_interpolation_points_within_grid(latlon_bounds, interval, optional_border_paths=[]):
	# The new interpolation grid: a new set of points with some chosen spacing (in degrees)
	# Respects path boundaries if they are provided (to cancel oceans, etc.)
	xarray=np.arange(latlon_bounds[0],latlon_bounds[1],interval);
	yarray=np.arange(latlon_bounds[2],latlon_bounds[3],interval);
	[X,Y]=np.meshgrid(xarray,yarray);
	x_for_interp=[]; y_for_interp=[]; 
	
	if len(optional_border_paths)>1:
		# We remove the interpolation points outside of CA and OR because they tend to explode. 
		for k in range(len(optional_border_paths)):
			include_path=path.Path(optional_border_paths[k]);
			for i in range(np.shape(X)[0]):
				for j in range(np.shape(X)[1]):
					if include_path.contains_point([X[i,j],Y[i,j]])==1:   # example: if the point is in CA or OR
						x_for_interp.append(X[i,j]);
						y_for_interp.append(Y[i,j]);
	else:
		for i in range(np.shape(X)[0]):
			for j in range(np.shape(X)[1]):
				x_for_interp.append(X[i,j]);
				y_for_interp.append(Y[i,j]);
	return [x_for_interp, y_for_interp];



def compute(vel_tuple,my_param_collection, ca_border, or_border):
	# Get the 1D lists of points we're going to interpolate (either a grid or a vector of GPS points)
	[x_for_interp, y_for_interp] = get_interpolation_points_within_grid(my_param_collection.bounds, 0.08, [ca_border, or_border]);

	# Interpolation functions that you can use on X-Y pairs
	f_east = interpolate.interp2d(vel_tuple.elon, vel_tuple.nlat, vel_tuple.e, kind=my_param_collection.type_of_interp);
	f_north = interpolate.interp2d(vel_tuple.elon, vel_tuple.nlat, vel_tuple.n, kind=my_param_collection.type_of_interp);

	# Evaluate the linear or cubic interpolation function at new points
	new_east=[]; new_north=[]; new_vertical=[];
	for i in range(len(x_for_interp)):
		new_east.append(f_east(x_for_interp[i],y_for_interp[i]))  # only want to give the functions one point at a time. 
		new_north.append(f_north(x_for_interp[i],y_for_interp[i]));
		new_vertical.append(0.0);  # there's no vertical deformation in this field by construction. 

	# Get the reference velocity in ENU
	velref_e, velref_n, velref_u = los_projection_math.get_point_enu_interp(my_param_collection.reference_point, f_east, f_north);

	# Transform each field into LOS fields with respect to a given pixel
	ascending_LOS_array=los_projection_math.simple_project_ENU_to_LOS(new_east,new_north,new_vertical,my_param_collection.ascending_flight_angle,my_param_collection.incidence_angle);
	ascending_LOS_reference=los_projection_math.simple_project_ENU_to_LOS(velref_e,velref_n,velref_u,my_param_collection.ascending_flight_angle,my_param_collection.incidence_angle)[0];
	los_ascending_velfield=Velfield(name=[], elon=x_for_interp, nlat=y_for_interp, e=ascending_LOS_array - ascending_LOS_reference,
		n=0*ascending_LOS_array, u=0*ascending_LOS_array, se=[], sn=[], su=[], first_epoch=[], last_epoch=[]);
	
	# Descending
	descending_LOS_array=los_projection_math.simple_project_ENU_to_LOS(new_east,new_north,new_vertical,my_param_collection.descending_flight_angle,my_param_collection.incidence_angle);
	descending_LOS_reference=los_projection_math.simple_project_ENU_to_LOS(velref_e,velref_n,velref_u,my_param_collection.descending_flight_angle,my_param_collection.incidence_angle)[0];
	los_descending_velfield=Velfield(name=[], elon=x_for_interp, nlat=y_for_interp, e=descending_LOS_array - descending_LOS_reference,
		n=0*descending_LOS_array, u=0*descending_LOS_array, se=[], sn=[], su=[], first_epoch=[], last_epoch=[]);

	# Packing up an object for returning
	eastnorthfield = Velfield(name=[], elon=x_for_interp, nlat=y_for_interp, e=new_east,
		n=new_north, u=0*descending_LOS_array, se=[], sn=[], su=[], first_epoch=[], last_epoch=[]);

	# Here I want to evaluate gradients at two hard-coded points. 
	e_def1, n_def1, u_def1 = los_projection_math.get_point_enu_interp(my_param_collection.point1, f_east, f_north);
	e_def2, n_def2, u_def2 = los_projection_math.get_point_enu_interp(my_param_collection.point2, f_east, f_north);
	ascending1=los_projection_math.simple_project_ENU_to_LOS(e_def1,n_def1,u_def1,my_param_collection.ascending_flight_angle,my_param_collection.incidence_angle);
	ascending2=los_projection_math.simple_project_ENU_to_LOS(e_def2,n_def2,u_def2,my_param_collection.ascending_flight_angle,my_param_collection.incidence_angle);
	descending1=los_projection_math.simple_project_ENU_to_LOS(e_def1,n_def1,u_def1,my_param_collection.descending_flight_angle,my_param_collection.incidence_angle);
	descending2=los_projection_math.simple_project_ENU_to_LOS(e_def2,n_def2,u_def2,my_param_collection.descending_flight_angle,my_param_collection.incidence_angle);
	evaluate_gradients(my_param_collection.point1, my_param_collection.point2, ascending1, ascending2, descending1, descending2);

	return los_ascending_velfield, los_descending_velfield, eastnorthfield;


def evaluate_gradients(point1, point2, ascending1, ascending2, descending1, descending2):
	mydistance = haversine.distance([point1[1],point1[0]],[point2[1],point2[0]]); 
	print("Gradient in ASCENDING track from point1 to point2: %f mm/yr in %f km " % (np.abs(ascending1-ascending2), mydistance) );
	print("Equal to: %f mm/yr per 100 km \n" % (100*np.abs(ascending1-ascending2)/mydistance) );
	print("Gradient in DESCENDING track from point1 to point2: %f mm/yr in %f km " % (np.abs(descending1-descending2), mydistance) );
	print("Equal to: %f mm/yr per 100 km \n" % (100*np.abs(descending1-descending2)/mydistance) );
	return;


# ---------- PLOTTING OUTPUTS --------------- # 
def outputs(vel_tuple, my_param_collection, ca_border, or_border, los_ascending_velfield, los_descending_velfield, eastnorthfield):

	# Making the plotting variables
	narray_x=np.array(eastnorthfield.elon);
	narray_y=np.array(eastnorthfield.nlat);
	narray_north=np.array(eastnorthfield.n)[:,0];
	narray_east=np.array(eastnorthfield.e)[:,0];
	asc_dataset=np.array(los_ascending_velfield.e)[:,0];
	des_dataset=np.array(los_descending_velfield.e)[:,0];

	# Looking for the max los deformation values, to use for color bars. 
	max_ascending_los, min_ascending_los = np.max(asc_dataset), np.min(asc_dataset);
	max_descending_los, min_descending_los = np.max(des_dataset), np.min(des_dataset);
	max_los_value = np.max([max_descending_los, max_ascending_los]);
	min_los_value = np.min([min_descending_los, min_ascending_los]);


	# Figure of interpolated velocities
	f1=plt.figure();
	ax=f1.add_subplot(1,1,1)
	ax.plot(eastnorthfield.elon,eastnorthfield.nlat,'.g');
	ax.quiver(eastnorthfield.elon, eastnorthfield.nlat, eastnorthfield.e, eastnorthfield.n, color='red',scale=500.0);
	f1 = my_plot_formatting(ax,my_param_collection.bounds,ca_border,or_border,"Interpolated GPS Velocity Field",vel_tuple);
	plt.savefig("Interpolated_field.png");
	plt.close();

	# Northward and Eastward velocities in a subplot
	f,axarr=plt.subplots(1,2,figsize=(12,6));
	h1=axarr[0].scatter(narray_x,narray_y,s=75,marker='s',c=narray_north,cmap='jet',edgecolors='face',vmin=-30,vmax=30);
	axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	my_plot_formatting(axarr[0],my_param_collection.bounds,ca_border,or_border,"North Velocity Interpolated from GPS",vel_tuple);
	h2=axarr[1].scatter(narray_x,narray_y,s=75,marker='s',c=narray_east,cmap='jet',edgecolors='face',vmin=-30,vmax=30);
	cbar=plt.colorbar(h1,ax=axarr[1]); cbar.set_label('mm/yr');
	axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	my_plot_formatting(axarr[1],my_param_collection.bounds,ca_border,or_border,"East Velocity Interpolated from GPS",vel_tuple);
	plt.savefig("East_North.png");
	plt.close();

	# ASCENDING AND DESCENDING VIEWING GEOMETRY
	f,axarr=plt.subplots(1,2,figsize=(12,6));
	axarr[0].scatter(narray_x,narray_y,s=75,marker='s',c=asc_dataset,cmap='jet',edgecolors='face',vmin=min_los_value, vmax=max_los_value);
	axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	quiver_point=[min(my_param_collection.bounds[0:2])+0.4, max(my_param_collection.bounds[2:])-0.4]
	axarr[0].quiver(quiver_point[0], quiver_point[1],np.cos(np.deg2rad(90-my_param_collection.ascending_flight_angle)),np.sin(np.deg2rad(90-my_param_collection.ascending_flight_angle)),scale=10);
	axarr[0].quiver(quiver_point[0], quiver_point[1],np.sin(np.deg2rad(90-my_param_collection.ascending_flight_angle)),-np.cos(np.deg2rad(90-my_param_collection.ascending_flight_angle)),scale=10,color='red');
	axarr[0].text(quiver_point[0]+0.05,quiver_point[1],"LOS",color='red',ha='right',va='top');
	axarr[0].plot(my_param_collection.point1[0],my_param_collection.point1[1],marker='s',color='black',markersize=5);
	axarr[0].plot(my_param_collection.point2[0],my_param_collection.point2[1],marker='s',color='black',markersize=5);
	axarr[0].plot(my_param_collection.reference_point[0],my_param_collection.reference_point[1],marker='s',color='black',markersize=5);
	my_plot_formatting(axarr[0],my_param_collection.bounds,ca_border,or_border,"Ascending GPS LOS Velocity (interp)",vel_tuple);
	
	h1=axarr[1].scatter(narray_x,narray_y,s=75,marker='s',c=des_dataset,cmap='jet',edgecolors='face',vmin=min_los_value,vmax=max_los_value);
	cbar=plt.colorbar(h1,ax=axarr[1]); cbar.set_label('mm/yr');
	quiver_point=[min(my_param_collection.bounds[0:2])+0.4, max(my_param_collection.bounds[2:])-0.4]
	axarr[1].quiver(quiver_point[0], quiver_point[1],np.cos(np.deg2rad(90-my_param_collection.descending_flight_angle)),np.sin(np.deg2rad(90-my_param_collection.descending_flight_angle)),scale=10);
	axarr[1].quiver(quiver_point[0], quiver_point[1],np.sin(np.deg2rad(90-my_param_collection.descending_flight_angle)),-np.cos(np.deg2rad(90-my_param_collection.descending_flight_angle)),scale=10,color='red');
	axarr[1].text(quiver_point[0]+0.05,quiver_point[1],"LOS",color='red',ha='left',va='bottom');
	axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	axarr[1].plot(my_param_collection.point1[0],my_param_collection.point1[1],marker='s',color='black',markersize=5);
	axarr[1].plot(my_param_collection.point2[0],my_param_collection.point2[1],marker='s',color='black',markersize=5);
	axarr[1].plot(my_param_collection.reference_point[0],my_param_collection.reference_point[1],marker='s',color='black',markersize=5);
	my_plot_formatting(axarr[1],my_param_collection.bounds,ca_border,or_border,"Descending GPS LOS Velocity (interp)",vel_tuple);

	# The sentinel descending scene. 
	xboxes=[-123.260666, -124.37, -123.98, -122.816429, -123.260666];
	yboxes=[39.701538, 39.86, 41.6, 41.437393, 39.701538];
	axarr[1].plot(xboxes, yboxes,'k');
	plt.savefig("LOS.png");

	return;


def my_plot_formatting(ax,bounds,ca_border,or_border,titlestring,vel_tuple):
	ax.set_xlim(bounds[0:2]);
	ax.set_ylim(bounds[2:]);
	ax.quiver(vel_tuple.elon, vel_tuple.nlat, vel_tuple.e, vel_tuple.n,scale=500.0);
	ax.plot(ca_border[:,0],ca_border[:,1],'k');
	ax.plot(or_border[:,0],or_border[:,1],'k');
	ax.plot(vel_tuple.elon, vel_tuple.nlat,'.');
	ax.set_title(titlestring);
	return ax;


if __name__=="__main__":
	do_interpolation();
