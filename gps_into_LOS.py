# Function that reads GPS field, 
# projects into LOS relative to a reflon and reflat
# And writes an output text file that can be used later for 
# plotting expected GPS LOS values. 

import numpy as np
import matplotlib.pyplot as plt 
import subprocess, sys
import datetime as dt 
import collections
from scipy import interpolate
import netcdf_read_write

Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);



def top_level_driver(config_params, rowref, colref):
	[gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, bounds, outfile]=configure(config_params);
	[gps_velfield, reflon, reflat] = inputs(gps_file, veldir, velfile, rowref, colref, bounds);
	[LOS_velfield] = compute(gps_velfield, reflon, reflat, flight_angle, look_angle, type_of_interp, bounds);
	outputs(gps_velfield, LOS_velfield, reflon, reflat, outfile);
	return;

# ----------------- CONFIGURE ----------------- # 
def configure(config_params):
	print("Starting gps_into_los.")
	gps_file=config_params.gps_file;
	veldir=config_params.ts_output_dir;
	velfile=veldir+'/vel.grd';
	flight_angle=config_params.flight_angle;
	look_angle=config_params.look_angle;
	outfile=veldir+'/gps_los.txt';
	bounds=[-125, -121, 38, 42.5]; # Good for Mendocino
	# bounds=[-125, -115, 35, 46]; # Good for WUS
	type_of_interp='linear';

	return [gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, bounds, outfile];


# ------------------ INPUTS --------------------- # 

def inputs(gps_file, veldir, velfile, rowref, colref, bounds):
	[gps_velfield]=read_pbo_vel_file(gps_file, bounds);
	[gps_velfield]=remove_duplicates(gps_velfield);
	[reflon, reflat]=generate_reflon_reflat(velfile, veldir, rowref, colref);
	return [gps_velfield, reflon, reflat];

def generate_reflon_reflat(velfile, veldir, rowref, colref):
	# In this part, I sometimes need to flip the x-axis of the input array to make sense with the geographic 
	# coordinates.
	# I suspect that for ascending orbits, this may not be necessary. 
	# Worth checking if it introduces bugs. 

	# Here we will use GMTSAR to geocode a small region including the reference pixel. 
	# We extract the latitude and longitude of the reference pixel. 
	refpoint_file='reference_point.grd';  # names only here. directory gets added later. 
	ref_ll_name = 'ref_ll';    # These are temporary files. 
	ref_ll = ref_ll_name+'.grd';

	[xdata, ydata,zdata]=netcdf_read_write.read_any_grd_xyz(velfile);

	colref=len(xdata)-colref;  # switching the x-axis direction. required for descending data, unsure about ascending. 
	# In general we can figure this out from the flight_angle. 

	rowarray=np.array([ydata[rowref],ydata[rowref+1]]);
	colarray=np.array([xdata[colref],xdata[colref+1]]);

	# plt.figure();
	# plt.contourf(xdata, ydata, zdata, cmap='jet');
	# plt.plot(xdata[colref],ydata[rowref],'.',markersize=20);
	# plt.gca().invert_xaxis();
	# plt.savefig('test1.eps');
	plt.figure();
	plt.imshow(zdata, vmin=-20, vmax=20, cmap='jet');
	plt.plot(len(xdata)-colref, rowref, '.',markersize=10, color='k');
	plt.savefig('refpoint.eps');

	zarray=np.array( [[0.0, 0.01], [0.01, 0.01 ]]);

	netcdf_read_write.produce_output_netcdf(colarray,rowarray,zarray, 'mm/yr', veldir+'/'+refpoint_file);
	netcdf_read_write.flip_if_necessary(veldir+'/'+refpoint_file);
	subprocess.call(['geocode_mod.csh',refpoint_file,ref_ll,ref_ll_name,veldir], shell=False);

	[xll, yll, zll]=netcdf_read_write.read_any_grd_variables(veldir+'/'+ref_ll, 'lon', 'lat', 'z');
	latref=yll[0];
	lonref=xll[0];
	print("\nReference Location is: " )
	print(lonref);
	print(latref);

	subprocess.call(['rm',veldir+'/'+ref_ll_name+'.png'],shell=False);
	subprocess.call(['rm',veldir+'/'+ref_ll_name+'.kml'],shell=False);

	return [lonref, latref];

def read_pbo_vel_file(infile, bounds):
# Meant for reading velocity files from the PBO/UNAVCO website. 
# Returns a Velfield collections object. 
	start=0;
	ifile=open(infile,'r');
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for line in ifile:
		if start==1:
			temp=line.split();
			name_temp=temp[0];
			nlat_temp=float(temp[7]);
			elon_temp=float(temp[8]);

			if elon_temp>180:
				elon_temp=elon_temp-360.0;

			if elon_temp>bounds[0] and elon_temp<bounds[1]:
				if nlat_temp>bounds[2] and nlat_temp<bounds[3]:
					name.append(name_temp);
					nlat.append(nlat_temp);
					elon.append(elon_temp);
					n.append(float(temp[19])*1000.0);
					e.append(float(temp[20])*1000.0);
					u.append(float(temp[21])*1000.0);
					sn.append(float(temp[22])*1000.0);
					se.append(float(temp[23])*1000.0);
					su.append(float(temp[24])*1000.0);
					t1=temp[-2];
					t2=temp[-1];
					first_epoch.append(dt.datetime.strptime(t1[0:8],'%Y%m%d'));
					last_epoch.append(dt.datetime.strptime(t2[0:8],'%Y%m%d'));
		if "*" in line:
			start=1;
	ifile.close();
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];

def remove_duplicates(velfield):
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	# These are the arrays that are growing. 
	
	for i in range(len(velfield.n)):
		is_duplicate = 0;
		for j in range(len(name)):
			# if abs(nlat[j]-velfield.nlat[i])<0.0005 and abs(elon[j]-velfield.elon[i])<0.0005:
			if velfield.name[i] in name:
				# we found a duplicate measurement. 
				is_duplicate = 1;
				# Right now assuming all entries at the same lat/lon have the same velocity values. 
			if velfield.name[i]=='P323':
				is_duplicate = 1; # A BAD STATION. 


		if is_duplicate == 0:
			name.append(velfield.name[i]);
			nlat.append(velfield.nlat[i]);
			elon.append(velfield.elon[i]);
			n.append(velfield.n[i]);
			sn.append(velfield.sn[i]);
			e.append(velfield.e[i]);
			se.append(velfield.se[i]);
			u.append(velfield.u[i]);
			su.append(velfield.su[i]);			
			first_epoch.append(velfield.first_epoch[i]);
			last_epoch.append(velfield.last_epoch[i]);

	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);	
	return [myVelfield];


# ------------------ COMPUTE --------------- # 

def compute(vel_tuple, reflon, reflat, flight_angle, look_angle, type_of_interp, bounds):


	# Scipy returns a function that you can use on a new set of x,y pairs. 
	f_east = interpolate.interp2d(vel_tuple.elon, vel_tuple.nlat, vel_tuple.e, kind=type_of_interp);
	f_north = interpolate.interp2d(vel_tuple.elon, vel_tuple.nlat, vel_tuple.n, kind=type_of_interp);


	new_east=[]; new_north=[]; new_vertical=[];
	# Evaluate the linear or cubic interpolation function at new points
	for i in range(len(vel_tuple.e)):
		new_east.append(f_east(vel_tuple.elon[i],vel_tuple.nlat[i]))  # only want to give the functions one point at a time. 
		new_north.append(f_north(vel_tuple.elon[i],vel_tuple.nlat[i]));
		new_vertical.append(0.0);  # there's no vertical deformation in this field by construction. 
	
	# Transform each field into the LOS fields 
	LOS_array=project_to_LOS(new_east,new_north,new_vertical,flight_angle,look_angle);

	# Take the reference point and transform its velocity into LOS. 
	velref_e=f_east(reflon,reflat);
	velref_n=f_north(reflon,reflat);
	LOS_reference=project_to_LOS(velref_e,velref_n,0,flight_angle,look_angle)[0];

	los_tuple=Velfield(name=vel_tuple.name, nlat=vel_tuple.nlat, elon=vel_tuple.elon, e=LOS_array-LOS_reference, n=0*LOS_array, u=0*LOS_array, 
		sn=vel_tuple.sn, se=vel_tuple.se, su=vel_tuple.su, first_epoch=vel_tuple.first_epoch, last_epoch=vel_tuple.last_epoch);

	return [los_tuple];



def project_to_LOS(U_e,U_n,U_u,flight_angle,incidence_angle):
	# Dlos = [U_n sin(phi) - U_e cos(phi)]*sin(lamda) + U_u cos(lamda)
	# [U_e, U_n, U_u] are the east, north, and up components of the deformation. 
	# phi   = azimuth of satellite heading vector, positive clockwise from north.
	# lamda = local incidence angle at the reflector.
	
	phi=np.deg2rad(flight_angle); 
	lamda=np.deg2rad(incidence_angle);

	if np.size(U_e)>1:  # processing a 1D array of values
		d_los = [0*i for i in U_e];
		for i in range(len(U_e)):
			d_los[i] = ( (U_n[i]*np.sin(phi) - U_e[i]*np.cos(phi) )*np.sin(lamda)  + U_u[i]*np.cos(lamda) );
	else:               # processing a single value
		d_los = (U_n*np.sin(phi) - U_e*np.cos(phi) )*np.sin(lamda)  + U_u*np.cos(lamda) ;

	# Flipping by negative one, because the convention is positive when moving away from the satellite. 
	d_los=np.array(d_los)*-1;

	return d_los



# ------------------- OUTPUTS ------------- # 
def outputs(gps_velfield, LOS_velfield, reflon, reflat, outfile):

	ofile=open(outfile,'w');
	for i in range(len(LOS_velfield.e)):
		ofile.write("%f %f %f %f %f %s \n" % (LOS_velfield.elon[i], LOS_velfield.nlat[i], gps_velfield.e[i], gps_velfield.n[i], LOS_velfield.e[i], LOS_velfield.name[i]) );
	ofile.write("%f %f %f %f %f %s\n" % (reflon, reflat, 0,0,0, 'reference') );
	ofile.close();
	print("Outputs printed to gps_los.txt");

	return;