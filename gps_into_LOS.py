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
	[gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, bounds_latlon, bounds_ref, outfile]=configure(config_params);
	[gps_velfield, gps_velfield_removed, reflon, reflat] = inputs(gps_file, veldir, velfile, rowref, colref, bounds_latlon, bounds_ref);
	[LOS_velfield] = compute(gps_velfield, gps_velfield_removed, reflon, reflat, flight_angle, look_angle, type_of_interp, bounds_ref);
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
	outfile=veldir+'/gps_ll_enu_los.txt';
	bounds_latlon=[-125, -121, 38, 42.5]; # Good for Mendocino
	# bounds=[-125, -115, 35, 46]; # Good for WUS
	nominal_boundary=0.4;
	bounds_ref=[-nominal_boundary, nominal_boundary, -nominal_boundary, nominal_boundary]; # in Longitude/Latitude units away from the reference pixel. 
	# I have found that it is very important to restrict the spline's domain. 
	# Otherwise it can get very unstable. 
	type_of_interp='linear';

	return [gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, bounds_latlon, bounds_ref, outfile];


# ------------------ INPUTS --------------------- # 

def inputs(gps_file, veldir, velfile, rowref, colref, bounds_latlon, bounds_ref):
	[gps_velfield]=read_unr_vel_file(gps_file, bounds_latlon);
	[gps_velfield]=remove_duplicates(gps_velfield);
	[reflon, reflat]=generate_reflon_reflat(velfile, veldir, rowref, colref);
	[gps_velfield_removed]=remove_farfield(gps_velfield, reflon, reflat, bounds_ref);
	return [gps_velfield, gps_velfield_removed, reflon, reflat];

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

	# Flipping the x-axis direction and the data itself. Required for descending data, unsure about ascending. 
	colref=len(xdata)-1-colref; 
	# rowref=len(ydata)-1-rowref;
	zdata = np.fliplr(zdata);
	# In general we can figure this out from the flight_angle. 

	print("\nHello! Your reference pixel is (row,col) = (%d, %d)" % (rowref, colref) );
	print("Its velocity is %.2f mm/yr\n" % zdata[rowref][colref]);
	print("Its azimuth is %.2f " % ydata[rowref])
	print("Its range is %.2f \n\n" % xdata[colref])

	rowarray=np.array([ydata[rowref],ydata[rowref+1]]);
	colarray=np.array([xdata[colref],xdata[colref+1]]);

	plt.figure();
	plt.imshow(zdata, vmin=-20, vmax=20, cmap='jet');
	plt.plot(colref,rowref, '.',markersize=10, color='k');
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

def read_unr_vel_file(infile,bounds):
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		if temp[0]=="#":
			continue;
		name_temp=temp[-1];
		nlat_temp=float(temp[1]);
		elon_temp=float(temp[0]);
		if elon_temp>bounds[0] and elon_temp<bounds[1]:
			if nlat_temp>bounds[2] and nlat_temp<bounds[3]:
				name.append(name_temp);
				nlat.append(nlat_temp);
				elon.append(elon_temp);
				n.append(float(temp[3]));
				e.append(float(temp[2]));
				u.append(float(temp[4]));
				sn.append(float(temp[6]));
				se.append(float(temp[5]));
				su.append(float(temp[7]));
				first_epoch.append(dt.datetime.strptime(temp[-3],'%Y%m%d'));
				last_epoch.append(dt.datetime.strptime(temp[-2],'%Y%m%d'));
	ifile.close();
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield]; 

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

def remove_farfield(velfield, reflon, reflat, bounds):
	# Expects the bounds here to be relative to the reference pixel (in degrees away)
	new_name=[]; new_nlat=[]; new_elon=[]; new_n=[]; new_e=[]; new_u=[]; new_sn=[]; new_se=[]; new_su=[]; new_first_epoch=[]; new_last_epoch=[];
	for i in range(len(velfield.name)):
		if velfield.nlat[i]>reflat+bounds[2] and velfield.nlat[i]<reflat+bounds[3]:
			if velfield.elon[i]>reflon+bounds[0] and velfield.elon[i]<reflon+bounds[1]:
				# The station is within the box.
				new_name.append(velfield.name[i]);
				new_nlat.append(velfield.nlat[i]);
				new_elon.append(velfield.elon[i]);
				new_n.append(velfield.n[i]);
				new_e.append(velfield.e[i]);
				new_u.append(velfield.u[i]);
				new_sn.append(velfield.sn[i]);
				new_se.append(velfield.se[i]);
				new_su.append(velfield.su[i]);
				new_first_epoch.append(velfield.first_epoch[i]);
				new_last_epoch.append(velfield.last_epoch[i]);
	print("Restricting velfield of %d GPS points to only nearby %d points. " % (len(velfield.name), len(new_name)) );
	myVelfield=Velfield(name=new_name, nlat=new_nlat, elon=new_elon, n=new_n, e=new_e, u=new_u, sn=new_sn, se=new_se, su=new_su, first_epoch=new_first_epoch, last_epoch=new_last_epoch);
	return [myVelfield];

# ------------------ COMPUTE --------------- # 

def compute(vel_tuple, vel_tuple_removed, reflon, reflat, flight_angle, look_angle, type_of_interp, bounds):

	# Scipy returns a function that you can use on a new set of x,y pairs. 
	f_east = interpolate.interp2d(vel_tuple_removed.elon, vel_tuple_removed.nlat, vel_tuple_removed.e, kind=type_of_interp);
	f_north = interpolate.interp2d(vel_tuple_removed.elon, vel_tuple_removed.nlat, vel_tuple_removed.n, kind=type_of_interp);

	plot_bounds=[reflon+bounds[0], reflon+bounds[1], reflat+bounds[2], reflat+bounds[3]]
	make_interp_checking_plot(vel_tuple_removed, plot_bounds, reflon, reflat, f_east, f_north);
	
	# Transform each field into the LOS fields 
	# Assuming zero vertical velocity. 
	U_u=np.zeros(np.shape(vel_tuple.e));  # only for 2D case
	LOS_array=project_to_LOS(vel_tuple.e,vel_tuple.n,vel_tuple.u,flight_angle,look_angle);

	# Take the reference point and transform its velocity into LOS. 
	velref_e=f_east(reflon,reflat);
	velref_n=f_north(reflon,reflat);
	LOS_reference=project_to_LOS(velref_e,velref_n,0,flight_angle,look_angle)[0];

	print("Reference e: %f" % velref_e)
	print("Reference n: %f" % velref_n)
	print("Reference LOS: " + str(LOS_reference))
	print("-----------------------")

	print("Nearby: Station lat lon LOS east north")
	for i in range(len(vel_tuple.name)):
		if vel_tuple.elon[i]>-124 and vel_tuple.elon[i]<-123.0:
			if vel_tuple.nlat[i]>40.0 and vel_tuple.nlat[i]<40.7:
				print("%s %f %f %f %f %f" % (vel_tuple.name[i], vel_tuple.nlat[i], vel_tuple.elon[i], LOS_array[i]-LOS_reference, vel_tuple.e[i], vel_tuple.n[i]) )


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
		ofile.write("%f %f %f %f %f %f %s \n" % (LOS_velfield.elon[i], LOS_velfield.nlat[i], gps_velfield.e[i], gps_velfield.n[i], gps_velfield.u[i], LOS_velfield.e[i], LOS_velfield.name[i]) );
	ofile.write("%f %f %f %f %f %f %s\n" % (reflon, reflat, 0,0,0,0, 'reference') );
	ofile.close();
	print("Outputs printed to %s" % outfile);

	return;



def make_interp_checking_plot(vel_tuple, bounds, reflon, reflat, f_east, f_north):
	# A debugging plot
	# I've found that interp2d doesn't always give reasonable numbers. 
	# Although sometimes it does. It's worth double checking. 
	# Compute the interpolation everywhere in the bounds.
	xarray=np.arange(bounds[0],bounds[1],0.20);
	yarray=np.arange(bounds[2],bounds[3],0.20);
	[X,Y]=np.meshgrid(xarray,yarray);

	new_x=[]; new_y=[]; new_east=[]; new_north=[]; new_vertical=[];
	for i in range(np.shape(X)[0]):
		for j in range(np.shape(X)[1]):
			new_x.append(X[i,j]);
			new_y.append(Y[i,j]);


	# Evaluate the linear or cubic interpolation function at new points
	for i in range(len(new_x)):
		new_east.append(f_east(new_x[i],new_y[i]))  # only want to give the functions one point at a time. 
		new_north.append(f_north(new_x[i],new_y[i]));
		new_vertical.append(0.0); # there's no vertical deformation in this field by construction. 


	velref_e=f_east(reflon,reflat);
	velref_n=f_north(reflon,reflat);

	# Making the plotting variables
	narray_x=np.array(new_x);
	narray_y=np.array(new_y);
	narray_north=np.array(new_north)[:,0];
	narray_east=np.array(new_east)[:,0];

	# Northward and Eastward velocities in a subplot
	f,axarr=plt.subplots(1,2,figsize=(12,6));
	h1=axarr[0].scatter(narray_x,narray_y,s=575,marker='s',c=narray_north,cmap='jet',edgecolors='face',vmin=-30,vmax=30);
	axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	axarr[0].quiver(reflon,reflat, velref_e, velref_n, color='purple',scale=500);
	my_plot_formatting(axarr[0],bounds,"North Velocity Interpolated from GPS",vel_tuple);
	h2=axarr[1].scatter(narray_x,narray_y,s=575,marker='s',c=narray_east,cmap='jet',edgecolors='face',vmin=-30,vmax=30);
	cbar=plt.colorbar(h1,ax=axarr[1]); cbar.set_label('mm/yr');
	axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white',scale=500.0);
	axarr[1].quiver(reflon,reflat, velref_e, velref_n, color='purple',scale=500);
	my_plot_formatting(axarr[1],bounds,"East Velocity Interpolated from GPS",vel_tuple);
	plt.savefig("East_North_GPS_interp.png");
	plt.close();

	return; 


def my_plot_formatting(ax,bounds,titlestring,vel_tuple):
	ax.set_xlim([bounds[0],bounds[1]]);
	ax.set_ylim([bounds[2],bounds[3]]);
	ax.quiver(vel_tuple.elon, vel_tuple.nlat, vel_tuple.e, vel_tuple.n,scale=500.0);
	ax.plot(vel_tuple.elon, vel_tuple.nlat,'.');
	ax.set_title(titlestring);
	return ax;