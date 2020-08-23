# Function that reads GPS field, 
# projects very simply into LOS, and subtracts a reflon and reflat
# Writes an output text file that can be plotted in GMT, etc. 
# For this driver, the reference point is a lat/lon pixel, not a particular GPS station 
# (sometimes the situation isn't ideal and you have to do that!)
# This script has been partly re-written but is still a bit dependent on earlier parts of my processing pipeline

import numpy as np
import matplotlib.pyplot as plt
import subprocess, sys
import datetime as dt
from scipy import interpolate
import netcdf_read_write
import gps_io_functions
import los_projection_tools


def top_level_driver(config_params, rowref, colref):
    [gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, coordbox_gps, coordbox_nearref,
     outfile] = configure(config_params);
    [gps_velfield, gps_velfield_removed, reflon, reflat] = inputs(gps_file, veldir, velfile, rowref, colref,
                                                                  coordbox_gps, coordbox_nearref);
    [LOS_velfield] = compute(gps_velfield, gps_velfield_removed, reflon, reflat, flight_angle, look_angle,
                             type_of_interp, coordbox_nearref);
    outputs(gps_velfield, LOS_velfield, reflon, reflat, outfile);
    return;


# ----------------- CONFIGURE ----------------- #
def configure(config_params):
    print("Starting gps_into_los.")
    gps_file = config_params.gps_file;
    veldir = config_params.ts_output_dir;
    velfile = veldir + '/vel.grd';
    flight_angle = config_params.flight_angle;
    look_angle = config_params.look_angle;
    outfile = veldir + '/gps_ll_enu_los.txt';
    coordbox_gps = [-125, -121, 38, 42.5];  # Good for Mendocino
    # bounds = [-125, -115, 35, 46]; # Good for WUS
    nominal_boundary = 0.4;
    bounds_ref = [-nominal_boundary, nominal_boundary, -nominal_boundary,
                  nominal_boundary];  # in Longitude/Latitude units away from the reference pixel.
    coordbox_nearref = [reflon + bounds_ref[0], reflon + bounds_ref[1], reflat + bounds_ref[2], reflat + bounds_ref[3]];
    # For interpolating around the reference, I have found that it is very important to restrict the spline's domain. Otherwise it can get very unstable.
    type_of_interp = 'linear';

    return [gps_file, veldir, velfile, flight_angle, look_angle, type_of_interp, coordbox_gps, coordbox_nearref,
            outfile];


# ------------------ INPUTS --------------------- # 

def inputs(gps_file, veldir, velfile, rowref, colref, coordbox_gps, coordbox_nearref):
    [reflon, reflat] = generate_reflon_reflat(velfile, veldir, rowref, colref);

    # The velocities within the latlon box.
    [gps_velfield] = gps_io_functions.read_unr_vel_file(gps_file);
    [gps_velfield] = gps_io_functions.remove_duplicates(gps_velfield);
    [gps_velfield] = gps_io_functions.clean_velfield(gps_velfield, coord_box=coordbox_gps);

    # A small range near the reference pixel for interpolating later.
    [gps_velfield_removed] = gps_io_functions.clean_velfield(gps_velfield, coord_box=coordbox_nearref);

    return [gps_velfield, gps_velfield_removed, reflon, reflat];


def generate_reflon_reflat(velfile, veldir, rowref, colref):
    # In this part, I sometimes need to flip the x-axis of the input array to make sense with the geographic coordinates.
    # I suspect that for ascending orbits, this may not be necessary.
    # Worth checking if it introduces bugs.

    # Here we will use GMTSAR to geocode a small region including the reference pixel.
    # We extract the latitude and longitude of the reference pixel.
    refpoint_file = 'reference_point.grd';  # names only here. directory gets added later.
    ref_ll_name = 'ref_ll';  # These are temporary files.
    ref_ll = ref_ll_name + '.grd';

    [xdata, ydata, zdata] = netcdf_read_write.read_any_grd_xyz(velfile);

    # Flipping the x-axis direction and the data itself. Required for descending data, unsure about ascending.
    # All of this will change with better grid referencing in the future.
    colref = len(xdata) - 1 - colref;
    # rowref = len(ydata)-1-rowref;
    zdata = np.fliplr(zdata);
    # In general we can figure this out from the flight_angle.

    print("\nHello! Your reference pixel is (row,col) = (%d, %d)" % (rowref, colref));
    print("Its velocity is %.2f mm/yr\n" % zdata[rowref][colref]);
    print("Its azimuth is %.2f " % ydata[rowref])
    print("Its range is %.2f \n\n" % xdata[colref])

    rowarray = np.array([ydata[rowref], ydata[rowref + 1]]);
    colarray = np.array([xdata[colref], xdata[colref + 1]]);

    plt.figure();
    plt.imshow(zdata, vmin=-20, vmax=20, cmap='jet');
    plt.plot(colref, rowref, '.', markersize=10, color='k');
    plt.savefig('refpoint.eps');

    zarray = np.array([[0.0, 0.01], [0.01, 0.01]]);

    netcdf_read_write.produce_output_netcdf(colarray, rowarray, zarray, 'mm/yr', veldir + '/' + refpoint_file);
    netcdf_read_write.flip_if_necessary(veldir + '/' + refpoint_file);
    subprocess.call(['geocode_mod.csh', refpoint_file, ref_ll, ref_ll_name, veldir], shell=False);

    [xll, yll, zll] = netcdf_read_write.read_any_grd_variables(veldir + '/' + ref_ll, 'lon', 'lat', 'z');
    latref = yll[0];
    lonref = xll[0];
    print("\nReference Location is: ", lonref, latref)

    subprocess.call(['rm', veldir + '/' + ref_ll_name + '.png'], shell=False);
    subprocess.call(['rm', veldir + '/' + ref_ll_name + '.kml'], shell=False);

    return [lonref, latref];


# ------------------ COMPUTE --------------- #

def compute(vel_tuple, vel_tuple_removed, reflon, reflat, flight_angle, look_angle, type_of_interp, bounds):
    # Scipy returns a function that you can use on a new set of x,y pairs.
    f_east = interpolate.interp2d(vel_tuple_removed.elon, vel_tuple_removed.nlat, vel_tuple_removed.e,
                                  kind=type_of_interp);
    f_north = interpolate.interp2d(vel_tuple_removed.elon, vel_tuple_removed.nlat, vel_tuple_removed.n,
                                   kind=type_of_interp);

    make_interp_checking_plot(vel_tuple_removed, bounds, reflon, reflat, f_east, f_north);

    # Take the reference point and transform its velocity into LOS.
    velref_e, velref_n, velref_u = los_projection_tools.get_point_enu_interp([reflon, reflat], f_east=f_east,
                                                                             f_north=f_north)
    LOS_reference = \
    los_projection_tools.simple_project_ENU_to_LOS(velref_e, velref_n, velref_u, flight_angle, look_angle)[0];

    # Transform GPS field into LOS field
    LOS_array = los_projection_tools.simple_project_ENU_to_LOS(vel_tuple.e, vel_tuple.n, vel_tuple.u, flight_angle,
                                                               look_angle);

    los_tuple = los_projection_tools.Velfield(name=vel_tuple.name, nlat=vel_tuple.nlat, elon=vel_tuple.elon,
                                              e=LOS_array - LOS_reference, n=0 * LOS_array, u=0 * LOS_array,
                                              sn=vel_tuple.sn, se=vel_tuple.se, su=vel_tuple.su,
                                              first_epoch=vel_tuple.first_epoch, last_epoch=vel_tuple.last_epoch);

    return [los_tuple];


# ------------------- OUTPUTS ------------- #
def outputs(gps_velfield, LOS_velfield, reflon, reflat, outfile):
    ofile = open(outfile, 'w');
    for i in range(len(LOS_velfield.e)):
        ofile.write("%f %f %f %f %f %f %s \n" % (
        LOS_velfield.elon[i], LOS_velfield.nlat[i], gps_velfield.e[i], gps_velfield.n[i], gps_velfield.u[i],
        LOS_velfield.e[i], LOS_velfield.name[i]));
    ofile.write("%f %f %f %f %f %f %s\n" % (reflon, reflat, 0, 0, 0, 0, 'reference'));
    ofile.close();
    print("Outputs printed to %s" % outfile);

    return;


def make_interp_checking_plot(vel_tuple, bounds, reflon, reflat, f_east, f_north):
    # A debugging plot
    # I've found that interp2d doesn't always give reasonable numbers.
    # Although sometimes it does. It's worth double checking.
    # Compute the interpolation everywhere in the bounds.
    xarray = np.arange(bounds[0], bounds[1], 0.20);
    yarray = np.arange(bounds[2], bounds[3], 0.20);
    [X, Y] = np.meshgrid(xarray, yarray);

    new_x = [];
    new_y = [];
    new_east = [];
    new_north = [];
    new_vertical = [];
    for i in range(np.shape(X)[0]):
        for j in range(np.shape(X)[1]):
            new_x.append(X[i, j]);
            new_y.append(Y[i, j]);

    # Evaluate the linear or cubic interpolation function at new points
    for i in range(len(new_x)):
        new_east.append(f_east(new_x[i], new_y[i]))  # only want to give the functions one point at a time.
        new_north.append(f_north(new_x[i], new_y[i]));
        new_vertical.append(0.0);  # there's no vertical deformation in this field by construction.

    velref_e = f_east(reflon, reflat);
    velref_n = f_north(reflon, reflat);

    # Making the plotting variables
    narray_x = np.array(new_x);
    narray_y = np.array(new_y);
    narray_north = np.array(new_north)[:, 0];
    narray_east = np.array(new_east)[:, 0];

    # Northward and Eastward velocities in a subplot
    f, axarr = plt.subplots(1, 2, figsize=(12, 6));
    h1 = axarr[0].scatter(narray_x, narray_y, s=575, marker='s', c=narray_north, cmap='jet', edgecolors='face',
                          vmin=-30, vmax=30);
    axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0);
    axarr[0].quiver(reflon, reflat, velref_e, velref_n, color='purple', scale=500);
    my_plot_formatting(axarr[0], bounds, "North Velocity Interpolated from GPS", vel_tuple);
    h2 = axarr[1].scatter(narray_x, narray_y, s=575, marker='s', c=narray_east, cmap='jet', edgecolors='face', vmin=-30,
                          vmax=30);
    cbar = plt.colorbar(h1, ax=axarr[1]);
    cbar.set_label('mm/yr');
    axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0);
    axarr[1].quiver(reflon, reflat, velref_e, velref_n, color='purple', scale=500);
    my_plot_formatting(axarr[1], bounds, "East Velocity Interpolated from GPS", vel_tuple);
    plt.savefig("East_North_GPS_interp.png");
    plt.close();

    print("Reference e: %f" % velref_e)
    print("Reference n: %f" % velref_n)
    print("Are these reasonable values?");

    return;


def my_plot_formatting(ax, bounds, titlestring, vel_tuple):
    ax.set_xlim([bounds[0], bounds[1]]);
    ax.set_ylim([bounds[2], bounds[3]]);
    ax.quiver(vel_tuple.elon, vel_tuple.nlat, vel_tuple.e, vel_tuple.n, scale=500.0);
    ax.plot(vel_tuple.elon, vel_tuple.nlat, '.');
    ax.set_title(titlestring);
    return ax;
