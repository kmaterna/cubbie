#!/usr/bin/env python
"""
This code performs cubic 2D interpolation on GPS velocities and then predicts a field of InSAR LOS velocities.
It assumes one look angle and flight vector (doesn't change incidence angle over the InSAR range)
It uses Scipy interpolate: cubic or linear
matplotlib.path: can replicate the matlab inpolygon() function.
Steps:
- Get interpolation points
- Interpolate based on function
- Project velocities into LOS
- Subtract reference point
- Report LOS gradient between points

Future tool-like versions of this should probably use PYGMT instead of oregon/ca border files.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as path
import collections
from scipy import interpolate
from gnss_timeseries_viewers.gps_tools import vel_functions
from . import los_projection_tools as los_proj
from . import file_io
from Tectonic_Utils.geodesy import haversine


param_object = collections.namedtuple("param_object", ['type_of_interp', 'asc_flight_angle', 'desc_flight_angle',
                                                       'inc_angle', 'bounds', 'inc', 'point1', 'point2',
                                                       'reference_point'])

# Define the global files that we are using
filedict = {'gps_input_file': "../../../../GEOPHYS_DATA/GPS_POS_DATA/PBO_Data/Velocities/NAM08_pbovelfile_feb2018.txt",
            'or_file': "../../../../Misc/Mapping_Resources/oregon_bdr",
            'ca_file': "../../../../Misc/Mapping_Resources/california_bdr"}


def do_interpolation():
    params = configure()
    [vel_tuple, ca_border, or_border] = inputs(filedict['gps_input_file'], filedict['or_file'], filedict['ca_file'])
    los_ascending, los_descending, eastnorthfield = compute(vel_tuple, params, ca_border, or_border)
    outputs(vel_tuple, params, ca_border, or_border, los_ascending, los_descending, eastnorthfield)


# -------------- CONFIG  --------------- # 
def configure():
    ascending_flight_angle = 360 - 14  # degrees from north (like strike)
    descending_flight_angle = 180 + 14  # degrees from north (like strike)
    incidence_angle = 30  # degrees from vertical (looking straight down is 0 degrees).
    inc = [0.08, 0.08]

    _Oregon_params = param_object(type_of_interp="cubic", asc_flight_angle=ascending_flight_angle,
                                  desc_flight_angle=descending_flight_angle, inc_angle=incidence_angle,
                                  bounds=[-125, -122, 41.5, 46.0], inc=inc, point1=[-123.0, 43.0],
                                  point2=[-123.0, 42.0], reference_point=[])

    _Mendocino_params = param_object(type_of_interp="cubic", asc_flight_angle=ascending_flight_angle,
                                     desc_flight_angle=descending_flight_angle, inc_angle=incidence_angle,
                                     bounds=[-125, -122, 39.0, 42.0], inc=inc, point1=[-124.0, 40.0],
                                     point2=[-124.0, 41.0], reference_point=[-123.4, 39.9])

    _Bay_Area_params = param_object(type_of_interp="cubic", asc_flight_angle=ascending_flight_angle,
                                    desc_flight_angle=descending_flight_angle, inc_angle=incidence_angle,
                                    bounds=[-124, -121.3, 36.8, 39.0], inc=inc, point1=[-122.3, 37.2],
                                    point2=[-121.5, 37.6], reference_point=[-121.5, 37.0])
    return _Mendocino_params  # Return the experiment you want.


# -------------- INPUTS  --------------- # 
def inputs(input_file, or_file, ca_file):
    """
    Input files. Should interpolate over the whole GPS field; don't chunk it smaller with clean_velfield.
    """
    vel_tuples = file_io.inputs_gps_pbo_like(input_file)
    ca_border = np.loadtxt(ca_file)
    or_border = np.loadtxt(or_file)
    return [vel_tuples, ca_border, or_border]


# ---------- COMPUTE --------------- #

def get_interp_points_within_grid(latlon_bounds, interval, optional_border_paths=None):
    """
    Establish new interpolation grid: set of points with chosen spacing (in degrees)
    Respects path boundaries if they are provided (to cancel oceans, etc.)
    :param latlon_bounds: list of [W, E, S, N] floats in degrees
    :param interval: list of two [xinc, yinc] in degrees
    :param optional_border_paths: list of paths, we restrict analysis to INSIDE of these paths
    :returns : [x, y] list of two 1-d arrays
    """
    xarray = np.arange(latlon_bounds[0], latlon_bounds[1], interval[0])
    yarray = np.arange(latlon_bounds[2], latlon_bounds[3], interval[1])
    [X, Y] = np.meshgrid(xarray, yarray)
    x_for_interp, y_for_interp = [], []

    if optional_border_paths is not None:
        # We remove interpolation points outside of CA and OR because they tend to explode.
        for k in range(len(optional_border_paths)):
            include_path = path.Path(optional_border_paths[k])
            for i in range(np.shape(X)[0]):
                for j in range(np.shape(X)[1]):
                    if include_path.contains_point((X[i, j], Y[i, j])) == 1:  # example: if point is in CA or OR
                        x_for_interp.append(X[i, j])
                        y_for_interp.append(Y[i, j])
    else:
        for i in range(np.shape(X)[0]):
            for j in range(np.shape(X)[1]):
                x_for_interp.append(X[i, j])
                y_for_interp.append(Y[i, j])
    return [x_for_interp, y_for_interp]


def compute(vel_tuple, params, ca_border, or_border):
    """
    Get 1D lists of points we're going to interpolate (either a grid or a vector of GPS points)
    """
    [x_for_interp, y_for_interp] = get_interp_points_within_grid(params.bounds, params.inc, [ca_border, or_border])

    elon_list = [item.elon for item in vel_tuple]
    nlat_list = [item.nlat for item in vel_tuple]
    e = [item.e for item in vel_tuple]
    n = [item.n for item in vel_tuple]
    # Interpolation functions that you can use on X-Y pairs
    f_east = interpolate.interp2d(elon_list, nlat_list, e, kind=params.type_of_interp)
    f_north = interpolate.interp2d(elon_list, nlat_list, n, kind=params.type_of_interp)

    # Evaluate the linear or cubic interpolation function at new points
    new_east, new_north, new_vertical = [], [], []
    los_ascending, los_descending, eastnorthfield = [], [], []
    for i in range(len(x_for_interp)):
        new_east.append(f_east(x_for_interp[i], y_for_interp[i]))  # only give functions one point at a time.
        new_north.append(f_north(x_for_interp[i], y_for_interp[i]))
        new_vertical.append(0.0)  # there's no vertical deformation in this field by construction.

    # Get the reference velocity in ENU
    velref_e, velref_n, velref_u = los_proj.get_point_enu_interp(params.reference_point, f_east, f_north)

    # Transform each field into LOS fields with respect to a given pixel
    ascending_LOS_array = los_proj.simple_project_ENU_to_LOS(new_east, new_north, new_vertical,
                                                             params.asc_flight_angle, params.inc_angle)
    ascending_LOS_reference = los_proj.simple_project_ENU_to_LOS(velref_e, velref_n, velref_u,
                                                                 params.asc_flight_angle, params.inc_angle)[0]

    # Same for Descending
    descending_LOS_array = los_proj.simple_project_ENU_to_LOS(new_east, new_north, new_vertical,
                                                              params.desc_flight_angle, params.inc_angle)
    descending_LOS_reference = los_proj.simple_project_ENU_to_LOS(velref_e, velref_n, velref_u,
                                                                  params.desc_flight_angle, params.inc_angle)[0]

    for i in range(len(x_for_interp)):
        item = vel_functions.Station_Vel(elon=x_for_interp[i], nlat=y_for_interp[i],
                                         e=ascending_LOS_array[i] - ascending_LOS_reference, n=0, u=0,
                                         se=0, sn=0, su=0, first_epoch=None, last_epoch=None, meas_type=None,
                                         name='', proccenter=None, refframe=None, subnetwork=None, survey=None)
        los_ascending.append(item)
        # Packing up an object for returning
        eastnorth_item = vel_functions.Station_Vel(elon=x_for_interp[i], nlat=y_for_interp[i],
                                                   e=new_east[i], n=new_north[i], u=0,
                                                   se=0, sn=0, su=0, first_epoch=None, last_epoch=None, meas_type=None,
                                                   name='', proccenter=None, refframe=None, subnetwork=None,
                                                   survey=None)
        eastnorthfield.append(eastnorth_item)
        # Packing up descending item
        item = vel_functions.Station_Vel(elon=x_for_interp[i], nlat=y_for_interp[i],
                                         e=descending_LOS_array[i] - descending_LOS_reference, n=0, u=0,
                                         se=0, sn=0, su=0, first_epoch=None, last_epoch=None, meas_type=None,
                                         name='', proccenter=None, refframe=None, subnetwork=None, survey=None)
        los_descending.append(item)

    # Here I want to evaluate gradients at two hard-coded points.
    e_def1, n_def1, u_def1 = los_proj.get_point_enu_interp(params.point1, f_east, f_north)
    e_def2, n_def2, u_def2 = los_proj.get_point_enu_interp(params.point2, f_east, f_north)
    ascending1 = los_proj.simple_project_ENU_to_LOS(e_def1, n_def1, u_def1, params.asc_flight_angle, params.inc_angle)
    ascending2 = los_proj.simple_project_ENU_to_LOS(e_def2, n_def2, u_def2, params.asc_flight_angle, params.inc_angle)
    descending1 = los_proj.simple_project_ENU_to_LOS(e_def1, n_def1, u_def1, params.desc_flight_angle,
                                                     params.inc_angle)
    descending2 = los_proj.simple_project_ENU_to_LOS(e_def2, n_def2, u_def2, params.desc_flight_angle,
                                                     params.inc_angle)
    evaluate_gradients(params.point1, params.point2, ascending1, ascending2, descending1, descending2)

    return los_ascending, los_descending, eastnorthfield


def evaluate_gradients(point1, point2, ascending1, ascending2, descending1, descending2):
    mydistance = haversine.distance([point1[1], point1[0]], [point2[1], point2[0]])
    print("Gradient in ASCENDING track from point1 to point2: %f mm/yr in %f km " % (np.abs(ascending1 - ascending2),
                                                                                     mydistance))
    print("Equal to: %f mm/yr per 100 km \n" % (100 * np.abs(ascending1 - ascending2) / mydistance))
    print("Gradient in DESCENDING track from point1 to point2: %f mm/yr in %f km " % (np.abs(descending1 - descending2),
                                                                                      mydistance))
    print("Equal to: %f mm/yr per 100 km \n" % (100 * np.abs(descending1 - descending2) / mydistance))
    return


# ---------- PLOTTING OUTPUTS --------------- # 
def outputs(vel_tuple, params, ca_border, or_border, los_ascending_velfield, los_descending_velfield, eastnorthfield):
    # Making the plotting variables
    narray_x = np.array([x.elon for x in eastnorthfield])
    narray_y = np.array([x.nlat for x in eastnorthfield])
    narray_north = np.array([x.n for x in eastnorthfield])[:, 0]
    narray_east = np.array([x.e for x in eastnorthfield])[:, 0]
    asc_dataset = np.array([x.e for x in los_ascending_velfield])[:, 0]
    des_dataset = np.array([x.e for x in los_descending_velfield])[:, 0]

    # Looking for the max los deformation values, to use for color bars.
    max_ascending_los, min_ascending_los = np.max(asc_dataset), np.min(asc_dataset)
    max_descending_los, min_descending_los = np.max(des_dataset), np.min(des_dataset)
    max_los_value = np.max([max_descending_los, max_ascending_los])
    min_los_value = np.min([min_descending_los, min_ascending_los])

    # Figure of interpolated velocities
    f1 = plt.figure()
    ax = f1.add_subplot(1, 1, 1)
    ax.plot(narray_x, narray_y, '.g')
    ax.quiver(narray_x, narray_y, narray_east, narray_north, color='red', scale=500.0)
    _f1 = my_plot_formatting(ax, params.bounds, ca_border, or_border, "Interpolated GPS Velocity Field", vel_tuple)
    plt.savefig("Interpolated_field.png")
    plt.close()

    # Northward and Eastward velocities in a subplot
    f, axarr = plt.subplots(1, 2, figsize=(12, 6))
    h1 = axarr[0].scatter(narray_x, narray_y, s=75, marker='s', c=narray_north, cmap='jet', edgecolors='face', vmin=-30,
                          vmax=30)
    axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0)
    my_plot_formatting(axarr[0], params.bounds, ca_border, or_border,
                       "North Velocity Interpolated from GPS", vel_tuple)
    _h2 = axarr[1].scatter(narray_x, narray_y, s=75, marker='s', c=narray_east, cmap='jet', edgecolors='face',
                           vmin=-30, vmax=30)
    cbar = plt.colorbar(h1, ax=axarr[1])
    cbar.set_label('mm/yr')
    axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0)
    my_plot_formatting(axarr[1], params.bounds, ca_border, or_border, "East Velocity Interpolated from GPS", vel_tuple)
    plt.savefig("East_North.png")
    plt.close()

    # ASCENDING AND DESCENDING VIEWING GEOMETRY
    f, axarr = plt.subplots(1, 2, figsize=(12, 6))
    axarr[0].scatter(narray_x, narray_y, s=75, marker='s', c=asc_dataset, cmap='jet', edgecolors='face',
                     vmin=min_los_value, vmax=max_los_value)
    axarr[0].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0)
    quiver_point = [min(params.bounds[0:2]) + 0.4, max(params.bounds[2:]) - 0.4]
    axarr[0].quiver(quiver_point[0], quiver_point[1],
                    np.cos(np.deg2rad(90 - params.asc_flight_angle)),
                    np.sin(np.deg2rad(90 - params.asc_flight_angle)), scale=10)
    axarr[0].quiver(quiver_point[0], quiver_point[1],
                    np.sin(np.deg2rad(90 - params.asc_flight_angle)),
                    -np.cos(np.deg2rad(90 - params.asc_flight_angle)), scale=10, color='red')
    axarr[0].text(quiver_point[0] + 0.05, quiver_point[1], "LOS", color='red', ha='right', va='top')
    axarr[0].plot(params.point1[0], params.point1[1], marker='s', color='black', markersize=5)
    axarr[0].plot(params.point2[0], params.point2[1], marker='s', color='black', markersize=5)
    axarr[0].plot(params.reference_point[0], params.reference_point[1], marker='s', color='black', markersize=5)
    my_plot_formatting(axarr[0], params.bounds, ca_border, or_border, "Ascending GPS LOS Velocity (interp)", vel_tuple)

    h1 = axarr[1].scatter(narray_x, narray_y, s=75, marker='s', c=des_dataset, cmap='jet', edgecolors='face',
                          vmin=min_los_value, vmax=max_los_value)
    cbar = plt.colorbar(h1, ax=axarr[1])
    cbar.set_label('mm/yr')
    quiver_point = [min(params.bounds[0:2]) + 0.4, max(params.bounds[2:]) - 0.4]
    axarr[1].quiver(quiver_point[0], quiver_point[1],
                    np.cos(np.deg2rad(90 - params.desc_flight_angle)),
                    np.sin(np.deg2rad(90 - params.desc_flight_angle)), scale=10)
    axarr[1].quiver(quiver_point[0], quiver_point[1],
                    np.sin(np.deg2rad(90 - params.desc_flight_angle)),
                    -np.cos(np.deg2rad(90 - params.desc_flight_angle)), scale=10, color='red')
    axarr[1].text(quiver_point[0] + 0.05, quiver_point[1], "LOS", color='red', ha='left', va='bottom')
    axarr[1].quiver(narray_x, narray_y, narray_east, narray_north, color='white', scale=500.0)
    axarr[1].plot(params.point1[0], params.point1[1], marker='s', color='black', markersize=5)
    axarr[1].plot(params.point2[0], params.point2[1], marker='s', color='black', markersize=5)
    axarr[1].plot(params.reference_point[0], params.reference_point[1], marker='s',
                  color='black', markersize=5)
    my_plot_formatting(axarr[1], params.bounds, ca_border, or_border,
                       "Descending GPS LOS Velocity (interp)", vel_tuple)

    # The sentinel descending scene.
    xboxes = [-123.260666, -124.37, -123.98, -122.816429, -123.260666]
    yboxes = [39.701538, 39.86, 41.6, 41.437393, 39.701538]
    axarr[1].plot(xboxes, yboxes, 'k')
    plt.savefig("LOS.png")
    return


def my_plot_formatting(ax, bounds, ca_border, or_border, titlestring, vel_tuple):
    narray_x = np.array([x.elon for x in vel_tuple])
    narray_y = np.array([x.nlat for x in vel_tuple])
    e = np.array([x.e for x in vel_tuple])
    n = np.array([x.n for x in vel_tuple])
    ax.set_xlim(bounds[0:2])
    ax.set_ylim(bounds[2:])
    ax.quiver(narray_x, narray_y, e, n, scale=500.0)
    ax.plot(ca_border[:, 0], ca_border[:, 1], 'k')
    ax.plot(or_border[:, 0], or_border[:, 1], 'k')
    ax.plot(narray_x, narray_y, '.')
    ax.set_title(titlestring)
    return ax


if __name__ == "__main__":
    do_interpolation()
