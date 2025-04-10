"""
Read and process a wrapped JPL UAVSAR interferogram from the UAVSAR website
"""

import numpy as np
import struct
from ..math_tools import phase_math
from Tectonic_Utils.geodesy import haversine, insar_vector_functions


def read_igram_data(data_file, ann_file, dtype='f', igram_type='ground'):
    """
    Data file for igrams is binary with real-complex float pairs
    Igram_type is ground or slant.

    :param data_file: string
    :param ann_file: string
    :param dtype: string, default 'f'
    :param igram_type: string, either 'ground' or 'slant'
    :returns: two 2d arrays, phase-amp
    """
    print("Reading %s-range file %s" % (igram_type, data_file))
    num_rows, num_cols = get_rows_cols(ann_file, igram_type)
    start_lon, start_lat, lon_inc, lat_inc = get_ground_range_corner_increment(ann_file)
    xarray = np.arange(start_lon, start_lon+lon_inc*num_cols, lon_inc)
    yarray = np.arange(start_lat, start_lat+lat_inc*num_rows, lat_inc)
    f = open(data_file, 'rb')
    final_shape = (num_rows, num_cols)
    num_data = final_shape[0] * final_shape[1] * 2  # 2 for real/complex
    rawnum = f.read()
    f.close()
    floats = np.array(struct.unpack(dtype * num_data, rawnum))
    real, imag = floats[::2], floats[1::2]
    real = real.reshape(final_shape)
    imag = imag.reshape(final_shape)
    phase, amp = phase_math.real_imag2phase_amp(real, imag)
    return xarray, yarray, phase, amp


def read_corr_data(data_file, ann_file, dtype='f', igram_type='ground'):
    """
    Data file is not a regular netcdf grd file
    dtype float works for corr, unwrapped, etc.
    igram_type is ground or slant
    """
    print("Reading %s-range file %s" % (igram_type, data_file))
    num_rows, num_cols = get_rows_cols(ann_file, igram_type)
    f = open(data_file, 'rb')
    final_shape = (num_rows, num_cols)
    num_data = final_shape[0] * final_shape[1]
    rawnum = f.read()
    f.close()
    floats = np.array(struct.unpack(dtype * num_data, rawnum))
    data = floats.reshape(final_shape)
    return data


def get_rows_cols(ann_file, igram_type):
    num_rows, num_cols = 0, 0
    for line in open(ann_file):
        if igram_type == 'ground':
            if 'Ground Range Data Latitude Lines' in line:
                num_rows = int(line.split('=')[1])
            if 'Ground Range Data Longitude Samples' in line:
                num_cols = int(line.split('=')[1])
        elif igram_type == 'slant':
            if 'Slant Range Data Azimuth Lines' in line:
                num_rows = int(line.split('=')[1])
            if 'Slant Range Data Range Samples' in line:
                num_cols = int(line.split('=')[1])
        else:
            print("Error! Igram type not recognized")
    return num_rows, num_cols


def get_ground_range_corner_increment(ann_file):
    start_lon, start_lat, lon_inc, lat_inc = 0, 0, 0, 0
    for line in open(ann_file):
        if 'Ground Range Data Starting Latitude' in line:
            start_lat = float(line.split('=')[1].split()[0])
        if 'Ground Range Data Starting Longitude' in line:
            start_lon = float(line.split('=')[1].split()[0])
        if 'Ground Range Data Latitude Spacing' in line:
            lat_inc = float(line.split('=')[1].split()[0])
        if 'Ground Range Data Longitude Spacing' in line:
            lon_inc = float(line.split('=')[1].split()[0])
    return start_lon, start_lat, lon_inc, lat_inc


def get_near_range_far_range_heading_angles(ann_file):
    near_range, far_range, heading_angle = 0, 0, 0
    for line in open(ann_file):
        if 'Average Look Angle in Near Range' in line:
            near_range = float(line.split('=')[1].split()[0])
        if 'Average Look Angle in Far Range' in line:
            far_range = float(line.split('=')[1].split()[0])
        if 'Peg Heading' in line:
            heading_angle = float(line.split('=')[1].split()[0])
    return near_range, far_range, heading_angle


def get_ground_range_left_corners(ann_file):
    """
    :param ann_file: string, filename
    :return: upper left lon, upper left lat, lower left lon, lower left lat
    """
    # upper left and lower left corners of the data in the track
    ul_lon, ul_lat, ll_lon, ll_lat = 0, 0, 0, 0
    for line in open(ann_file):
        if 'Approximate Upper Left Longitude' in line:
            ul_lon = float(line.split('=')[1].split()[0])
        if 'Approximate Upper Left Latitude' in line:
            ul_lat = float(line.split('=')[1].split()[0])
        if 'Approximate Lower Left Longitude' in line:
            ll_lon = float(line.split('=')[1].split()[0])
        if 'Approximate Lower Left Latitude' in line:
            ll_lat = float(line.split('=')[1].split()[0])
    return ul_lon, ul_lat, ll_lon, ll_lat


# ------------ JPL UAVSAR IGRAM FORMATS -------------- #
# A set of tools designed for handling of ground-range igrams
# from the JPL website for UAVSAR individual igram products

def read_los_rdr_geo_from_ground_ann_file(ann_file, x_axis, y_axis):
    """
    Make los.rdr.geo given .ann file from JPL website's UAVSAR interferograms and the ground-range sample points.
    x-axis and y-axis are the x and y arrays where los vectors will be extracted on a corresponding grid.

    :param ann_file: filename, string
    :param x_axis: 1d array
    :param y_axis: 1d array
    :return: 2d array of incidence angles, 2d array of azimuths
    """
    near_angle, far_angle, heading = get_near_range_far_range_heading_angles(ann_file)
    heading_cartesian = insar_vector_functions.bearing_to_cartesian(heading)  # CCW from east
    print("Heading is %f degrees CW from north" % heading)
    print("Cartesian Heading is %f" % heading_cartesian)
    # Get the upper and lower left corners, so we can compute the length of the across-track extent in km
    ul_lon, ul_lat, ll_lon, ll_lat = get_ground_range_left_corners(ann_file)

    cross_track_max = haversine.distance((ll_lat, ll_lon), (ul_lat, ul_lon))  # in km

    # Get the azimuth angle for the pixels looking up to the airplane
    # My own documentation says CCW from north, even though that's really strange.
    azimuth = heading_cartesian - 90  # 90 degrees to the right of the airplane heading
    # (for the look vector from ground to plane)
    azimuth = insar_vector_functions.cartesian_to_ccw_from_north(azimuth)  # degrees CCW from North
    print("azimuth from ground to plane is:", azimuth)

    [X, Y] = np.meshgrid(x_axis, y_axis)
    (ny, nx) = np.shape(X)
    grid_az = azimuth * np.ones(np.shape(X))
    grid_inc = np.zeros(np.shape(X))
    xtp = np.zeros(np.shape(X))

    print("Computing incidence angles for all pixels")
    for i in range(ny):
        for j in range(nx):
            xtp[i, j] = cross_track_pos(X[i, j], Y[i, j], ll_lon, ll_lat, heading_cartesian)  # DIFFERENT FOR ASC AND DESC

    grid_inc = incidence_angle_trig(xtp, cross_track_max, near_angle, far_angle)

    # Finally, write the 2 bands for los.rdr.geo
    # isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", 'los.rdr.geo')
    return grid_inc, grid_az


def cross_track_pos(target_lon, target_lat, nearrange_lon, nearrange_lat, heading_cartesian):
    """
    Get cross-track position of point in a coordinate system centered at (nearrange_lon, nearrange_lat)
    with given heading of a plane and coordinates of one near-range point

    :param target_lon: float, longitude
    :param target_lat: float, latitude
    :param nearrange_lon: float, longitude
    :param nearrange_lat: float, latitude
    :param heading_cartesian: float, heading in a cartesian system from the x-axis
    :return: value across the track, in km
    """
    tupleA = (nearrange_lat, nearrange_lon)
    tupleB = (target_lat, target_lon)
    distance = haversine.distance(tupleB, tupleA)
    compass_bearing = haversine.calculate_initial_compass_bearing(tupleA, tupleB)  # this comes CW from north
    theta = insar_vector_functions.bearing_to_cartesian(compass_bearing)  # angle of position vector, cartesian coords
    # heading_cartesian is the angle between east unit vector and the flight direction
    x0 = distance * np.cos(np.deg2rad(theta))
    y0 = distance * np.sin(np.deg2rad(theta))  # in the east-north coordinate system
    x_prime, y_prime = insar_vector_functions.rotate_vector_by_angle(x0, y0, heading_cartesian)
    return y_prime


def incidence_angle_trig(xtp, cross_track_max, near_inc_angle, far_inc_angle):
    """
    Using the incidence angles (to the vertical) at the upper and lower corners of the track,
    what's the incidence angle at some location in between (xtp=cross-track-position)?
    near_angle is the incidence angle between the viewing geometry and the vertical at the near-range.
    nearcomp is the complement of that angle.
    This function is kind of like linear interpolation, but a little bit curved
    It solves an equation I derived on paper from the two near-range and far-range triangles in July 2020
    Now operates on numpy arrays.

    :param xtp: the cross-track position of the target point, float or 2d array
    :param cross_track_max: the cross-track position of the farthest range point, float
    :param near_inc_angle: the incidence angle at the near-range, float
    :param far_inc_angle: the incidence angle at the far-range, float
    :return: same data type as xtp, float or 2d array
    """
    nearcomp = np.deg2rad(np.subtract(90, near_inc_angle))
    farcomp = np.deg2rad(np.subtract(90, far_inc_angle))  # angles from ground to satellite
    h = (np.tan(nearcomp) * np.tan(farcomp) * cross_track_max) / (np.tan(nearcomp) - np.tan(farcomp))
    angle_to_horizontal = np.rad2deg(np.arctan(h / (xtp + (h / np.tan(nearcomp)))))
    return np.subtract(90, angle_to_horizontal)
