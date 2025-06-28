"""
Read and process wrapped JPL UAVSAR interferograms from UAVSAR website
"""

import numpy as np
import struct
from ..math_tools import phase_math
from tectonic_utils.geodesy import haversine, insar_vector_functions


def read_igram_data(data_file, ann_file, dtype='f', igram_type='ground'):
    """
    Data file for igrams is binary with real-complex float pairs
    Igram_type is ground or slant.

    :param data_file: string
    :param ann_file: string
    :param dtype: string, default 'f' for float
    :param igram_type: string, either 'ground' or 'slant'
    :returns: x, y, two 2d arrays, phase-amp
    """
    print("Reading %s-range file %s" % (igram_type, data_file))
    num_rows, num_cols = get_rows_cols(ann_file, igram_type)
    start_lon, start_lat, lon_inc, lat_inc = get_ground_range_corner_increment(ann_file)
    xarray = np.arange(start_lon, start_lon+lon_inc*num_cols, lon_inc)
    yarray = np.arange(start_lat, start_lat+lat_inc*num_rows, lat_inc)
    if len(yarray) > num_rows:
        yarray = yarray[0:num_rows]
    if len(xarray) > num_rows:
        xarray = xarray[0:num_cols]
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
    if np.shape(phase) != (len(yarray), len(xarray)):
        raise ValueError("shape of Phase doesn't match shape of axes: ", len(xarray), len(yarray), np.shape(phase))
    if np.shape(amp) != (len(yarray), len(xarray)):
        raise ValueError("shape of Amp doesn't match shape of axes: ", len(xarray), len(yarray), np.shape(amp))
    return xarray, yarray, phase, amp


def read_corr_data(data_file, ann_file, dtype='f', igram_type='ground'):
    """
    Read coherence into a 2D array from UAVSAR data

    :param data_file: filename, binary file of coherence from UAVSAR platform
    :param ann_file: string, filename of annotation file from UAVSAR platform
    :param dtype: default 'f'
    :param igram_type: string, default 'ground'
    :return: 2d array of coherence values
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


def get_rows_cols(ann_file, igram_type='ground'):
    """
    :param ann_file: string, filename
    :param igram_type: string, default 'ground'
    :return: two ints, the size of the data product in nrows, ncols
    """
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
    """
    :param ann_file: string, filename
    :return: four floats, the coordinates of the starting corner for the track, and the lon/lat spacing increments
    """
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


def get_near_range_far_range_incidence_angles(ann_file):
    """
    :param ann_file: string, filename
    :return: two floats, the incidence angle of the near-range and far-range pixels
    """
    near_range, far_range = 0, 0
    for line in open(ann_file):
        if 'Average Look Angle in Near Range' in line:
            near_range = float(line.split('=')[1].split()[0])
        if 'Average Look Angle in Far Range' in line:
            far_range = float(line.split('=')[1].split()[0])
    return near_range, far_range


def get_heading_angle(ann_file):
    """
    :param ann_file: string, filename
    :return: float, the angle of the airplane's heading, in degrees CW from north
    """
    heading_angle = 0
    for line in open(ann_file):
        if 'Peg Heading' in line:
            heading_angle = float(line.split('=')[1].split()[0])
    return heading_angle


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


def get_ground_range_right_corners(ann_file):
    """
    :param ann_file: string, filename
    :return: upper right lon, upper right lat, lower right lon, lower right lat
    """
    # upper right and lower right corners of the data in the track
    ur_lon, ur_lat, lr_lon, lr_lat = 0, 0, 0, 0
    for line in open(ann_file):
        if 'Approximate Upper Right Longitude' in line:
            ur_lon = float(line.split('=')[1].split()[0])
        if 'Approximate Upper Right Latitude' in line:
            ur_lat = float(line.split('=')[1].split()[0])
        if 'Approximate Lower Right Longitude' in line:
            lr_lon = float(line.split('=')[1].split()[0])
        if 'Approximate Lower Right Latitude' in line:
            lr_lat = float(line.split('=')[1].split()[0])
    return ur_lon, ur_lat, lr_lon, lr_lat


def get_four_ground_range_corners(ann_file):
    ul_lon, ul_lat, ll_lon, ll_lat = get_ground_range_left_corners(ann_file)
    ur_lon, ur_lat, lr_lon, lr_lat = get_ground_range_right_corners(ann_file)
    four_corners = [(ul_lon, ul_lat), (ll_lon, ll_lat), (lr_lon, lr_lat), (ur_lon, ur_lat), (ul_lon, ul_lat)]
    return four_corners


# ------------ JPL UAVSAR IGRAM FORMATS -------------- #
# A set of tools designed for handling of ground-range igrams
# from the JPL website for UAVSAR individual igram products

def read_los_rdr_geo_from_ground_ann_file(ann_file, x_axis, y_axis):
    """
    Make los.rdr.geo given .ann file from JPL website's UAVSAR interferograms and the ground-range sample points.
    x-axis and y-axis are the x and y arrays where los vectors will be extracted on a corresponding grid.
    The azimuth returned is the flight direction, CW from North, in degrees.

    :param ann_file: filename, string
    :param x_axis: 1d array, longitude values
    :param y_axis: 1d array, latitude values
    :return: 2d array of incidence angles, 2d array of azimuths, shape matching ((len(y), len(x))
    """
    near_angle, far_angle = get_near_range_far_range_incidence_angles(ann_file)
    heading = get_heading_angle(ann_file)  # in degrees CW from north
    heading_cart = insar_vector_functions.bearing_to_cartesian(heading)  # CCW from east
    print("Heading is %f degrees CW from north" % heading)
    print("Cartesian Heading is %f CCW from East" % heading_cart)

    # Get coords for upper and lower corners at start of the track, to compute length of across-track extent in km
    if 0 < heading < 180:
        far_start_lon, far_start_lat, near_start_lon, near_start_lat = get_ground_range_left_corners(ann_file)
    else:
        near_start_lon, near_start_lat, far_start_lon, far_start_lat = get_ground_range_right_corners(ann_file)
    cross_track_max = haversine.distance((near_start_lat, near_start_lon),
                                         (far_start_lat, far_start_lon))  # in km

    # This code is designed to make an azimuth angle consistent with ISCE format, for the rdr.los.geo file
    # Get azimuth angle for the pixels looking up to the airplane, consistent with ISCE I guess
    # My own documentation says CCW from north, even though that's really strange.
    azimuth = heading_cart - 90  # 90 degrees to the right of the airplane heading
    # (for the look vector from ground to plane)
    azimuth = insar_vector_functions.cartesian_to_ccw_from_north(azimuth)  # degrees CCW from North
    print("azimuth from ground to plane is: %s degrees CCW from north" % azimuth)

    [X, Y] = np.meshgrid(x_axis, y_axis)
    _grid_az_isce = azimuth * np.ones(np.shape(X))
    grid_az_flight = heading * np.ones(np.shape(X))

    # Compute the incidence angle for every pixel
    print("Computing incidence angles for all pixels")
    lat_flat, lon_flat = Y.reshape(-1), X.reshape(-1)
    nr_corner_lon = near_start_lon * np.ones(np.shape(lon_flat))  # corner at the near-range
    nr_corner_lat = near_start_lat * np.ones(np.shape(lat_flat))  # corner at the near-range
    target_coords = np.stack((lat_flat, lon_flat), axis=1)  # Nx2 array of coordinates
    corner_coords = np.stack((nr_corner_lat, nr_corner_lon), axis=1)  # Nx2 array of coordinates

    distance = haversine.distance_vectorized(target_coords, corner_coords)
    compass_bearing = haversine.calculate_initial_compass_bearing_vectorized(corner_coords, target_coords)  # CW from N
    theta = insar_vector_functions.bearing_to_cartesian(compass_bearing)  # angle of position vector, cartesian coords
    # heading_cartesian is the angle between east unit vector and the flight direction
    x0 = distance * np.cos(np.deg2rad(theta))
    y0 = distance * np.sin(np.deg2rad(theta))  # in the east-north coordinate system
    _, xtp = insar_vector_functions.rotate_vector_by_angle(x0, y0, heading_cart)
    xtp = xtp.reshape(np.shape(X))

    grid_inc = incidence_angle_trig(xtp, cross_track_max, near_angle, far_angle)

    # Finally, write the 2 bands for los.rdr.geo
    # isce_read_write.write_isce_unw(grid_inc, grid_az_isce, nx, ny, "FLOAT", 'los.rdr.geo')
    return grid_inc, grid_az_flight


def cross_track_pos(target_lon, target_lat, nearrange_lon, nearrange_lat, heading_cartesian):
    """
    Get cross-track position of point in a coordinate system centered at (nearrange_lon, nearrange_lat)
    with given heading of a plane and coordinates of one near-range point.
    Not used anymore for performance reasons.

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
