# A set of useful trigonometric and angle functions
# Especially good for working with look vectors
import numpy as np


def bearing_to_cartesian(heading):
    # Bearing or heading: CW from North
    # Cartesian: CCW from East
    return 90 - heading;


def cartesian_to_heading(cartesian_angle):
    # Take a cartesian angle in degrees (CCW from east)
    # Convert to heading angle in degrees (CW from north)
    return 90 - cartesian_angle;


def complement_angle(angle):
    return 90 - angle;


def cartesian_to_ccw_from_north(angle):
    return angle - 90;


def rotate_vector_by_angle(x0, y0, theta):
    # theta in degrees CCW from East, like mathematical definition.
    x_prime = x0 * np.cos(np.deg2rad(-theta)) - y0 * np.sin(np.deg2rad(-theta));
    y_prime = x0 * np.sin(np.deg2rad(-theta)) + y0 * np.cos(np.deg2rad(-theta));
    return x_prime, y_prime;


def normalize_look_vector(lkve, lkvn, lkvu):
    east_sq = np.square(lkve)
    north_sq = np.square(lkvn)
    up_sq = np.square(lkvu)
    sumarray = np.add(east_sq, north_sq)
    sumarray = np.add(sumarray, up_sq);
    magnitude = np.sqrt(sumarray);
    norm_lkve = np.divide(lkve, magnitude)
    norm_lkvn = np.divide(lkvn, magnitude)
    norm_lkvu = np.divide(lkvu, magnitude)
    return norm_lkve, norm_lkvn, norm_lkvu;


def calc_azimuth_incidence_from_lkv(lkve, lkvn, lkvu):
    # lkve, lkvn, lkvu describe vector from plane to ground
    # Convention: Azimuth angle measured from North in Anti-clockwise direction, in degrees, from ground to plane (I think)
    # Convention: Incidence angle measured from vertical at target (aka degrees from vertical at satellite) (always +ve), in degrees
    east_sq = np.square(lkve);
    north_sq = np.square(lkvn);
    sumarray = np.add(east_sq, north_sq);
    magnitude = np.sqrt(sumarray);
    azimuth_standard = np.arctan2(-lkvn, -lkve);
    azimuth_standard = np.rad2deg(azimuth_standard);
    azimuth = np.add(azimuth_standard, -90);

    incidence = np.arctan2(magnitude, -lkvu);
    incidence = np.rad2deg(incidence);
    return azimuth, incidence;


def look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u):
    """
    lkv_e, lkv_n, lkv_u are the components of the look vector from ground to satellite
    incidence angle is angle between look vector and vertical in degrees
    Flight angle is clockwise from north in degrees
    """
    unit_lkv = [lkv_e, lkv_n, lkv_u]
    unit_lkv = unit_lkv / np.linalg.norm(unit_lkv);
    vert_vector = [0, 0, 1]
    dotproduct = np.dot(unit_lkv, vert_vector);
    incidence_angle = np.rad2deg(np.arccos(dotproduct));

    lkv_horiz_angle = np.arctan2(lkv_n,
                                 lkv_e);  # the cartesian angle of the horizontal look vector (negative small # for DESC)
    heading_deg = cartesian_to_heading(np.rad2deg(lkv_horiz_angle));
    flight_angle = heading_deg + 90;  # satellite flies 90 degrees away from look vector direction
    return [flight_angle, incidence_angle];
