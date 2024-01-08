"""
July 2020
Read and process a wrapped JPL UAVSAR interferogram from the UAVSAR website
"""

import numpy as np
import struct


def read_igram_data(data_file, ann_file, dtype='f', igram_type='ground', return_type='phase_amp'):
    """
    Data file for igrams is binary with real-complex float pairs
    Igram_type is ground or slant
    return_type is phase_amp or real_imag
    """
    print("Reading %s-range file %s" % (igram_type, data_file))
    num_rows, num_cols = get_rows_cols(ann_file, igram_type)
    f = open(data_file, 'rb')
    final_shape = (num_rows, num_cols)
    num_data = final_shape[0] * final_shape[1] * 2  # 2 for real/complex
    rawnum = f.read()
    f.close()
    floats = np.array(struct.unpack(dtype * num_data, rawnum))
    real = floats[::2]
    imag = floats[1::2]
    phase = np.arctan2(imag, real)
    amp = np.sqrt(np.multiply(imag, imag) + np.multiply(real, real))

    phase = phase.reshape(final_shape)
    amp = amp.reshape(final_shape)
    real = real.reshape(final_shape)
    imag = imag.reshape(final_shape)
    if return_type == "real_imag":
        return real, imag
    else:
        return phase, amp


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


def get_nearrange_farrange_heading_angles(ann_file):
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
