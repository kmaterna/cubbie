"""
Grid tools in Python
"""

import numpy as np

def clip_array_by_bbox(x, y, array1, bbox, verbose=True):
    """
    Clip an array into a bounding box that's smaller than the original grid. Like grdcut, but python API.

    :param x: 1d array, like lons
    :param y: 1d array, like lats
    :param array1: 2d array, shape is len(y) x len(x)
    :param bbox: [W, E, S, N]
    :param verbose: bool
    """
    original_bbox = (np.min(x), np.max(x), np.min(y), np.max(y));
    if bbox[0] < original_bbox[0] or bbox[1] > original_bbox[1] or \
            bbox[2] < original_bbox[2] or bbox[3] > original_bbox[3]:  # sanity check
        print("ERROR! invalid bounding box provided.")
    if verbose:
        print("  Original grid and bbox: ", np.shape(array1), original_bbox);
    find_w = np.argmin(np.abs(x-bbox[0]));
    find_e = np.argmin(np.abs(x-bbox[1]));
    find_s = np.argmin(np.abs(y-bbox[2]));
    find_n = np.argmin(np.abs(y-bbox[3]));
    new_x = np.array(x[find_w:find_e]);
    new_y = np.array(y[find_n:find_s]);
    new_array1 = array1[find_n:find_s, find_w:find_e].copy();  # saving off a subset of the array
    revised_bbox = (np.min(new_x), np.max(new_x), np.min(new_y), np.max(new_y));
    if verbose:
        print("  Revised grid and bbox: ", np.shape(new_array1), revised_bbox);
    return new_x, new_y, new_array1;


def mismatching_array_sizes(array_tuple):
    """
    Return 0 if the array sizes match.
    """
    num_arrays = len(array_tuple);
    for i in range(1, num_arrays):
        if np.shape(array_tuple[0]) != np.shape(array_tuple[i]):
            return 1;
    return 0;
