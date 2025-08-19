"""
LiCSBAS files come in two formats.
The first is standard geotiff -- use gdal_read_write_convert.read_geotiff
The second is a binary data dump.  Use the following functions.
"""

import numpy as np
from . import isce_read_write


def read_licsbas_file(bin_filename, xlen, ylen, bbox, verbose=False):
    """
    :param bin_filename: string, filename
    :param xlen: integer, size of x-array, likely read this from the LICSBAS logs
    :param ylen: integer, size of y-array, likely read this from the LICSBAS logs
    :param bbox: list of 4 integers, [W, E, S, N] in degrees, likely read this from the LICSBAS logs
    :param verbose: default False
    :return: xarray (1d), yarray (1d), zdata (2d)
    """

    data = np.fromfile(bin_filename, dtype=np.float32)
    data = data.reshape((ylen, xlen))
    yinc = (bbox[3] - bbox[2]) / ylen
    xinc = (bbox[1] - bbox[0]) / xlen

    if verbose:
        print("Reading file %s" % bin_filename)
        print("  Resulting size: (%d, %d)" % (ylen, xlen))
        print("  Increments: %.5f, %.5f" % (xinc, yinc))

    xarr, yarr = isce_read_write.get_xarray_yarray_from_shape(bbox[0], bbox[2], xinc, -yinc, xlen, ylen)
    return xarr, yarr, data
