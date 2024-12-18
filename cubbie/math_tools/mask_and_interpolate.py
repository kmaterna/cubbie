import numpy as np
import scipy.interpolate


def make_coherence_mask(cor, threshold):
    """Build a mask: 1 if above coherence threshold, nan if below coherence threshold.

    :param cor: 2D array of raster data, between 0 and 1, holding correlation data
    :param threshold: float, used for the cutoff
    :returns: mask, 2D raster, either ones or nans
    """
    print("Making coherence mask.")
    mask = np.ones(np.shape(cor))
    mask[np.isnan(cor)] = np.nan
    mask[cor < threshold] = np.nan
    return mask


def apply_coherence_mask(data, mask, is_complex=False, is_float32=False, mask_value=np.nan):
    """
    A future version of this function should probably check the type of the input data
    and return the same type that came in.

    :param data: 2d raster
    :param mask: 2d raster
    :param is_complex: bool, default False
    :param is_float32: bool, default False
    :param mask_value: default is np.nan
    :returns: 2d raster with mask applied, in the same format as the original data
    """
    if np.shape(data) != np.shape(mask):
        raise ValueError("Error! Shape of data ("+str(np.shape(data))+") and shape of mask (" +
                         str(np.shape(mask))+") do not match")
    if is_float32:
        if is_complex == 1:
            masked = np.complex32(np.multiply(data, mask))
        else:
            masked = np.float32(np.multiply(data, mask))
            masked[np.isnan(masked)] = mask_value
    else:
        if is_complex == 1:
            masked = np.complex64(np.multiply(data, mask))
        else:
            masked = np.float64(np.multiply(data, mask))
            masked[np.isnan(masked)] = mask_value
    return masked


def interpolate_2d(data_array, is_complex=False):
    """
    A function to perform scipy linear interpolation scheme on 2d raster data, for custom workflows.
    Removes nans and replaces them with an interpolated version.

    :param data_array: 2d raster data
    :param is_complex: bool, default is False
    :returns: an interpolated 2d raster
    """
    print("Performing 2d interpolation")
    if is_complex:
        data_array = np.angle(data_array)

    ymax, xmax = np.shape(data_array)
    yarray, xarray = range(ymax), range(xmax)

    x_interps, y_interps, z_interps = [], [], []
    xy_interps, xy_targets = [], []

    # Time to get rid of nan's.
    for i in range(len(yarray)):
        for j in range(len(xarray)):
            # Get the real values for use in interpolating
            if not np.isnan(data_array[i][j]):
                x_interps.append(xarray[j])
                y_interps.append(yarray[i])
                xy_interps.append([xarray[j], yarray[i]])
                z_interps.append(data_array[i][j])
            # Collect the points where we interpolate
            else:
                xy_targets.append([xarray[j], yarray[i]])

    # f = scipy.interpolate.interp2d(y_interps, x_interps, z_interps)
    z_targets = scipy.interpolate.griddata(xy_interps, z_interps, xy_targets, method='linear', fill_value=1)

    # Fill in gaps using the interpolated value of z
    smoothdata = np.copy(data_array)
    for i in range(len(xy_targets)):
        idxx = xarray.index(xy_targets[i][0])
        idxy = yarray.index(xy_targets[i][1])
        if is_complex:
            r = 1
            smoothdata[idxy][idxx] = np.cmath.rect(r, z_targets[i])
        else:
            smoothdata[idxy][idxx] = z_targets[i]

    return smoothdata
