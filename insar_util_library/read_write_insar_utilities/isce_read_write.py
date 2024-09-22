"""
A set of functions that read and write vrt gdal grid files compatible with ISCE.
"""

import numpy as np
import matplotlib.pyplot as plt
import struct


# ----------- READING FUNCTIONS ------------- #

def read_isce_1d_arrays(filename):
    """
    Read only the axes of an ISCE file

    :param filename: string, name of isce file
    :return: 1d array for x-axis, 1d array for y-axis
    """
    firstlon, firstlat, dlon, dlat, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml')
    xarray, yarray = get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, nlon, nlat)
    return xarray, yarray


def read_complex_data(gdal_filename):
    """
    Read isce SLC data into a 2D array where each element is a complex number.

    :param gdal_filename: string, name of file
    :returns: 2d raster of complex values
    """
    from osgeo import gdal  # GDAL support for reading virtual files
    print("Reading file %s " % gdal_filename)
    ds = gdal.Open(gdal_filename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    _xmin, _xmax, _ymin, _ymax = get_xmin_xmax_xinc_from_geotransform(transform, slc)

    # put all zero values to nan
    slc = flush_zeros_to_nans(slc)

    return slc


def read_scalar_data(gdal_filename, band=1, flush_zeros=True):
    """
    Read an isce data file.
    band = 1 for most scalar fields, like coherence.
    band = 2 for some unwrapped phase files.

    :param gdal_filename: string, filename
    :param band: int representing the band of information, default is 1
    :param flush_zeros: default True
    :returns: x-axis 1d array, y-axis 1d array, data 2d raster array
    """
    from osgeo import gdal  # GDAL support for reading virtual files
    print("Reading file %s " % gdal_filename)
    if ".unw" in gdal_filename and ".unw." not in gdal_filename and band == 1:
        print("WARNING: We usually read band=2 for snaphu unwrapped files. Are you sure you want band 1 ????")
    ds = gdal.Open(gdal_filename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    _xmin, _xmax, _ymin, _ymax = get_xmin_xmax_xinc_from_geotransform(transform, data)
    xarray, yarray = read_isce_1d_arrays(gdal_filename)

    # put all zero values to nan
    if flush_zeros:
        data = flush_zeros_to_nans(data)

    return xarray, yarray, data


def read_phase_data(gdal_filename):
    """
    Start with a complex quantity, and return only the phase of that quantity.

    :param gdal_filename: string, name of file
    :returns: 2d raster data representing phase only
    """
    slc = read_complex_data(gdal_filename)
    phasearray = np.angle(slc)
    return phasearray


def read_amplitude_data(gdal_filename):
    """
    Start with a complex quantity, and return only the amplitude of that quantity.

    :param gdal_filename: string, name of file
    :returns: 2d raster data representing amplitude only
    """
    slc = read_complex_data(gdal_filename)
    amparray = np.absolute(slc)
    return amparray


def read_scalar_data_no_isce(filename, nx, ny):
    """
    Take float32 numbers from binary file into 2d array without using ISCE.

    :param filename: string
    :param nx: size of x-axis, int
    :param ny: size of y-axis, int
    :returns: 2d array of floats, of size (ny, nx)
    """
    final_shape = (ny, nx)
    num_data = nx * ny
    print("Reading file %s into %d x %d array" % (filename, ny, nx))
    with open(filename, 'rb') as f:
        raw_num = f.read()
    floats = np.array(struct.unpack('f' * num_data, raw_num))
    scalar_field = floats.reshape(final_shape)
    return scalar_field


def read_phase_data_no_isce(filename, nx, ny):
    """
    Read phase data from binary file into 2d array without using ISCE.

    :param filename: string
    :param nx: size of x-axis, int
    :param ny: size of y-axis, int
    :returns: 2d array of phase values, floats, of size (ny, nx)
    """
    final_shape = (ny, nx)
    num_data = nx * ny * 2
    print("Reading file %s into %d x %d array" % (filename, ny, nx))
    f = open(filename, 'rb')
    rawnum = f.read()
    floats = np.array(struct.unpack('f' * num_data, rawnum))
    f.close()
    real = floats[::2]
    imag = floats[1::2]
    phase = np.arctan2(imag, real)
    phase = phase.reshape(final_shape)
    return phase


def read_isce_unw_geo(filename):
    """
    Read isce unwrapped geocoded product, which has two datasets interleaved: amp and unwrapped phase
    Return x and y axes too, in lon/lat
    """
    firstlon, firstlat, dlon, dlat, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml')
    twox_data = read_scalar_data_no_isce(filename, nlon, nlat*2)   # separate the unw phase layer
    unw_data = twox_data[nlat:, :]   # unw_phase is the second layer
    (y, x) = np.shape(unw_data)
    xarray, yarray = get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, x, y)
    return xarray, yarray, unw_data


def read_isce_unw_geo_single(filename):
    """
    Read isce unwrapped geocoded single-band, which has one dataset

    :param filename: string
    :returns: 1D_xarray, 1D_yarray, 2D_data_array
    """
    firstlon, firstlat, dlon, dlat, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml')
    unw_data = read_scalar_data_no_isce(filename, nlon, nlat)
    (y, x) = np.shape(unw_data)
    xarray, yarray = get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, x, y)
    return xarray, yarray, unw_data


def read_isce_unw_geo_alternative(filename):
    """
    Read a custom isce unwrapped geocoded product, which has two copies of unwrapped phase
    Uses a format found in some unwrapped files
    Return x and y axes too, in lon/lat
    :param filename: string
    """
    firstlon, firstlat, dlon, dlat, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml')
    twox_data = read_scalar_data_no_isce(filename, nlon*2, nlat)   # separate the unw phase layer
    # unw_data = twox_data[:, 0:nlon]   # unw_phase is the second layer  # THIS COULD BE HAPPENING
    unw_data = twox_data[:, nlon:]   # unw_phase is the second layer  # THIS COULD BE HAPPENING

    (y, x) = np.shape(unw_data)
    xarray, yarray = get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, x, y)
    return xarray, yarray, unw_data


# ----------- WRITING FUNCTIONS ------------- #

def write_isce_data(data, nx, ny, dtype, filename, firstlat=None, firstlon=None, deltalon=None, deltalat=None,
                    xmin=None, xmax=None):
    """
    Write ISCE data into a single-band file with given filename with associated .vrt and .xml
    If DTYPE=="FLOAT": write scalar data (float32)
    IF DTYPE=="CFLOAT": write complex data (float32 + j*float32)
    """
    from isce.components import isceobj
    print("Writing data as file %s " % filename)
    out = isceobj.createImage()
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.setInterleavedScheme('BIP')  # 'BIP'/ 'BIL' / ‘BSQ’
    out.setAccessMode('READ')
    out.setDataType(dtype)
    if firstlon is not None:  # Special options that aren't usually used.
        out.setFirstLongitude(firstlon)
    if firstlat is not None:
        out.setFirstLatitude(firstlat)
    if deltalon is not None:
        out.setDeltaLongitude(deltalon)
    if deltalat is not None:
        out.setDeltaLatitude(deltalat)
    if xmin is not None:
        out.setXmin(xmin)
    if xmax is not None:
        out.setXmax(xmax)
    out.renderHdr()
    data.tofile(filename)  # write file out
    return


def write_isce_unw(data1, data2, nx, ny, dtype, filename, firstlat=None, firstlon=None, deltalon=None, deltalat=None,
                   xmin=None, xmax=None):
    """
    ISCE uses band=2 for unwrapped phase, .unw files.
    Write to float32
    """
    from isce.components import isceobj
    print("Writing data as file %s " % filename)
    out = isceobj.Image.createUnwImage()
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.imageType = 'unw'
    out.bands = 2
    out.scheme = "BIL"
    out.setAccessMode('read')
    if firstlon is not None:  # Special options that aren't usually used.
        out.setFirstLongitude(firstlon)
    if firstlat is not None:
        out.setFirstLatitude(firstlat)
    if deltalon is not None:
        out.setDeltaLongitude(deltalon)
    if deltalat is not None:
        out.setDeltaLatitude(deltalat)
    if xmin is not None:
        out.setXmin(xmin)
    if xmax is not None:
        out.setXmax(xmax)
    out.setDataType(dtype)
    out.renderHdr()
    data_to_file_2_bands(data1, data2, filename)  # dump the data into a binary file
    return


def data_to_file_2_bands(data1, data2, filename):
    data1 = np.float32(data1)  # we should be consistent about float types here.
    data2 = np.float32(data2)
    data = np.hstack((data1, data2))  # establishing two bands
    data.tofile(filename)
    return


def data_to_file_1_bands(data1, filename):
    data1 = np.float32(data1)  # we should be consistent about float types here.
    data1.tofile(filename)
    return


def plot_scalar_data(gdal_filename, band=1, title="", colormap='gray', aspect=1,
                     datamin=None, datamax=None, draw_colorbar=True, colorbar_orientation="horizontal", background=None,
                     outname=None):
    """
    Visualize scalar data within an ISCE file, assuming band 1 unless otherwise specified

    :param gdal_filename: string, filename of isce file
    :param band: int, default is 1
    :param title: string, default is ""
    :param colormap: string representing color map; default is 'gray'
    :param aspect: default 1
    :param datamin: float, minimum data value in color scale, default None
    :param datamax: float, maximum data value in color scale, default None
    :param draw_colorbar: bool, default True
    :param colorbar_orientation: string, default 'horizontal'
    :param background: default None
    :param outname: string, default None
    :return:
    """
    from osgeo import gdal
    ds = gdal.Open(gdal_filename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    # getting the min max of the axes
    # Note: this assumes that the transform is north-up
    # There are transform[2] and transform[4] for other projections (not implemented).
    xmin, xmax, ymin, ymax = get_xmin_xmax_xinc_from_geotransform(transform, data)

    # put all zero values to nan and do not plot nan
    if background is None:
        data = flush_zeros_to_nans(data)

    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, vmin=datamin, vmax=datamax, cmap=colormap, extent=(xmin, xmax, ymin, ymax))
    ax.set_title(title)
    if draw_colorbar is not None:
        _cbar = fig.colorbar(cax, orientation=colorbar_orientation)
    ax.set_aspect(aspect)
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname)
    return


def plot_complex_data(gdal_filename, title="", aspect=1, band=1, colormap='rainbow',
                      datamin=None, datamax=None, draw_colorbar=None, colorbar_orientation="horizontal", outname=None):
    from osgeo import gdal
    ds = gdal.Open(gdal_filename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    xmin, xmax, ymin, ymax = get_xmin_xmax_xinc_from_geotransform(transform, slc)

    # put all zero values to nan and do not plot nan
    slc = flush_zeros_to_nans(slc)

    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(1, 2, 1)
    cax1 = ax.imshow(np.abs(slc), vmin=datamin, vmax=datamax, cmap='gray', extent=(xmin, xmax, ymin, ymax))
    ax.set_title(title + " (amplitude)")
    if draw_colorbar is not None:
        _cbar1 = fig.colorbar(cax1, orientation=colorbar_orientation)
    ax.set_aspect(aspect)

    ax = fig.add_subplot(1, 2, 2)
    cax2 = ax.imshow(np.angle(slc), cmap=colormap, extent=(xmin, xmax, ymin, ymax))
    ax.set_title(title + " (phase [rad])")
    if draw_colorbar is not None:
        _cbar2 = fig.colorbar(cax2, orientation=colorbar_orientation)
    ax.set_aspect(aspect)
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname)
    return


# ----------- UTILITY FUNCTIONS ------------- #

def flush_zeros_to_nans(data_array):
    try:
        data_array[data_array == 0] = np.nan
    except:
        pass
    return data_array


def get_xmin_xmax_xinc_from_geotransform(transform, dataset):
    """
    Get min/max of transform axes
    """
    firstx = transform[0]
    firsty = transform[3]
    deltay = transform[5]
    deltax = transform[1]
    lastx = firstx + dataset.shape[1] * deltax
    lasty = firsty + dataset.shape[0] * deltay
    ymin = np.min([lasty, firsty])
    ymax = np.max([lasty, firsty])
    xmin = np.min([lastx, firstx])
    xmax = np.max([lastx, firstx])
    return xmin, xmax, ymin, ymax


def get_xmin_xmax_xinc_from_xml(xml_file):
    isce_xml = isce_xml_parser(xml_file)

    coord_lon = get_property(isce_xml, 'coordinate1')
    coord_lat = get_property(isce_xml, 'coordinate2')
    dlat = coord_lat['delta']
    dlon = coord_lon['delta']
    nlon = int(coord_lon['size'])
    nlat = int(coord_lat['size'])
    firstlat = coord_lat['startingvalue']
    firstlon = coord_lon['startingvalue']
    xmin = firstlon
    xmax = coord_lon['startingvalue'] + (nlon * coord_lon['delta'])
    return firstlon, firstlat, dlon, dlat, xmin, xmax, nlon, nlat


def get_xarray_yarray_from_xml(filename):
    """
    Get 1d axis-arrays from isce xml data

    :param filename: string, isce filename
    :returns: 1d x-array, 1d y-array
    """
    firstlon, firstlat, dlon, dlat, xmin, xmax, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename)
    xarray, yarray = get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, nlon, nlat)
    return xarray, yarray


def get_xarray_yarray_from_shape(firstlon, firstlat, dlon, dlat, x, y):
    """
    Build x and y arrays of latitude and longitude from the elements of the transform.

    :param firstlon: beginning of the x-axis
    :param firstlat: beginning of the y-axis
    :param dlon: the increment of the x-axis
    :param dlat: the increment of the y-axis
    :param x: int, size of the eventual x-array
    :param y: int, size of the eventual y-array
    :returns: 1d array representing x-axis, 1d array representing y-axis
    """
    xarray = np.arange(firstlon, firstlon+x*dlon, dlon)
    yarray = np.arange(firstlat, firstlat+y*dlat, dlat)
    return xarray[0:x], yarray[0:y]


def isce_xml_parser(filename):
    import xml.etree.ElementTree as et
    root = et.parse(filename).getroot()
    return root


def type_convert(value):
    for t in (float, int, str):
        try:
            return t(value)
        except ValueError:
            continue
    raise ValueError('Could not convert value')


def get_property(root, name):
    name = name.lower()
    values = {}

    for child in root.iter():
        child_name = child.get('name')
        if isinstance(child_name, str):
            child_name = child_name.lower()
        if child_name == name.lower():
            if child.tag == 'property':
                return type_convert(child.find('value').text)
            elif child.tag == 'component':
                values = {}
                for prop in child.iter('property'):
                    values[prop.get('name')] = type_convert(prop.find('value').text)
    return values
