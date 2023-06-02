"""
A set of functions that read and write vrt gdal grid files compatible with ISCE.
"""

import numpy as np
import matplotlib.pyplot as plt
import struct


# ----------- READING FUNCTIONS ------------- #


def read_complex_data(GDALfilename):
    """
    Read isce SLC data into a 2D array where each element is a complex number.
    """
    from osgeo import gdal  # GDAL support for reading virtual files
    print("Reading file %s " % GDALfilename);
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    _xmin, _xmax, _ymin, _ymax = get_xmin_xmax_xinc_from_geotransform(transform, slc);

    # put all zero values to nan
    try:
        slc[slc == 0] = np.nan
    except:
        pass

    return slc;


def read_scalar_data(GDALfilename, band=1, flush_zeros=True):
    """
    Read an isce data file.
    band = 1 for most scalar fields, like coherence.
    band = 2 for some unwrapped phase files.
    """
    from osgeo import gdal  # GDAL support for reading virtual files
    print("Reading file %s " % GDALfilename);
    if ".unw" in GDALfilename and ".unw." not in GDALfilename and band == 1:
        print("WARNING: We usually read band=2 for snaphu unwrapped files. Are you sure you want band 1 ????");
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    _xmin, _xmax, _ymin, _ymax = get_xmin_xmax_xinc_from_geotransform(transform, data);

    # put all zero values to nan
    if flush_zeros:
        try:
            data[data == 0] = np.nan
        except:
            pass

    return data;


def read_phase_data(GDALfilename):
    """
    Start with a complex quantity, and return only the phase of that quantity.
    """
    slc = read_complex_data(GDALfilename);
    phasearray = np.angle(slc);
    return phasearray;


def read_amplitude_data(GDALfilename):
    """
    Start with a complex quantity, and return only the amplitude of that quantity.
    """
    slc = read_complex_data(GDALfilename);
    amparray = np.absolute(slc);
    return amparray;


def read_scalar_data_no_isce(filename, nx, ny):
    """ Take float32 numbers from binary file into 2d array """
    final_shape = (ny, nx);
    num_data = nx * ny;
    print("Reading file %s into %d x %d array" % (filename, ny, nx));
    f = open(filename, 'rb')
    rawnum = f.read();
    f.close();
    floats = np.array(struct.unpack('f' * num_data, rawnum));
    scalar_field = floats.reshape(final_shape);
    return scalar_field;


def read_phase_data_no_isce(filename, nx, ny):
    final_shape = (ny, nx);
    num_data = nx * ny * 2;
    print("Reading file %s into %d x %d array" % (filename, ny, nx));
    f = open(filename, 'rb')
    rawnum = f.read();
    floats = np.array(struct.unpack('f' * num_data, rawnum))
    f.close();
    real = floats[::2];
    imag = floats[1::2];
    phase = np.arctan2(imag, real);
    phase = phase.reshape(final_shape);
    return phase;


def read_isce_unw_geo(filename):
    """
    Read isce unwrapped geocoded product, which has two datasets interleaved: amp and unwrapped phase
    Return x and y axes too, in lon/lat
    """
    firstLon, firstLat, dE, dN, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml');
    twox_data = read_scalar_data_no_isce(filename, nlon, nlat*2);   # separate the unw phase layer
    unw_data = twox_data[nlat:, :];   # unw_phase is the second layer
    (y, x) = np.shape(unw_data);
    xarray, yarray = get_xarray_yarray_from_shape(firstLon, firstLat, dE, dN, x, y);
    return xarray, yarray, unw_data;


def read_isce_unw_geo_single(filename):
    """
    Read isce unwrapped geocoded single-band, which has one dataset

    :returns: 1D_xarray, 1D_yarray, 2D_data_array
    """
    firstLon, firstLat, dE, dN, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml');
    unw_data = read_scalar_data_no_isce(filename, nlon, nlat);
    (y, x) = np.shape(unw_data);
    xarray, yarray = get_xarray_yarray_from_shape(firstLon, firstLat, dE, dN, x, y);
    return xarray, yarray, unw_data;


def read_isce_unw_geo_alternative(filename):
    """
    Read a custom isce unwrapped geocoded product, which has two copies of unwrapped phase
    Uses a format found in some unwrapped files
    Return x and y axes too, in lon/lat
    """
    firstLon, firstLat, dE, dN, _, _, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename+'.xml');
    twox_data = read_scalar_data_no_isce(filename, nlon*2, nlat);   # separate the unw phase layer
    unw_data = twox_data[:, 0:nlon];   # unw_phase is the second layer
    unw_data = twox_data[:, nlon:];   # unw_phase is the second layer

    (y, x) = np.shape(unw_data);
    xarray, yarray = get_xarray_yarray_from_shape(firstLon, firstLat, dE, dN, x, y);
    return xarray, yarray, unw_data;


# ----------- WRITING FUNCTIONS ------------- #

def write_isce_data(data, nx, ny, dtype, filename,
                    firstLat=None, firstLon=None, deltaLon=None, deltaLat=None, Xmin=None, Xmax=None):
    """
    Write ISCE data into a single-band file with given filename with associated .vrt and .xml
    If DTYPE=="FLOAT": write scalar data (float32)
    IF DTYPE=="CFLOAT": write complex data (float32 + j*float32)
    """
    from isce.components import isceobj
    print("Writing data as file %s " % filename);
    out = isceobj.createImage()
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.setInterleavedScheme('BIP')  # 'BIP'/ 'BIL' / ‘BSQ’
    out.setAccessMode('READ')
    out.setDataType(dtype)
    if firstLon is not None:  # Special options that aren't usually used. 
        out.setFirstLongitude(firstLon);
    if firstLat is not None:
        out.setFirstLatitude(firstLat);
    if deltaLon is not None:
        out.setDeltaLongitude(deltaLon);
    if deltaLat is not None:
        out.setDeltaLatitude(deltaLat);
    if Xmin is not None:
        out.setXmin(Xmin);
    if Xmax is not None:
        out.setXmax(Xmax);
    out.renderHdr()
    data.tofile(filename)  # write file out
    return


def write_isce_unw(data1, data2, nx, ny, dtype, filename,
                   firstLat=None, firstLon=None, deltaLon=None, deltaLat=None, Xmin=None, Xmax=None):
    """
    ISCE uses band=2 for unwrapped phase, .unw files.
    Write to float32
    """
    from isce.components import isceobj
    print("Writing data as file %s " % filename);
    out = isceobj.Image.createUnwImage()
    out.setFilename(filename)
    out.setWidth(nx)
    out.setLength(ny)
    out.imageType = 'unw'
    out.bands = 2
    out.scheme = "BIL"
    out.setAccessMode('read')
    if firstLon is not None:  # Special options that aren't usually used. 
        out.setFirstLongitude(firstLon);
    if firstLat is not None:
        out.setFirstLatitude(firstLat);
    if deltaLon is not None:
        out.setDeltaLongitude(deltaLon);
    if deltaLat is not None:
        out.setDeltaLatitude(deltaLat);
    if Xmin is not None:
        out.setXmin(Xmin);
    if Xmax is not None:
        out.setXmax(Xmax);
    out.setDataType(dtype)
    out.renderHdr()
    data_to_file_2_bands(data1, data2, filename);  # dump the data into a binary file
    return;


def data_to_file_2_bands(data1, data2, filename):
    data1 = np.float32(data1);  # we should be consistent about float types here.
    data2 = np.float32(data2);
    data = np.hstack((data1, data2));  # establishing two bands
    data.tofile(filename)
    return;


def data_to_file_1_bands(data1, filename):
    data1 = np.float32(data1);  # we should be consistent about float types here.
    data1.tofile(filename)
    return;


def plot_scalar_data(GDALfilename, band=1, title="", colormap='gray', aspect=1,
                     datamin=None, datamax=None, draw_colorbar=True, colorbar_orientation="horizontal", background=None,
                     outname=None):
    from osgeo import gdal
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    data = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    # getting the min max of the axes
    # Note: this assumes that the transform is north-up
    # There are transform[2] and transform[4] for other projections (not implemented).
    xmin, xmax, ymin, ymax = get_xmin_xmax_xinc_from_geotransform(transform, data);

    # put all zero values to nan and do not plot nan
    if background is None:
        try:
            data[data == 0] = np.nan
        except:
            pass

    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(111)
    cax = ax.imshow(data, vmin=datamin, vmax=datamax, cmap=colormap, extent=[xmin, xmax, ymin, ymax])
    ax.set_title(title)
    if draw_colorbar is not None:
        cbar = fig.colorbar(cax, orientation=colorbar_orientation)
    ax.set_aspect(aspect)
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname);
    return;


def plot_complex_data(GDALfilename, title="", aspect=1, band=1, colormap='rainbow',
                      datamin=None, datamax=None, draw_colorbar=None, colorbar_orientation="horizontal", outname=None):
    from osgeo import gdal
    ds = gdal.Open(GDALfilename, gdal.GA_ReadOnly)
    slc = ds.GetRasterBand(band).ReadAsArray()
    transform = ds.GetGeoTransform()
    ds = None

    xmin, xmax, ymin, ymax = get_xmin_xmax_xinc_from_geotransform(transform, slc);

    # put all zero values to nan and do not plot nan
    try:
        slc[slc == 0] = np.nan
    except:
        pass

    fig = plt.figure(figsize=(18, 16))
    ax = fig.add_subplot(1, 2, 1)
    cax1 = ax.imshow(np.abs(slc), vmin=datamin, vmax=datamax, cmap='gray', extent=[xmin, xmax, ymin, ymax])
    ax.set_title(title + " (amplitude)")
    if draw_colorbar is not None:
        _cbar1 = fig.colorbar(cax1, orientation=colorbar_orientation)
    ax.set_aspect(aspect)

    ax = fig.add_subplot(1, 2, 2)
    cax2 = ax.imshow(np.angle(slc), cmap=colormap, extent=[xmin, xmax, ymin, ymax])
    ax.set_title(title + " (phase [rad])")
    if draw_colorbar is not None:
        _cbar2 = fig.colorbar(cax2, orientation=colorbar_orientation)
    ax.set_aspect(aspect)
    if outname is None:
        plt.show()
    else:
        fig.savefig(outname);

    return;


# ----------- UTILITY FUNCTIONS ------------- #

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
    return xmin, xmax, ymin, ymax;


def get_xmin_xmax_xinc_from_xml(xml_file):
    isce_xml = ISCEXMLParser(xml_file)

    coord_lon = getProperty(isce_xml, 'coordinate1')
    coord_lat = getProperty(isce_xml, 'coordinate2')
    dN = coord_lat['delta']
    dE = coord_lon['delta']
    nlon = int(coord_lon['size'])
    nlat = int(coord_lat['size'])
    firstLat = coord_lat['startingvalue']
    firstLon = coord_lon['startingvalue']
    xmin = firstLon;
    xmax = coord_lon['startingvalue'] + (nlon * coord_lon['delta'])
    return firstLon, firstLat, dE, dN, xmin, xmax, nlon, nlat;


def get_xarray_yarray_from_xml(filename):
    """
    Another function to get arrays from isce xml data
    """
    firstLon, firstLat, dE, dN, xmin, xmax, nlon, nlat = get_xmin_xmax_xinc_from_xml(filename);
    xarray, yarray = get_xarray_yarray_from_shape(firstLon, firstLat, dE, dN, nlon, nlat);
    return xarray, yarray;


def get_xarray_yarray_from_shape(firstLon, firstLat, dE, dN, x, y):
    """
    Building x and y arrays of latitude and longitude
    """
    xarray = np.arange(firstLon, firstLon+x*dE, dE);
    yarray = np.arange(firstLat, firstLat+y*dN, dN);
    return xarray[0:x], yarray[0:y];


def ISCEXMLParser(filename):
    import xml.etree.ElementTree as ET
    root = ET.parse(filename).getroot()
    return root;


def type_convert(value):
    for t in (float, int, str):
        try:
            return t(value)
        except ValueError:
            continue
    raise ValueError('Could not convert value')


def getProperty(root, name):
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
