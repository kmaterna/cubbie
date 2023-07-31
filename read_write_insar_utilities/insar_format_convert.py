"""
A set of functions that read and write vrt gdal grid files compatible with ISCE.
"""

import numpy as np
from . import isce_read_write
from Tectonic_Utils.read_write import netcdf_read_write


def isce_to_grd(isce_name, grdname):
    """ Convert an isce scalar file into a grdfile """
    _, _, data = isce_read_write.read_scalar_data(isce_name);
    firstLon, firstLat, dE, dN, xmin, xmax, nlon, nlat = isce_read_write.get_xmin_xmax_xinc_from_xml(isce_name+".xml");
    (y, x) = np.shape(data);
    xarr = np.arange(firstLon, firstLon+x*dE, dE);
    yarr = np.arange(firstLat, firstLat+y*dN, dN);
    xarr = xarr[0:np.shape(data)[1]];
    yarr = yarr[0:np.shape(data)[0]];
    if dN < 0:
        data = np.flipud(data);
        yarr = np.flipud(yarr);
    netcdf_read_write.write_netcdf4(xarr, yarr, data, grdname);
    return;

def convert_intf_phase(infilename, outfilename):
    """
    Writes a single-band floating point number representing the phase from a CFLOAT32

    :param infilename: name of isce file, such as filt_fine.int
    :param outfilename: string
    """
    slc = isce_read_write.read_complex_data(infilename);
    phase = np.angle(slc);
    ny, nx = np.shape(phase);
    isce_read_write.write_isce_data(phase, nx, ny, "FLOAT", outfilename);
    return;

def extract_unw_phase(infilename, outfilename):
    """
    Write a single-band floating point number representing unwrapped phase from geo_unw.

    :param infilename: name of unw file, string
    :param outfilename: name of single-band file, string
    """
    _, _, unw = isce_read_write.read_scalar_data(infilename, band=2);  # reading second band
    ny, nx = np.shape(unw);
    isce_read_write.write_isce_data(unw, nx, ny, "FLOAT", outfilename);
    return;
