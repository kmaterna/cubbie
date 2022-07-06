

"""
A set of functions that read and write vrt gdal grid files compatible with ISCE.
"""

import numpy as np
from . import isce_read_write
from Tectonic_Utils.read_write import netcdf_read_write


def isce_to_grd(isce_scalar_name, grdname):
    """ Convert an isce scalar file into a grdfile """
    data = isce_read_write.read_scalar_data(isce_scalar_name);
    firstLon, firstLat, dE, dN, xmin, xmax, nlon, nlat = isce_read_write.get_xmin_xmax_xinc_from_xml(isce_scalar_name+".xml");
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

