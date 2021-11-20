"""
A set of python scripts
What is the range/azimuth and row/col of a particular geographic coordinate?
This will use trans.dat, and some inputs of an example .grd file and geographic coords.
"""

import numpy as np
import subprocess
from Tectonic_Utils.read_write import netcdf_read_write


def get_nearest_row_col(example_grd, ra, az):
    [xdata, ydata, _] = netcdf_read_write.read_netcdf4(example_grd);
    col_idx = (np.abs(xdata - ra)).argmin()  # xdata is columns
    row_idx = (np.abs(ydata - az)).argmin()  # ydata is rows
    print(ydata[row_idx]);
    print(xdata[col_idx]);
    return [row_idx, col_idx];


def get_ra_from_ll(trans_dat, example_grd, lon, lat):
    """
    INPUTS: trans_dat: name of file
            example_grd: name of file
            lon: array or single value
            lat: array or single value
    RETURNS: a list of ra/az that matches the dimensions of the input lists.
    If the range and azimuth are outside of the range of the file, will return nan.
    """

    print("converting ll to ra")
    ll_temp = "geo_temp.txt"
    ra_temp = "ra_temp.txt"
    returnasarray = 1;

    # The case of a single float
    if type(lon) == np.float:
        returnasarray = 0;
        lon = [lon];
        lat = [lat];

    # Write the point we're going to convert
    ofile = open(ll_temp, 'w');
    for i in range(len(lon)):
        ofile.write("%f %f 0\n" % (lon[i], lat[i]));
    if len(lon) < 3:
        ofile.write("%f %f 0\n" % (lon[0] + 0.1, lat[0]));
        ofile.write("%f %f 0\n" % (lon[0] + 0.1, lat[0] + 0.1));
        ofile.write("%f %f 0\n" % (lon[0], lat[0] + 0.1));
    ofile.close();

    # Get the range and azimuth bouding box of this swath.
    ra_box = subprocess.check_output(['gmt', 'grdinfo', example_grd, '-I-']);
    ra_box = ra_box.decode('utf-8');
    ra_box = ra_box.split('\n')[0];
    range0 = float(ra_box.split('/')[0][2:]);
    range1 = float(ra_box.split('/')[1]);
    az0 = float(ra_box.split('/')[2]);
    az1 = float(ra_box.split('/')[3]);
    # print(range0, range1, az0, az1);

    # Here we implement the guts of proj_ll2ra_ascii.csh.
    # In Feb 2020, I needed to change the increment in the gmtinfo
    # I also needed to change the surface interpolation to -nn for a more stable interpolation near edges
    subprocess.call(['gmt', 'gmtconvert', trans_dat, '-o3,4,0', '-bi5d', '-bo3f', '>', 'llr'], shell=False);
    subprocess.call(['gmt', 'gmtconvert', trans_dat, '-o3,4,1', '-bi5d', '-bo3f', '>', 'lla'], shell=False);
    gmtrange = subprocess.check_output(['gmt', 'gmtinfo', ll_temp, '-I0.0007'], shell=False);
    gmtrange = gmtrange.decode('utf-8');
    gmtrange = gmtrange.split('\n')[0];
    errorcode = subprocess.call(['gmt', 'surface', 'llr', gmtrange, '-bi3f', '-I0.005', '-T.50', '-Gllr.grd', '-V'],
                                shell=False);
    if errorcode != 0:
        return [np.nan, np.nan];
    errorcode = subprocess.call(['gmt', 'surface', 'lla', gmtrange, '-bi3f', '-I0.005', '-T.50', '-Glla.grd', '-V'],
                                shell=False);
    if errorcode != 0:
        return [np.nan, np.nan];
    subprocess.call(['gmt', 'grdtrack', ll_temp, '-N', '-nn', '-Gllr.grd', '>', 'llpr'], shell=False);
    subprocess.call(['gmt', 'grdtrack', 'llpr', '-N', '-nn', '-Glla.grd', '>', 'llpra'], shell=False);
    subprocess.call("awk '{print $4,$5,$3}' < llpra > " + ra_temp, shell=True);
    subprocess.call(['rm', 'llr', 'lla', 'llpr', 'llpra', 'llr.grd', 'lla.grd'], shell=False);
    subprocess.call(['rm', 'gmt.history'], shell=False);

    # Here we read the results
    # We determine whether they're inside the range and azimuth box
    # We return the right values.
    ra_return = [];
    az_return = [];
    [ra, az, _] = np.loadtxt(ra_temp, unpack=True)
    for i in range(len(az)):
        if range0 <= ra[i] <= range1 and az0 <= az[i] <= az1:
            ra_return.append(ra[i]);
            az_return.append(az[i]);
        else:
            ra_return.append(np.nan);
            az_return.append(np.nan);

    if len(lon) <= 2:  # removing the points that were arbitrarily added
        ra_return = ra_return[0:len(lon)];
        az_return = az_return[0:len(lon)];

    if len(lon) == 1 and returnasarray == 1:  # returning a one-element array if that's what came in.
        ra_return = [ra_return[0]];
        az_return = [az_return[0]];
    if len(lon) == 1 and returnasarray == 0:  # returning a float if that's what came in.
        ra_return = ra_return[0];
        az_return = az_return[0];

    return [ra_return, az_return];


def get_ll_from_ra(trans_dat, ra, az):
    # Works on a single point
    print("converting ra to ll")
    ll_temp = "geo_temp.txt"
    ofile = open("ra_temp.txt", 'w');
    ofile.write("%f %f 0\n" % (ra, az));
    ofile.write("%f %f 0\n" % (ra + 100, az + 150));
    ofile.write("%f %f 0\n" % (ra + 150, az));
    ofile.close();
    subprocess.call(['rm', 'raln', 'raln.grd', 'ralt', 'ralt.grd'], shell=False);  # just to be clean.
    print("calling proj_ra2ll_ascii.csh " + trans_dat + ' ra_temp.txt ' + ll_temp);
    subprocess.call(['proj_ra2ll_ascii.csh', trans_dat, 'ra_temp.txt', ll_temp], shell=False);
    subprocess.call(['rm', 'gmt.history'], shell=False);
    [lon, lat, _] = np.loadtxt(ll_temp, unpack=True);
    lon = lon[0];
    lat = lat[0];
    return [lon, lat];


def get_ll_from_row_col(row, col, example_grd, trans_dat):
    [xdata, ydata, _] = netcdf_read_write.read_netcdf4(example_grd)
    ra = xdata[col];
    az = ydata[row];  # check this
    [lon, lat] = get_ll_from_ra(trans_dat, ra, az);
    return [lon, lat];


if __name__ == "__main__":
    trans_dat = "topo/trans.dat";
    example_grd = "stacking/unwrapped/2015157_2018177_unwrap.grd"
    lon, lat = [-116.572], [35.321];  # P617. Corresponding RA maybe equal to: ra = 10296.8986328, az = 5984.77251953
    [ra, az] = get_ra_from_ll(trans_dat, example_grd, lon, lat);
    [lon1, lat1] = get_ll_from_ra(trans_dat, ra, az);
    [row, col] = get_nearest_row_col(example_grd, ra, az);
    print(row, col)
