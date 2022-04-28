"""
A few functions that help project into and out of LOS
Mostly trig
Example: Dlos = [U_n sin(phi) - U_e cos(phi)]*sin(lamda) + U_u cos(lamda)
We use Station_Vels hold GPS velocity fields and their LOS-projected versions
In the LOS-projected case, we use 'e' as LOS velocity and other columns are zeros
"""

import numpy as np
from Tectonic_Utils.geodesy import insar_vector_functions
from GNSS_TimeSeries_Viewers.gps_tools import gps_io_functions


def closest_index(lst, K):
    # Given list lst, and target K, which element is closest?
    index = min(range(len(lst)), key=lambda i: abs(lst[i] - K));
    deg_distance = lst[index] - K;
    return index, deg_distance;


def simple_project_ENU_to_LOS(U_e, U_n, U_u, flight_angle, incidence_angle):
    """
    Dlos = [U_n sin(phi) - U_e cos(phi)]*sin(lamda) + U_u cos(lamda), Fialko 2001
    [U_e, U_n, U_u] are the east, north, and up components of the deformation.
    phi   = azimuth of satellite heading vector, positive clockwise from north.
    lamda = local incidence angle at the reflector (usually angle from the vertical).
    Works for single values and for 1D arrays of all arguments
    """

    if np.size(U_e) > 1:  # processing a 1D array of values
        d_los = np.zeros(np.shape(U_e));
        if np.size(flight_angle) == 1:
            for i in range(len(U_e)):
                d_los[i] = insar_vector_functions.def3D_into_LOS(U_e[i], U_n[i], U_u[i], flight_angle, incidence_angle);
        else:
            fa_array = [flight_angle for _x in U_e];  # bumping up the flight and incidence angles into 1D arrays
            ia_array = [incidence_angle for _x in U_e];
            for i in range(len(U_e)):
                d_los[i] = insar_vector_functions.def3D_into_LOS(U_e[i], U_n[i], U_u[i], fa_array[i], ia_array[i]);

    else:  # processing a single value
        d_los = insar_vector_functions.def3D_into_LOS(U_e, U_n, U_u, flight_angle, incidence_angle);
    return d_los


def get_point_enu_interp(reference_point, f_east=None, f_north=None, f_up=None):
    """
    Useful for reference points, for example
    reference point is [lon, lat]
    Interpolated functions are required if the reference isn't located at a gps station.
    """
    velref_e = f_east(reference_point[0], reference_point[1]);
    velref_n = f_north(reference_point[0], reference_point[1]);
    if f_up is None:
        velref_u = 0;
    else:
        velref_u = f_up(reference_point[0], reference_point[1]);
    return velref_e, velref_n, velref_u;


def get_point_enu_veltuple(vel_tuple, reference_point_coords=None, reference_pt_name=None, zero_vertical=False):
    """
    If the reference is co-located with a GPS station, we find it within the tuple.
    We either use name or lat/lon
    """
    velref_e, velref_n, velref_u, reflon, reflat = 0, 0, 0, 0, 0;
    if reference_pt_name is not None:
        for station_vel in vel_tuple:
            if station_vel.name == reference_pt_name:
                velref_e = station_vel.e;
                velref_n = station_vel.n;
                if zero_vertical:
                    velref_u = 0;
                else:
                    velref_u = station_vel.u;
                reflon = station_vel.elon;
                reflat = station_vel.nlat;
        return velref_e, velref_n, velref_u, reflon, reflat;
    else:
        for station_vel in vel_tuple:
            # should I make these tolerances?
            if station_vel.elon == reference_point_coords[0] and station_vel.nlat == reference_point_coords[1]:
                velref_e = station_vel.e;
                velref_n = station_vel.n;
                if zero_vertical:
                    velref_u = 0;
                else:
                    velref_u = station_vel.u;
        return velref_e, velref_n, velref_u;


def paired_gps_geocoded_insar(gps_los_velfield, xarray, yarray, LOS_array, window_pixels=5):
    """
    Extract paired LOS InSAR and GPS velocities at pixels given by lons/lats in gps_los_velfield.
    Average over a window of pixels whose width is given by an input parameter.
    Returns non-nan values.
    """
    insar_los_array, gps_los_array, lonarray, latarray = [], [], [], [];
    distance_tolerance = 0.01;  # degrees (approximately 1 km)
    for station_vel in gps_los_velfield:
        xi, deg_distance_x = closest_index(xarray, station_vel.elon);
        yi, deg_distance_y = closest_index(yarray, station_vel.nlat);
        if deg_distance_x < distance_tolerance and deg_distance_y < distance_tolerance:
            target_array = LOS_array[yi - window_pixels:yi + window_pixels, xi - window_pixels:xi + window_pixels]
            # default tolerance is about 1 km, and it shows whether we are accidentally outside of InSAR domain
            InSAR_LOS_value = np.nanmean(target_array);
            if ~np.isnan(InSAR_LOS_value):
                insar_los_array.append(InSAR_LOS_value);
                gps_los_array.append(station_vel.e);  # the LOS velocity
                lonarray.append(station_vel.elon);
                latarray.append(station_vel.nlat);
    return np.array(insar_los_array), np.array(gps_los_array), np.array(lonarray), np.array(latarray);


# A little bit of IO for working with GPS velocity fields in LOS geometries
# -------------------------------- # 
def input_gps_as_los(filename):
    print("Reading file %s " % filename);
    gps_velfield = [];
    [elon, nlat, los_vel, name] = np.loadtxt(filename, usecols=(0, 1, 5, 6), unpack=True,
                                             dtype={'names': ('elon', 'nlat', 'los_vel', 'station_name'),
                                                    'formats': (float, float, float, 'U4')});
    for i in range(len(elon)):
        gps_station_as_los = gps_io_functions.Station_Vel(name=name[i], elon=elon[i], nlat=nlat[i], e=los_vel[i], n=0,
                                                          u=0, se=0, sn=0, su=0, first_epoch=0, last_epoch=0,
                                                          refframe=0, proccenter=0, subnetwork=0, survey=0,
                                                          meas_type='los');
        gps_velfield.append(gps_station_as_los);
    return [gps_velfield];


def output_gps_as_los(gps_velfield, LOS_velfield, outfile):
    ofile = open(outfile, 'w');
    ofile.write("# lon lat gpsE gpsN gpsU LOS name\n");
    for i in range(len(LOS_velfield)):
        ofile.write("%f %f %f %f %f %f %s \n" %
                    (LOS_velfield[i].elon, LOS_velfield[i].nlat, gps_velfield[i].e, gps_velfield[i].n,
                     gps_velfield[i].u, LOS_velfield[i].e, LOS_velfield[i].name));
    ofile.close();
    print("-->Outputs printed to %s" % outfile);
    return;
