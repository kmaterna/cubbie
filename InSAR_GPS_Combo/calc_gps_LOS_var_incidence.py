"""
Python script to convert GNSS velocities into LOS relative velocities
Using a grid of look vector components.
The incidence angle of the look vector varies across the scene!
The reference pixel must be a GPS station in the Velfield
"""
import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from gnss_timeseries_viewers.gps_tools import vel_functions, file_io, gps_objects
from . import los_projection_tools


def top_level_driver(config_dict, insar_data_struct):
    """Now using insar inputs from InSAR CGM dictionary"""
    [gps_velfield] = inputs_gps(config_dict["gps_filename"], config_dict["coordbox_gps"]);
    [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_insar(insar_data_struct);
    [LOS_velfield] = compute(gps_velfield, config_dict["reference_gps"], xarray, yarray, lkv_east, lkv_north, lkv_up);
    los_projection_tools.output_gps_as_los(gps_velfield, LOS_velfield, config_dict["outdir"] + config_dict["outfile"]);
    return;


# ------------------ INPUTS --------------------- #
def inputs_gps(gps_file, coordbox_gps):
    # The velocities within the lat/lon box.
    if '.vel' in gps_file:
        [gps_velfield] = file_io.io_nota.read_gamit_velfile(gps_file);
    elif '_human_' in gps_file:
        [gps_velfield] = file_io.io_other.read_humanread_vel_file(gps_file);
    else:
        [gps_velfield] = file_io.io_nota.read_pbo_vel_file(gps_file);
    gps_velfield = vel_functions.remove_duplicates(gps_velfield);
    gps_velfield = vel_functions.clean_velfield(gps_velfield, max_horiz_sigma=2, max_vert_sigma=5,
                                                coord_box=coordbox_gps);
    return [gps_velfield];


def inputs_insar(insar_data_struct):
    if isinstance(insar_data_struct, dict):
        [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_lkv_cgm_dict(insar_data_struct);
    else:
        [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_lkv(insar_data_struct);
    return [xarray, yarray, lkv_east, lkv_north, lkv_up];


def inputs_lkv(look_vector_files):
    print("-->Reading files ", look_vector_files);
    [_, _, lkv_e] = netcdf_read_write.read_netcdf4(look_vector_files[0])
    [_, _, lkv_n] = netcdf_read_write.read_netcdf4(look_vector_files[1]);
    [xarray, yarray, lkv_u] = netcdf_read_write.read_netcdf4(look_vector_files[2]);
    return [xarray, yarray, lkv_e, lkv_n, lkv_u];


def inputs_lkv_cgm_dict(insar_dict):
    xarray, yarray = insar_dict["lon"], insar_dict["lat"];
    lkv_e, lkv_n, lkv_u = insar_dict["lkv_E"], insar_dict["lkv_N"], insar_dict["lkv_U"];
    return [xarray, yarray, lkv_e, lkv_n, lkv_u];


# ------------------ COMPUTE --------------- #
def compute(gps_velfield, reference_gps, xarray, yarray, lkv_east, lkv_north, lkv_up):
    """
    Get look angle for each point
    Transform GPS field into LOS field (all are arguments are single values)
    """

    # Take the reference point and transform its velocity into LOS.
    ref_e, ref_n, ref_u, reflon, reflat = los_projection_tools.get_point_enu_veltuple(gps_velfield,
                                                                                      reference_pt_name=reference_gps);
    lkv_e_ref, lkv_n_ref, lkv_u_ref = get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up,
                                                                      reflon, reflat);
    [flight_angle_ref, look_angle_ref] = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e_ref,
                                                                                                    lkv_n_ref,
                                                                                                    lkv_u_ref);

    LOS_reference = los_projection_tools.simple_project_ENU_to_LOS(ref_e, ref_n, ref_u, flight_angle_ref,
                                                                   look_angle_ref);
    print("ref_e, ref_n, ref_u, ref_LOS (mm/yr): ", ref_e, ref_n, ref_u, LOS_reference);

    # Now compute relative LOS for each GPS station
    LOS_velstations = [];
    for item in gps_velfield:
        lkv_e, lkv_n, lkv_u = get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up,
                                                              item.elon, item.nlat);
        if np.isnan(lkv_e):
            one_station = gps_objects.Station_Vel(name=item.name, nlat=item.nlat, elon=item.elon, e=np.nan, n=0, u=0,
                                                  sn=item.sn, se=item.se, su=item.su, first_epoch=item.first_epoch,
                                                  last_epoch=item.last_epoch, refframe=0, proccenter=0, subnetwork=0,
                                                  survey=0, meas_type='gnss');
        else:
            [flight_angle_i, look_angle_i] = insar_vector_functions.look_vector2flight_incidence_angles(lkv_e,
                                                                                                        lkv_n,
                                                                                                        lkv_u);
            LOS_array_i = los_projection_tools.simple_project_ENU_to_LOS(item.e, item.n, item.u,
                                                                         flight_angle_i, look_angle_i);
            one_station = gps_objects.Station_Vel(name=item.name, nlat=item.nlat, elon=item.elon,
                                                  e=LOS_array_i - LOS_reference, n=0, u=0, sn=item.sn, se=item.se,
                                                  su=item.su, first_epoch=item.first_epoch,
                                                  last_epoch=item.last_epoch, refframe=0, proccenter=0,
                                                  subnetwork=0, survey=0, meas_type='gnss');
        LOS_velstations.append(one_station);

    return [LOS_velstations];


def get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up, target_lon, target_lat, tol=0.1):
    """
    Find target_lon and target_lat in a regular geocoded grid
    tol is in degrees
    """
    xi, distance_x = los_projection_tools.closest_index(xarray, target_lon);
    yi, distance_y = los_projection_tools.closest_index(yarray, target_lat);
    if abs(distance_x) > tol or abs(distance_y) > tol:
        lkv_e, lkv_n, lkv_u = np.nan, np.nan, np.nan;
    else:
        lkv_e, lkv_n, lkv_u = lkv_east[yi, xi], lkv_north[yi, xi], lkv_up[yi, xi];
    return [lkv_e, lkv_n, lkv_u];
