# Python script to convert GNSS velocities into LOS relative velocities
# Using a grid of look vector components. 
# The incidence angle of the look vector varies across the scene! 
# The reference pixel must be a GPS station in the Velfield
from Tectonic_Utils.read_write import netcdf_read_write
import gps_io_functions
import gps_vel_functions
import los_projection_tools
from math_tools import lkv_trig_math


def top_level_driver(config_params):
    [gps_file, look_vector_files, reference_gps, coordbox_gps, outfile] = configure(config_params);
    [gps_velfield] = inputs_gps(gps_file, coordbox_gps);
    [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_lkv(look_vector_files);

    [LOS_velfield] = compute(gps_velfield, reference_gps, xarray, yarray, lkv_east, lkv_north, lkv_up);
    los_projection_tools.output_gps_as_los(gps_velfield, LOS_velfield, outfile);
    return;


# ----------------- CONFIGURE ----------------- #
def configure(config_params):
    print("-->Starting to project gps into LOS using variable look angles.");
    print("Config parameters: ", config_params);
    gps_file = config_params["gps_file"];
    look_vector_files = config_params["look_vector_files"];
    reference_gps = config_params["reference_gps"];
    outfile = config_params["outfile"];
    coordbox_gps = config_params["coordbox_gps"];
    return [gps_file, look_vector_files, reference_gps, coordbox_gps, outfile];


# ------------------ INPUTS --------------------- #
def inputs_gps(gps_file, coordbox_gps):
    # The velocities within the lat/lon box.
    if '.vel' in gps_file:
        [gps_velfield] = gps_io_functions.read_gamit_velfile(gps_file);
    else:
        [gps_velfield] = gps_io_functions.read_pbo_vel_file(gps_file);
    gps_velfield = gps_vel_functions.remove_duplicates(gps_velfield);
    gps_velfield = gps_vel_functions.clean_velfield(gps_velfield, max_horiz_sigma=2, max_vert_sigma=5,
                                                    coord_box=coordbox_gps);
    return [gps_velfield];


def inputs_lkv(look_vector_files):
    print("-->Reading files ", look_vector_files);
    [_, _, lkv_e] = netcdf_read_write.read_netcdf4(look_vector_files[0])
    [_, _, lkv_n] = netcdf_read_write.read_netcdf4(look_vector_files[1]);
    [xarray, yarray, lkv_u] = netcdf_read_write.read_netcdf4(look_vector_files[2]);
    return [xarray, yarray, lkv_e, lkv_n, lkv_u];


# ------------------ COMPUTE --------------- #
def compute(gps_velfield, reference_gps, xarray, yarray, lkv_east, lkv_north, lkv_up):
    # Get the look angle for each point
    # Transform GPS field into LOS field (all are arguments are single values)

    # Take the reference point and transform its velocity into LOS.
    velref_e, velref_n, velref_u, reflon, reflat = los_projection_tools.get_point_enu_veltuple(gps_velfield, reference_point_name=reference_gps);
    lkv_e_ref, lkv_n_ref, lkv_u_ref = get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up,
                                                                      reflon, reflat);
    [flight_angle_ref, look_angle_ref] = lkv_trig_math.look_vector2flight_incidence_angles(lkv_e_ref, lkv_n_ref,
                                                                                           lkv_u_ref);
    LOS_reference = los_projection_tools.simple_project_ENU_to_LOS(velref_e, velref_n, velref_u, flight_angle_ref,
                                                                   look_angle_ref);

    # Now compute relative LOS for each GPS station
    LOS_velstations = [];
    for item in gps_velfield:
        lkv_e, lkv_n, lkv_u = get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up,
                                                              item.elon, item.nlat);
        [flight_angle_i, look_angle_i] = lkv_trig_math.look_vector2flight_incidence_angles(lkv_e, lkv_n, lkv_u);
        LOS_array_i = los_projection_tools.simple_project_ENU_to_LOS(item.e, item.n,
                                                                     item.u, flight_angle_i, look_angle_i);
        one_station = los_projection_tools.Velfield(name=item.name, nlat=item.nlat, elon=item.elon,
                                                    e=LOS_array_i - LOS_reference,
                                                    n=0, u=0, sn=item.sn, se=item.se, su=item.su,
                                                    first_epoch=item.first_epoch, last_epoch=item.last_epoch);
        LOS_velstations.append(one_station);

    return [LOS_velstations];


def get_lookvectors_by_nearest_grid(xarray, yarray, lkv_east, lkv_north, lkv_up, target_lon, target_lat):
    # Find target_lon and target_lat in a regular geocoded grid
    xi, _ = los_projection_tools.closest_index(xarray, target_lon);
    yi, _ = los_projection_tools.closest_index(yarray, target_lat);
    lkv_e, lkv_n, lkv_u = lkv_east[yi, xi], lkv_north[yi, xi], lkv_up[yi, xi];
    return [lkv_e, lkv_n, lkv_u];
