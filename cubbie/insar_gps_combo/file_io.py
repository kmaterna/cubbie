# A little bit of IO for working with GPS velocity fields in LOS geometries
# -------------------------------- #


import numpy as np
from tectonic_utils.read_write import netcdf_read_write
from gnss_timeseries_viewers.gps_tools import vel_functions, file_io
from gnss_timeseries_viewers.gps_tools.file_io import io_nota, io_other, io_magnet_unr


def input_gps_as_los(filename):
    """
    Read a LOS velocity field into a list of StationVel objects, placing the los in the East component.

    :param filename: string, name of file that contains LOS velocities
    :return: velocity field
    """
    print("Reading file %s " % filename)
    gps_velfield = []
    [elon, nlat, los_vel, name] = np.loadtxt(filename, usecols=(0, 1, 5, 6), unpack=True,
                                             dtype={'names': ('elon', 'nlat', 'los_vel', 'station_name'),
                                                    'formats': (float, float, float, 'U4')})
    for i in range(len(elon)):
        gps_station_as_los = vel_functions.Station_Vel(name=name[i], elon=elon[i], nlat=nlat[i], e=los_vel[i], n=0, u=0,
                                                       se=0, sn=0, su=0, first_epoch=0, last_epoch=0, refframe=0,
                                                       proccenter=0, subnetwork=0, survey=0, meas_type='los')
        gps_velfield.append(gps_station_as_los)
    print("  --> Read in %d LOS velocities." % len(gps_velfield))
    return gps_velfield


def inputs_gps_pbo_like(gps_file, coordbox_gps=(-180, 180, -90, 90)):
    """
    Read a pbo/nota type input file. Optionally filter based upon a desired bounding box.

    :param gps_file: string, filename
    :param coordbox_gps: tuple containing (W, E, S, N)
    :return: list of StationVel objects
    """
    # The velocities within the lat/lon box.
    if '.vel' in gps_file:
        [gps_velfield] = io_nota.read_gamit_velfile(gps_file)
    elif '_human_' in gps_file:
        [gps_velfield] = io_other.read_humanread_vel_file(gps_file)
    else:
        [gps_velfield] = file_io.io_nota.read_pbo_vel_file(gps_file)
    gps_velfield = vel_functions.remove_duplicates(gps_velfield)
    gps_velfield = vel_functions.clean_velfield(gps_velfield, max_horiz_sigma=2, max_vert_sigma=5,
                                                coord_box=coordbox_gps)
    return gps_velfield


def inputs_gps_magent_like(gps_file, coordbox_gps):
    # The velocities within the latlon box.
    [gps_velfield] = io_magnet_unr.read_unr_vel_file(gps_file)
    [gps_velfield] = vel_functions.remove_duplicates(gps_velfield)
    [gps_velfield] = vel_functions.clean_velfield(gps_velfield, coord_box=coordbox_gps)
    return gps_velfield


def input_insar_grdfile(geocoded_insar_file, extreme_val=1e20):
    """Read a geocoded grd file with InSAR data. Clean it of spurious values and replace them with np.nan.

    :param geocoded_insar_file: string
    :param extreme_val: float, default 1e20
    :returns: [xarray_1d, yarray_1d, data_array_2d]
    """
    [xarray, yarray, LOS_array] = netcdf_read_write.read_any_grd(geocoded_insar_file)
    LOS_array[np.where(LOS_array > extreme_val)] = np.nan  # Filter spurious values from InSAR array
    if np.nanmean(xarray) > 180:
        xarray = np.subtract(xarray, 360)  # some files come in with lon=244 instead of -115.  Fixing that.
    return [xarray, yarray, LOS_array]


def input_insar_cgm_dict(insar_dict, extreme_val=1e20):
    """
    :param insar_dict: a dictionary of CGM format
    :param extreme_val: float, default 1e20
    :return: xarray_1d, yarray_1d, data_array_2d
    """
    LOS_array = insar_dict['velocities']
    LOS_array[np.where(LOS_array > extreme_val)] = np.nan  # Filter spurious values from InSAR array
    return [insar_dict['lon'], insar_dict['lat'], LOS_array]


def inputs_insar_data(insar_struct):
    if isinstance(insar_struct, dict):
        [xarray, yarray, LOS_array] = input_insar_cgm_dict(insar_struct)  # if CGM format
    else:
        [xarray, yarray, LOS_array] = input_insar_grdfile(insar_struct[3])  # if just a regular structure
    return [xarray, yarray, LOS_array]


def input_insar_lkv(insar_data_struct):
    if isinstance(insar_data_struct, dict):
        [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_lkv_cgm_dict(insar_data_struct)  # if CGM format
    else:
        [xarray, yarray, lkv_east, lkv_north, lkv_up] = inputs_lkv_grd_files(insar_data_struct)  # if regular structure
    return [xarray, yarray, lkv_east, lkv_north, lkv_up]


def inputs_lkv_grd_files(look_vector_files):
    """
    Read three grd files with look vector information from GMTSAR.

    :param look_vector_files: a list of three strings, filenames for lkv_e.grd, lkv_n.grd, lkv_u.grd
    :return: xarray_1d, yarray_1d, lkv_e_2d, lkv_n_2d, lkv_u_2d
    """
    print("-->Reading files ", look_vector_files)
    [_, _, lkv_e] = netcdf_read_write.read_netcdf4(look_vector_files[0])
    [_, _, lkv_n] = netcdf_read_write.read_netcdf4(look_vector_files[1])
    [xarray, yarray, lkv_u] = netcdf_read_write.read_netcdf4(look_vector_files[2])
    return [xarray, yarray, lkv_e, lkv_n, lkv_u]


def inputs_lkv_cgm_dict(insar_dict):
    xarray, yarray = insar_dict["lon"], insar_dict["lat"]
    lkv_e, lkv_n, lkv_u = insar_dict["lkv_E"], insar_dict["lkv_N"], insar_dict["lkv_U"]
    return [xarray, yarray, lkv_e, lkv_n, lkv_u]


def output_gps_as_los(gps_velfield, LOS_velfield, outfile):
    """
    :param gps_velfield: a list of objects
    :param LOS_velfield: a list of objects, matching in length to the gps_velfield
    :param outfile: string, name of outfile
    """
    ofile = open(outfile, 'w')
    ofile.write("# lon lat gpsE gpsN gpsU LOS name\n")
    for i in range(len(LOS_velfield)):
        ofile.write("%f %f %f %f %f %f %s \n" %
                    (LOS_velfield[i].elon, LOS_velfield[i].nlat, gps_velfield[i].e, gps_velfield[i].n,
                     gps_velfield[i].u, LOS_velfield[i].e, LOS_velfield[i].name))
    ofile.close()
    print("-->Outputs printed to %s" % outfile)
    return
