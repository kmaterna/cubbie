"""
June 2020
A series of functions to geocode isce images and los.rdr.geo in various formats
Including the UAVSAR stacks
And the UAVSAR ground range igram format from the JPL website
"""


import numpy as np
import matplotlib.pyplot as plt
import sys, glob, os
from subprocess import call
from Tectonic_Utils.read_write import netcdf_read_write as rwr
from Tectonic_Utils.geodesy import haversine, insar_vector_functions
from ..read_write_insar_utilities import isce_read_write, jpl_uav_read_write
from . import unwrapping_isce_custom


# ------------ UTILITY FUNCTIONS -------------- #

def cut_resampled_grid(outdir, filename, variable, config_params):
    # This is for metadata like lon, lat, and lookvector
    # Given an isce file and a set of bounds to cut the file,
    # Produce the isce data and gmtsar netcdf that match each pixel.
    _, _, temp = isce_read_write.read_scalar_data(os.path.join(outdir, filename));
    print("Shape of the " + variable + " file: ", np.shape(temp));
    xbounds = [float(config_params.xbounds.split(',')[0]), float(config_params.xbounds.split(',')[1])];
    ybounds = [float(config_params.ybounds.split(',')[0]), float(config_params.ybounds.split(',')[1])];
    cut_grid = unwrapping_isce_custom.cut_grid(temp, xbounds, ybounds, fractional=True, buffer_rows=3);
    print("Shape of the cut lon file: ", np.shape(cut_grid));
    nx = np.shape(cut_grid)[1];
    ny = np.shape(cut_grid)[0];
    isce_read_write.write_isce_data(cut_grid, nx, ny, "FLOAT", os.path.join(outdir, 'cut_' + variable + '.gdal'));
    rwr.produce_output_netcdf(np.array(range(0, nx)), np.array(range(0, ny)), cut_grid,
                              "degrees", os.path.join(outdir, 'cut_' + variable + '.nc'));
    return;


# ------------ GEOCODING FUNCTIONS FOR UAVSAR STACKS -------------- #
# Based on stacks of 3D netcdf's from the time series processing

def gmtsar_nc_stack_2_isce_stack(ts_file, output_dir, bands=2):
    # Decompose a 3D time series object into a series of slices
    # Write the slices into isce unwrapped format.
    os.makedirs(output_dir, exist_ok=True);
    tdata, xdata, ydata, zdata = rwr.read_3D_netcdf(ts_file);
    for i in range(np.shape(zdata)[0]):
        os.makedirs(os.path.join(output_dir, "scene_"+str(i)), exist_ok=True);
        temp = zdata[i, :, :];

        # Write data out in isce format
        ny, nx = np.shape(temp);
        name = "ts_slice_" + str(i);
        filename = os.path.join(output_dir, "scene_" + str(i), name + ".unw");
        temp = np.float32(temp);
        isce_read_write.write_isce_unw(temp, temp, nx, ny, "FLOAT", filename);

        isce_read_write.plot_scalar_data(filename, band=bands, colormap='rainbow', datamin=-50, datamax=200,
                                         aspect=1 / 5, outname=os.path.join(output_dir, "scene_" + str(i),
                                                                            "isce_unw_band.png"));
    return;


def geocode_UAVSAR_stack(config_params, geocoded_folder):
    # The goals here for UAVSAR:
    # Load lon/lat grids and look vector grids
    # Resample and cut the grids appropriately
    # Write pixel-wise metadata out in the output folder
    # All these grids have only single band.
    os.makedirs(geocoded_folder, exist_ok=True);
    llh_array = np.fromfile(config_params.llh_file, dtype=np.float32);  # this is a vector.
    lkv_array = np.fromfile(config_params.lkv_file, dtype=np.float32);
    lat = llh_array[np.arange(0, len(llh_array), 3)];  # ordered array opened from the provided UAVSAR files
    lon = llh_array[np.arange(1, len(llh_array), 3)];
    hgt = llh_array[np.arange(2, len(llh_array), 3)];
    lkv_e = lkv_array[np.arange(0, len(lkv_array), 3)]
    lkv_n = lkv_array[np.arange(1, len(lkv_array), 3)]
    lkv_u = lkv_array[np.arange(2, len(lkv_array), 3)]
    example_igram = glob.glob("../Igrams/????????_????????/*.int")[0];
    phase_array = isce_read_write.read_phase_data(example_igram);
    print("Shape of the interferogram: ", np.shape(phase_array));

    # Determine the shape of the llh array
    # assuming there's a giant gap somewhere in the lat array
    # that can tell us how many elements are in the gridded array
    typical_gap = abs(lat[1] - lat[0]);
    for i in range(1, len(lat)):
        if abs(lat[i] - lat[i - 1]) > 100 * typical_gap:
            print(lat[i] - lat[i - 1]);
            print("There are %d columns in the lon/lat arrays" % i);
            llh_pixels_range = i;
            break;
    llh_pixels_azimuth = int(len(lon) / llh_pixels_range);
    print("llh_pixels_azimuth: ", llh_pixels_azimuth);
    print("llh_pixels_range: ", llh_pixels_range);

    # We turn the llh data into 2D arrays.
    # The look vector is in meters from the aircraft to the ground.
    lat_array = np.reshape(lat, (llh_pixels_azimuth, llh_pixels_range));
    lon_array = np.reshape(lon, (llh_pixels_azimuth, llh_pixels_range));
    lkve_array = np.reshape(lkv_e, (llh_pixels_azimuth, llh_pixels_range));
    lkvn_array = np.reshape(lkv_n, (llh_pixels_azimuth, llh_pixels_range));
    lkvu_array = np.reshape(lkv_u, (llh_pixels_azimuth, llh_pixels_range));
    lkve_array, lkvn_array, lkvu_array = insar_vector_functions.normalize_vector(lkve_array, lkvn_array, lkvu_array);
    azimuth, incidence = insar_vector_functions.calc_rdr_azimuth_incidence_from_lkv_plane_down(lkve_array, lkvn_array,
                                                                                               lkvu_array);

    # # write the data into a GDAL format.
    isce_read_write.write_isce_data(lon_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "lon_total.gdal"));
    isce_read_write.write_isce_data(lat_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "lat_total.gdal"));
    isce_read_write.write_isce_data(azimuth, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "azimuth_total.gdal"));
    isce_read_write.write_isce_data(incidence, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "incidence_total.gdal"));

    # Resampling in GDAL to match the interferogram sampling
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/lon_total.gdal',
          geocoded_folder + '/lon_igram_res.tif'], shell=False);
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/lat_total.gdal',
          geocoded_folder + '/lat_igram_res.tif'], shell=False);
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/incidence_total.gdal',
          geocoded_folder + '/incidence_igram_res.tif'], shell=False);
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/azimuth_total.gdal',
          geocoded_folder + '/azimuth_igram_res.tif'], shell=False);

    # Cut the data, and quality check.
    # Writing the cut lon/lat into new files.
    cut_resampled_grid(geocoded_folder, "lon_igram_res.tif", "lon", config_params);
    cut_resampled_grid(geocoded_folder, "lat_igram_res.tif", "lat", config_params);
    cut_resampled_grid(geocoded_folder, "incidence_igram_res.tif", "incidence", config_params);
    cut_resampled_grid(geocoded_folder, "azimuth_igram_res.tif", "azimuth", config_params);

    isce_read_write.plot_scalar_data(os.path.join(geocoded_folder, 'cut_lat.gdal'),
                                     colormap='rainbow', aspect=1 / 4,
                                     outname=os.path.join(geocoded_folder, 'cut_lat_geocoded.png'));
    _, _, cut_lon = isce_read_write.read_scalar_data(os.path.join(geocoded_folder, 'cut_lon.gdal'));
    _, _, cut_lat = isce_read_write.read_scalar_data(os.path.join(geocoded_folder, 'cut_lat.gdal'));
    W, E = np.min(cut_lon), np.max(cut_lon);
    S, N = np.min(cut_lat), np.max(cut_lat);

    # This last thing may not work when finding the reference pixel, only when geocoding at the very last.
    # Double checking the shape of the interferogram data (should match!)
    _, _, signalspread = isce_read_write.read_scalar_data(os.path.join(config_params.ts_output_dir,
                                                                       'signalspread_cut.nc'));
    print("For comparison, shape of cut data is: ", np.shape(signalspread));

    return W, E, S, N;


def create_isce_stack_unw_geo(geocoded_dir, w, e, s, n):
    # With pixel-wise lat and lon and lookvector information,
    # Can we make isce geocoded unwrapped .unw.geo / .unw.geo.xml
    # geocodeGdal.py -l cut_lat.gdal -L cut_lon.gdal -f cut_something.gdal -b "S N W E"
    # After that, the BIL arrangement can be switched to BSQ,
    # So I need to make an adjustment
    folders = glob.glob(os.path.join(geocoded_dir, "scene*"));
    i = 0;
    for folder_i in folders:
        # Run the geocode command.
        # This places the geocoded .unw.geo into each sub-directory.
        datafile = glob.glob(os.path.join(folder_i, "*.unw"));
        datafile = datafile[0]
        command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " +\
                  datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(e) + "\" -x 0.00025 -y 0.00025"
        print(command);
        print("\n");
        call(command, shell=True);

        # Unfortunately, after geocodeGdal, the files end up BSQ instead of BIL.  This is necessary to reshape them.
        # For making this more streamlined, I should definitely use a regular isce_write function in the future.
        filename = datafile + ".geo"
        isce_read_write.plot_scalar_data(filename, colormap='rainbow', datamin=-50, datamax=200,
                                         outname='test_after_geocode.png', band=2);

        print("DANGER! PLEASE FIGURE OUT A SIMPLE WRITE FUNCTION FOR THIS");
        i = i + 1;
        sys.exit(0);
    return;


def create_isce_stack_rdr_geo(geocoded_dir, w, e, s, n):
    # Create a geocoded azimuth and geocoded incidence file
    # Then concatenate them into a two-band-file (los.rdr.geo)
    # Then update the xml metadata.
    print("Creating los.rdr.geo")
    datafile = os.path.join(geocoded_dir, "cut_azimuth.gdal")
    command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " + \
              datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(e) + "\" -x 0.00025 -y 0.00025"
    print(command + "\n");
    call(command, shell=True);
    datafile = os.path.join(geocoded_dir, "cut_incidence.gdal")
    command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " + \
              datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(s) + "\" -x 0.00025 -y 0.00025"
    call(command, shell=True);
    _, _, grid_inc = isce_read_write.read_scalar_data(os.path.join(geocoded_dir, "cut_incidence.gdal.geo"),
                                                      flush_zeros=False);
    _, _, grid_az = isce_read_write.read_scalar_data(os.path.join(geocoded_dir, "cut_azimuth.gdal.geo"),
                                                     flush_zeros=False);
    ny, nx = np.shape(grid_inc);
    filename = os.path.join(geocoded_dir, "los.rdr.geo")
    isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", filename);
    return;


def inspect_isce(geocoded_dir):
    # What progress was made?  Plot things.
    folders = glob.glob(os.path.join(geocoded_dir, "scene*"));
    for folder_i in folders:
        datafile = glob.glob(os.path.join(folder_i, "*.unw.geo"));
        datafile = datafile[0];
        _, _, grid = isce_read_write.read_scalar_data(datafile, flush_zeros=False);
        print("Statistics:")
        print("shape: ", np.shape(grid))
        print("max: ", np.nanmax(grid))
        print("min: ", np.nanmin(grid))
        isce_read_write.plot_scalar_data(datafile, colormap="rainbow", datamin=-50, datamax=200,
                                         outname=os.path.join(folder_i, "geocoded_data.png"));
    return;


def fix_hacky_BSQ_BIL_problem(geocoded_directory, mynum):
    # August 2020
    # This script is meant to un-do something that happened before on NoMachine
    # The .unw.geo files ended up BSQ instead of BIL
    # So we need to fix it.
    # At the end, the new .unw.geo and xml should be properly formatted
    # If we fix the end of the isce_geocode script for nomachine, than this should never be necessary.

    # Find the files
    unw_file = geocoded_directory + 'ts_slice_' + mynum + '.unw.geo';
    unw_xml = unw_file + '.xml';
    unw_file_final = os.path.join(geocoded_directory + 'BIL_correct', 'ts_slice_' + mynum + '.unw.geo');

    # Read the problematic bands and get ready to package them into a real geocoded file.
    _, _, data_top = isce_read_write.read_scalar_data(unw_file, band=1);
    # I'm not even sure how I'm allowed to read band 2 (xml says 1 band).
    _, _, data_bottom = isce_read_write.read_scalar_data(unw_file, band=2);  # xml is clearly wrong.
    data = np.vstack((data_top, data_bottom));  # each of these has a duplicate row by accident.
    data_surviving = np.zeros(np.shape(data_top));
    for i in range(np.shape(data)[0]):
        counter = int(np.floor(i / 2.0));
        data_surviving[counter, :] = data[i, :];
    (ny, nx) = np.shape(data_surviving);

    firstLon, firstLat, dE, dN, xmin, xmax = isce_read_write.get_xmin_xmax_xinc_from_xml(unw_xml);

    isce_read_write.write_isce_unw(data_surviving, data_surviving, nx, ny, "FLOAT", unw_file_final, firstLat=firstLat,
                                   firstLon=firstLon, deltaLon=dE, deltaLat=dN, Xmin=xmin, Xmax=xmax);
    return;


# ------------ JPL UAVSAR IGRAM FORMATS -------------- # 
# A set of tools designed for handling of ground-range igrams
# from the JPL website for UAVSAR individual igram products

def cross_track_pos(target_lon, target_lat, nearrange_lon, nearrange_lat, heading_cartesian):
    # Given the heading of a plane and the coordinates of one near-range point
    # Get the cross-track position of point in a coordinate system centered
    # at (nearrange_lon, nearrange_lat) with given heading
    distance = haversine.distance((target_lat, target_lon), (nearrange_lat, nearrange_lon));
    compass_bearing = haversine.calculate_initial_compass_bearing((nearrange_lat, nearrange_lon),
                                                                  (target_lat, target_lon));  # this comes CW from north
    theta = insar_vector_functions.bearing_to_cartesian(compass_bearing);  # angle of position vector, cartesian coords
    # heading_cartesian is the angle between the east unit vector and the flight direction
    x0 = distance * np.cos(np.deg2rad(theta));
    y0 = distance * np.sin(np.deg2rad(theta));  # in the east-north coordinate systeem
    x_prime, y_prime = insar_vector_functions.rotate_vector_by_angle(x0, y0, heading_cartesian);
    return y_prime;


def incidence_angle_trig(xtp, cross_track_max, near_inc_angle, far_inc_angle):
    # Using the incidence angles (to the vertical) at the upper and lower corners of the track,
    # what's the incidence angle at some location in between (xtp=cross-track-position)?
    # near_angle is the incidence angle between the viewing geometry and the vertical at the near-range.
    # nearcomp is the complement of that angle.
    # This function is kind of like linear interpolation, but a little bit curved
    # It solves an equation I derived on paper from the two near-range and far-range triangles in July 2020
    nearcomp = np.deg2rad(insar_vector_functions.complement_angle(near_inc_angle));
    farcomp = np.deg2rad(insar_vector_functions.complement_angle(far_inc_angle));  # angles from ground to satellite
    h = (np.tan(nearcomp) * np.tan(farcomp) * cross_track_max) / (np.tan(nearcomp) - np.tan(farcomp));
    angle_to_horizontal = np.rad2deg(np.arctan(h / (xtp + (h / np.tan(nearcomp)))));
    return insar_vector_functions.complement_angle(angle_to_horizontal);


def get_geocoded_axes_from_ann(ann_file, cut_rowcol, looks_x, looks_y):
    # Given .ann file and cutting/multilooking scheme, give us the ground-range points of the final pixels
    # in two east-and-north axes
    # cut_rowcol is an array specifying our cut range.
    # Example: [2500, 5100, 7800, 13000] where 0-1 are rows and 2-3 are cols.
    # looks_x and looks_y were used in filtering.
    num_rows, num_cols = jpl_uav_read_write.get_rows_cols(ann_file, 'ground');
    start_lon, start_lat, lon_inc, lat_inc = jpl_uav_read_write.get_ground_range_corner_increment(ann_file);
    x_orig = [start_lon + i * lon_inc for i in range(0, num_cols)];
    y_orig = [start_lat + i * lat_inc for i in range(0, num_rows)];
    x_cut = x_orig[cut_rowcol[2]: cut_rowcol[3]];
    y_cut = y_orig[cut_rowcol[0]: cut_rowcol[1]];  # implement the grid cut

    # next, implement the multilooking
    x_filt = [];
    y_filt = [];
    counter = np.arange(0, len(x_cut), looks_x)
    for i in range(len(counter)):
        region = np.mean(x_cut[counter[i]:counter[i] + looks_x])
        x_filt.append(region);
    counter = np.arange(0, len(y_cut), looks_y);
    for i in range(len(counter)):
        region = np.mean(y_cut[counter[i]:counter[i] + looks_y])
        y_filt.append(region);

    return x_filt, y_filt;


def write_unwrapped_ground_range_displacements(ground_range_phase_file, output_file, x_axis, y_axis, wavelength):
    # Given a file with ground range pixels in unwrapped phase,
    # Multiply by wavelength
    # Write the response into a unw.geo file with special xml
    lon_inc = x_axis[1] - x_axis[0];
    lat_inc = y_axis[1] - y_axis[0];

    [_, _, unw] = rwr.read_netcdf4(ground_range_phase_file);

    plt.figure(figsize=(11, 7), dpi=300)
    X, Y = np.meshgrid(x_axis, y_axis);
    plt.pcolormesh(X, Y, unw, cmap='jet', vmin=0, vmax=20);
    plt.colorbar();
    plt.savefig('unwrapped_geocoded_phase.png');

    # CONVERT TO MM using the wavelength of UAVSAR
    unw = np.multiply(unw, wavelength / (4 * np.pi));
    (ny, nx) = np.shape(unw);

    # ISCE UNW.GEO (IN MM)
    isce_read_write.write_isce_unw(unw, unw, nx, ny, "FLOAT", output_file,
                                   firstLat=max(y_axis), firstLon=min(x_axis), deltaLon=lon_inc, deltaLat=lat_inc,
                                   Xmin=min(x_axis), Xmax=max(x_axis));  # 2 bands, floats
    return;


def create_los_rdr_geo_from_ground_ann_file(ann_file, x_axis, y_axis):
    # Make los.rdr.geo given .ann file from JPL website's UAVSAR interferograms and the ground-range sample points.
    # x-axis and y-axis are the x and y arrays where los vectors will be extracted on a corresponding grid.
    near_angle, far_angle, heading = jpl_uav_read_write.get_nearrange_farrange_heading_angles(ann_file);
    heading_cartesian = insar_vector_functions.bearing_to_cartesian(heading);  # CCW from east
    print("Heading is %f degrees CW from north" % heading);
    print("Cartesian Heading is %f" % heading_cartesian)
    # Get the upper and lower left corners, so we can compute the length of the across-track extent in km
    ul_lon, ul_lat, ll_lon, ll_lat = jpl_uav_read_write.get_ground_range_left_corners(ann_file);

    cross_track_max = haversine.distance((ll_lat, ll_lon), (ul_lat, ul_lon));  # in km

    # Get the azimuth angle for the pixels looking up to the airplane
    # My own documentation says CCW from north, even though that's really strange.
    azimuth = heading_cartesian - 90;  # 90 degrees to the right of the airplane heading
    # (for the look vector from ground to plane)
    azimuth = insar_vector_functions.cartesian_to_ccw_from_north(azimuth);  # degrees CCW from North
    print("azimuth from ground to plane is:", azimuth)

    [X, Y] = np.meshgrid(x_axis, y_axis);
    (ny, nx) = np.shape(X);
    grid_az = azimuth * np.ones(np.shape(X));
    grid_inc = np.zeros(np.shape(X));

    print("Computing incidence angles for all pixels")
    for i in range(ny):
        for j in range(nx):
            xtp = cross_track_pos(X[i, j], Y[i, j], ll_lon, ll_lat,
                                  heading_cartesian);  # THIS WILL HAVE TO CHANGE FOR ASCENDING AND DESCENDING
            inc = incidence_angle_trig(xtp, cross_track_max, near_angle, far_angle);
            grid_inc[i, j] = inc;

    # Finally, write the 2 bands for los.rdr.geo
    isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", 'los.rdr.geo');

    return;
