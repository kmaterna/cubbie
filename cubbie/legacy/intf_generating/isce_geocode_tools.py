"""
June 2020
A series of functions to geocode isce images and los.rdr.geo in various formats
Including the UAVSAR stacks
And the UAVSAR ground range igram format from the JPL website
"""


import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from subprocess import call
from tectonic_utils.read_write import netcdf_read_write as rwr
from tectonic_utils.geodesy import insar_vector_functions
from ...read_write_insar_utilities import isce_read_write, jpl_uav_read_write
from . import unwrapping_isce_custom


# ------------ UTILITY FUNCTIONS -------------- #

def cut_resampled_grid(outdir, filename, variable, config_params):
    # This is for metadata like lon, lat, and lookvector
    # Given an isce file and a set of bounds to cut the file,
    # Produce the isce data and gmtsar netcdf that match each pixel.
    _, _, temp = isce_read_write.read_scalar_data(os.path.join(outdir, filename))
    print("Shape of the " + variable + " file: ", np.shape(temp))
    xbounds = [float(config_params.xbounds.split(',')[0]), float(config_params.xbounds.split(',')[1])]
    ybounds = [float(config_params.ybounds.split(',')[0]), float(config_params.ybounds.split(',')[1])]
    cut_grid = unwrapping_isce_custom.cut_grid(temp, xbounds, ybounds, fractional=True, buffer_rows=3)
    print("Shape of the cut lon file: ", np.shape(cut_grid))
    nx = np.shape(cut_grid)[1]
    ny = np.shape(cut_grid)[0]
    isce_read_write.write_isce_data(cut_grid, nx, ny, "FLOAT", os.path.join(outdir, 'cut_' + variable + '.gdal'))
    rwr.produce_output_netcdf(np.array(range(0, nx)), np.array(range(0, ny)), cut_grid,
                              "degrees", os.path.join(outdir, 'cut_' + variable + '.nc'))
    return


# ------------ GEOCODING FUNCTIONS FOR UAVSAR STACKS -------------- #
# Based on stacks of 3D netcdf's from the time series processing

def gmtsar_nc_stack_2_isce_stack(ts_file, output_dir, bands=2):
    # Decompose a 3D time series object into a series of slices
    # Write the slices into isce unwrapped format.
    os.makedirs(output_dir, exist_ok=True)
    tdata, xdata, ydata, zdata = rwr.read_3D_netcdf(ts_file)
    for i in range(np.shape(zdata)[0]):
        os.makedirs(os.path.join(output_dir, "scene_"+str(i)), exist_ok=True)
        temp = zdata[i, :, :]

        # Write data out in isce format
        ny, nx = np.shape(temp)
        name = "ts_slice_" + str(i)
        filename = os.path.join(output_dir, "scene_" + str(i), name + ".unw")
        temp = np.float32(temp)
        isce_read_write.write_isce_unw(temp, temp, nx, ny, "FLOAT", filename)

        isce_read_write.plot_scalar_data(filename, band=bands, colormap='rainbow', datamin=-50, datamax=200,
                                         aspect=1 / 5, outname=os.path.join(output_dir, "scene_" + str(i),
                                                                            "isce_unw_band.png"))
    return


def geocode_UAVSAR_stack(config_params, geocoded_folder):
    # The goals here for UAVSAR:
    # Load lon/lat grids and look vector grids
    # Resample and cut the grids appropriately
    # Write pixel-wise metadata out in the output folder
    # All these grids have only single band.
    os.makedirs(geocoded_folder, exist_ok=True)
    llh_array = np.fromfile(config_params.llh_file, dtype=np.float32)  # this is a vector.
    lkv_array = np.fromfile(config_params.lkv_file, dtype=np.float32)
    lat = llh_array[np.arange(0, len(llh_array), 3)]  # ordered array opened from the provided UAVSAR files
    lon = llh_array[np.arange(1, len(llh_array), 3)]
    _hgt = llh_array[np.arange(2, len(llh_array), 3)]
    lkv_e = lkv_array[np.arange(0, len(lkv_array), 3)]
    lkv_n = lkv_array[np.arange(1, len(lkv_array), 3)]
    lkv_u = lkv_array[np.arange(2, len(lkv_array), 3)]
    example_igram = glob.glob("../Igrams/????????_????????/*.int")[0]
    phase_array = isce_read_write.read_phase_data(example_igram)
    print("Shape of the interferogram: ", np.shape(phase_array))

    # Determine the shape of the llh array
    # assuming there's a giant gap somewhere in the lat array
    # that can tell us how many elements are in the gridded array
    typical_gap = abs(lat[1] - lat[0])
    llh_pixels_range = 1
    for i in range(1, len(lat)):
        if abs(lat[i] - lat[i - 1]) > 100 * typical_gap:
            print(lat[i] - lat[i - 1])
            print("There are %d columns in the lon/lat arrays" % i)
            llh_pixels_range = i
            break
    llh_pixels_azimuth = int(len(lon) / llh_pixels_range)
    print("llh_pixels_azimuth: ", llh_pixels_azimuth)
    print("llh_pixels_range: ", llh_pixels_range)

    # We turn the llh data into 2D arrays.
    # The look vector is in meters from the aircraft to the ground.
    lat_array = np.reshape(lat, (llh_pixels_azimuth, llh_pixels_range))
    lon_array = np.reshape(lon, (llh_pixels_azimuth, llh_pixels_range))
    lkve_array = np.reshape(lkv_e, (llh_pixels_azimuth, llh_pixels_range))
    lkvn_array = np.reshape(lkv_n, (llh_pixels_azimuth, llh_pixels_range))
    lkvu_array = np.reshape(lkv_u, (llh_pixels_azimuth, llh_pixels_range))
    lkve_array, lkvn_array, lkvu_array = insar_vector_functions.normalize_vector(lkve_array, lkvn_array, lkvu_array)
    azimuth, incidence = insar_vector_functions.calc_rdr_azimuth_incidence_from_lkv_plane_down(lkve_array, lkvn_array,
                                                                                               lkvu_array)

    # # write the data into a GDAL format.
    isce_read_write.write_isce_data(lon_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "lon_total.gdal"))
    isce_read_write.write_isce_data(lat_array, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "lat_total.gdal"))
    isce_read_write.write_isce_data(azimuth, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "azimuth_total.gdal"))
    isce_read_write.write_isce_data(incidence, llh_pixels_range, llh_pixels_azimuth, "FLOAT",
                                    os.path.join(geocoded_folder, "incidence_total.gdal"))

    # Resampling in GDAL to match the interferogram sampling
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/lon_total.gdal',
          geocoded_folder + '/lon_igram_res.tif'], shell=False)
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/lat_total.gdal',
          geocoded_folder + '/lat_igram_res.tif'], shell=False)
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/incidence_total.gdal',
          geocoded_folder + '/incidence_igram_res.tif'], shell=False)
    call(['gdalwarp', '-ts', str(np.shape(phase_array)[1]), str(np.shape(phase_array)[0]),
          '-r', 'bilinear', '-to', 'SRC_METHOD=NO_GEOTRANSFORM',
          '-to', 'DST_METHOD=NO_GEOTRANSFORM', geocoded_folder + '/azimuth_total.gdal',
          geocoded_folder + '/azimuth_igram_res.tif'], shell=False)

    # Cut the data, and quality check.
    # Writing the cut lon/lat into new files.
    cut_resampled_grid(geocoded_folder, "lon_igram_res.tif", "lon", config_params)
    cut_resampled_grid(geocoded_folder, "lat_igram_res.tif", "lat", config_params)
    cut_resampled_grid(geocoded_folder, "incidence_igram_res.tif", "incidence", config_params)
    cut_resampled_grid(geocoded_folder, "azimuth_igram_res.tif", "azimuth", config_params)

    isce_read_write.plot_scalar_data(os.path.join(geocoded_folder, 'cut_lat.gdal'),
                                     colormap='rainbow', aspect=1 / 4,
                                     outname=os.path.join(geocoded_folder, 'cut_lat_geocoded.png'))
    _, _, cut_lon = isce_read_write.read_scalar_data(os.path.join(geocoded_folder, 'cut_lon.gdal'))
    _, _, cut_lat = isce_read_write.read_scalar_data(os.path.join(geocoded_folder, 'cut_lat.gdal'))
    W, E = np.min(cut_lon), np.max(cut_lon)
    S, N = np.min(cut_lat), np.max(cut_lat)

    # This last thing may not work when finding the reference pixel, only when geocoding at the very last.
    # Double-checking the shape of the interferogram data (should match!)
    _, _, signalspread = isce_read_write.read_scalar_data(os.path.join(config_params.ts_output_dir,
                                                                       'signalspread_cut.nc'))
    print("For comparison, shape of cut data is: ", np.shape(signalspread))

    return W, E, S, N


def create_isce_stack_unw_geo(geocoded_dir, w, e, s, n):
    # With pixel-wise lat and lon and lookvector information,
    # Can we make isce geocoded unwrapped .unw.geo / .unw.geo.xml
    # geocodeGdal.py -l cut_lat.gdal -L cut_lon.gdal -f cut_something.gdal -b "S N W E"
    # After that, the BIL arrangement can be switched to BSQ,
    # So I need to make an adjustment
    folders = glob.glob(os.path.join(geocoded_dir, "scene*"))
    for folder_i in folders:
        # Run the geocode command.
        # This places the geocoded .unw.geo into each subdirectory.
        datafile = glob.glob(os.path.join(folder_i, "*.unw"))
        datafile = datafile[0]
        command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " +\
                  datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(e) + "\" -x 0.00025 -y 0.00025"
        print(command)
        print("\n")
        call(command, shell=True)

        # Unfortunately, after geocodeGdal, the files end up BSQ instead of BIL.  This is necessary to reshape them.
        # For making this more streamlined, I should definitely use a regular isce_write function in the future.
        filename = datafile + ".geo"
        isce_read_write.plot_scalar_data(filename, colormap='rainbow', datamin=-50, datamax=200,
                                         outname='test_after_geocode.png', band=2)

        print("DANGER! PLEASE FIGURE OUT A SIMPLE WRITE FUNCTION FOR THIS")
        sys.exit(0)
    return


def create_isce_stack_rdr_geo(geocoded_dir, w, e, s, n):
    # Create a geocoded azimuth and geocoded incidence file
    # Then concatenate them into a two-band-file (los.rdr.geo)
    # Then update the xml metadata.
    print("Creating los.rdr.geo")
    datafile = os.path.join(geocoded_dir, "cut_azimuth.gdal")
    command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " + \
              datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(e) + "\" -x 0.00025 -y 0.00025"
    print(command + "\n")
    call(command, shell=True)
    datafile = os.path.join(geocoded_dir, "cut_incidence.gdal")
    command = "geocodeGdal.py -l " + geocoded_dir + "/cut_lat.gdal -L " + geocoded_dir + "/cut_lon.gdal " + "-f " + \
              datafile + " -b \"" + str(s) + " " + str(n) + " " + str(w) + " " + str(s) + "\" -x 0.00025 -y 0.00025"
    call(command, shell=True)
    _, _, grid_inc = isce_read_write.read_scalar_data(os.path.join(geocoded_dir, "cut_incidence.gdal.geo"),
                                                      flush_zeros=False)
    _, _, grid_az = isce_read_write.read_scalar_data(os.path.join(geocoded_dir, "cut_azimuth.gdal.geo"),
                                                     flush_zeros=False)
    ny, nx = np.shape(grid_inc)
    filename = os.path.join(geocoded_dir, "los.rdr.geo")
    isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", filename)
    return


def inspect_isce(geocoded_dir):
    # What progress was made?  Plot things.
    folders = glob.glob(os.path.join(geocoded_dir, "scene*"))
    for folder_i in folders:
        datafile = glob.glob(os.path.join(folder_i, "*.unw.geo"))
        datafile = datafile[0]
        _, _, grid = isce_read_write.read_scalar_data(datafile, flush_zeros=False)
        print("Statistics:")
        print("shape: ", np.shape(grid))
        print("max: ", np.nanmax(grid))
        print("min: ", np.nanmin(grid))
        isce_read_write.plot_scalar_data(datafile, colormap="rainbow", datamin=-50, datamax=200,
                                         outname=os.path.join(folder_i, "geocoded_data.png"))
    return


def fix_hacky_BSQ_BIL_problem(geocoded_directory, mynum):
    # August 2020
    # This script is meant to un-do something that happened before on NoMachine
    # The .unw.geo files ended up BSQ instead of BIL
    # So we need to fix it.
    # At the end, the new .unw.geo and xml should be properly formatted
    # If we fix the end of the isce_geocode script for nomachine, than this should never be necessary.

    # Find the files
    unw_file = geocoded_directory + 'ts_slice_' + mynum + '.unw.geo'
    unw_xml = unw_file + '.xml'
    unw_file_final = os.path.join(geocoded_directory + 'BIL_correct', 'ts_slice_' + mynum + '.unw.geo')

    # Read the problematic bands and get ready to package them into a real geocoded file.
    _, _, data_top = isce_read_write.read_scalar_data(unw_file, band=1)
    # I'm not even sure how I'm allowed to read band 2 (xml says 1 band).
    _, _, data_bottom = isce_read_write.read_scalar_data(unw_file, band=2)  # xml is clearly wrong.
    data = np.vstack((data_top, data_bottom))  # each of these has a duplicate row by accident.
    data_surviving = np.zeros(np.shape(data_top))
    for i in range(np.shape(data)[0]):
        counter = int(np.floor(i / 2.0))
        data_surviving[counter, :] = data[i, :]
    (ny, nx) = np.shape(data_surviving)

    firstLon, firstLat, dE, dN, xmin, xmax = isce_read_write.get_xmin_xmax_xinc_from_xml(unw_xml)

    isce_read_write.write_isce_unw(data_surviving, data_surviving, nx, ny, "FLOAT", unw_file_final, firstlat=firstLat,
                                   firstlon=firstLon, deltalon=dE, deltalat=dN, xmin=xmin, xmax=xmax)
    return


# ------------ JPL UAVSAR IGRAM FORMATS -------------- # 
# A set of tools designed for handling of ground-range igrams
# from the JPL website for UAVSAR individual igram products

def get_geocoded_axes_from_ann(ann_file, cut_rowcol, looks_x, looks_y):
    # Given .ann file and cutting/multilooking scheme, give us the ground-range points of the final pixels
    # in two east-and-north axes
    # cut_rowcol is an array specifying our cut range.
    # Example: [2500, 5100, 7800, 13000] where 0-1 are rows and 2-3 are cols.
    # looks_x and looks_y were used in filtering.
    num_rows, num_cols = jpl_uav_read_write.get_rows_cols(ann_file, 'ground')
    start_lon, start_lat, lon_inc, lat_inc = jpl_uav_read_write.get_ground_range_corner_increment(ann_file)
    x_orig = [start_lon + i * lon_inc for i in range(0, num_cols)]
    y_orig = [start_lat + i * lat_inc for i in range(0, num_rows)]
    x_cut = x_orig[cut_rowcol[2]: cut_rowcol[3]]
    y_cut = y_orig[cut_rowcol[0]: cut_rowcol[1]]  # implement the grid cut

    # next, implement the multilooking
    x_filt, y_filt = [], []
    counter = np.arange(0, len(x_cut), looks_x)
    for i in range(len(counter)):
        region = np.mean(x_cut[counter[i]:counter[i] + looks_x])
        x_filt.append(region)
    counter = np.arange(0, len(y_cut), looks_y)
    for i in range(len(counter)):
        region = np.mean(y_cut[counter[i]:counter[i] + looks_y])
        y_filt.append(region)

    return x_filt, y_filt


def write_unwrapped_ground_range_displacements(ground_range_phase_file, output_file, x_axis, y_axis, wavelength):
    # Given a file with ground range pixels in unwrapped phase,
    # Multiply by wavelength
    # Write the response into a unw.geo file with special xml
    lon_inc = x_axis[1] - x_axis[0]
    lat_inc = y_axis[1] - y_axis[0]

    [_, _, unw] = rwr.read_netcdf4(ground_range_phase_file)

    plt.figure(figsize=(11, 7), dpi=300)
    X, Y = np.meshgrid(x_axis, y_axis)
    plt.pcolormesh(X, Y, unw, cmap='jet', vmin=0, vmax=20)
    plt.colorbar()
    plt.savefig('unwrapped_geocoded_phase.png')

    # CONVERT TO MM using the wavelength of UAVSAR
    unw = np.multiply(unw, wavelength / (4 * np.pi))
    (ny, nx) = np.shape(unw)

    # ISCE UNW.GEO (IN MM)
    isce_read_write.write_isce_unw(unw, unw, nx, ny, "FLOAT", output_file, firstlat=max(y_axis), firstlon=min(x_axis),
                                   deltalon=lon_inc, deltalat=lat_inc, xmin=min(x_axis),
                                   xmax=max(x_axis))  # 2 bands, floats
    return
