import sys
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
import stacking_utilities
import readmytupledata
import stack_corr
import netcdf_read_write as rwr
import isce_read_write
import unwrapping_isce_custom
import isce_geocode_tools
import stacking_configparser


# --------------- STEP 1: Make corrections and TS prep ------------ #
def make_corrections_isce(config_params):
    if config_params.startstage > 1:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 1:  # if we're ending before, we don't do this.
        return;
    print("Start Stage 1 - optional atm corrections");

    # For ISCE, we might want to re-make all the interferograms and unwrap them in custom fashion.
    if config_params.custom_unwrapping:
        custom_params = stacking_configparser.read_config_isce(config_params.config_file)
        unwrapping_isce_custom.main_function(custom_params.rlks, custom_params.alks, custom_params.filt,
                                             custom_params.xbounds, custom_params.ybounds,
                                             custom_params.cor_cutoff_mask);

    # WE ALSO MAKE THE SIGNAL SPREAD FOR FULL IMAGES
    intf_file_tuples = stacking_utilities.get_list_of_intf_all(config_params);
    corr_files = [x[3] for x in intf_file_tuples];
    stack_corr.drive_signal_spread_isce(corr_files, 0.5, config_params.ts_output_dir, "signalspread_full.nc");
    print("End Stage 1 - optional atm corrections\n");
    return;


# --------------- UTILITIES ------------ # 
def get_100p_pixels_get_ref(filenameslist, ref_idx, outdir):
    # This function helps you manually choose a reference pixel.
    # You might have to run through this function a number of times
    # To select your boxes and your eventual reference pixel.
    # I pick one in a stable area, outside of the deformation, ideally in a desert.

    print("Finding the pixels that are 100 percent coherent");
    # Finding pixels that are completely non-nan
    mydata = readmytupledata.reader_isce(filenameslist, band=1);
    # if it's straight from SNAPHU, band = 2
    # if it's manually processed first, band = 1
    total_pixels = np.shape(mydata.zvalues)[1] * np.shape(mydata.zvalues)[2];
    total_images = np.shape(mydata.zvalues)[0];
    xvalues = range(np.shape(mydata.zvalues)[2]);
    yvalues = range(np.shape(mydata.zvalues)[1]);
    count = 0;
    ypixels_good = [];
    xpixels_good = [];
    ypixels_options = [];
    xpixels_options = [];

    for i in range(np.shape(mydata.zvalues)[1]):
        for j in range(np.shape(mydata.zvalues)[2]):
            oneslice = mydata.zvalues[:, i, j];
            if np.sum(~np.isnan(oneslice)) == total_images:  # if we have perfect coherence
                count = count + 1;
                xpixels_good.append(j);
                ypixels_good.append(i);
                # Here we will adjust parameters until we find a reference pixel that we like.
                if 280 < j < 300 and 2500 < i < 2700:
                    xpixels_options.append(j);
                    ypixels_options.append(i);

    idx_lucky = 710;
    xref = xpixels_options[idx_lucky];
    yref = ypixels_options[idx_lucky];

    print("%d of %d (%f percent) are totally coherent. " % (count, total_pixels, 100 * (count / total_pixels)));
    print(np.shape(mydata.zvalues));
    print("%d pixels are good options for the reference pixel. " % (len(xpixels_options)));

    # Make a plot that shows where those pixels are
    # a=rwr.read_grd(outdir+'/signalspread_cut.nc');
    plt.figure();
    # plt.imshow(a,aspect=1/4, cmap='rainbow');
    plt.plot(xpixels_good, ypixels_good, '.', color='k');
    plt.plot(xpixels_options, ypixels_options, '.', color='g');
    plt.plot(xref, yref, '.', color='r');
    plt.savefig(outdir + '/best_pixels.png');
    plt.close();

    print("Based on 100p pixels, selecting reference pixel at row/col %d, %d " % (yref, xref));
    print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.");

    return yref, xref;


def from_lonlat_get_rowcol(config_params):
    # Given a ref loc, get the geocoded grids and find the nearest point.
    # Return its row and column.
    # Step 1: Geocode properly based on a sample interferogram grid (done)
    # Step 2: extract nearest pixel (code in Brawley repo)
    isce_geocode_tools.geocode_UAVSAR_stack(config_params, 'geocoded');

    reflon = float(config_params.ref_loc.split(',')[0]);
    reflat = float(config_params.ref_loc.split(',')[1]);
    # Next we get the nearest pixel from the rasters
    raster_lon = isce_read_write.read_scalar_data(config_params.ts_output_dir + "/cut_lon.gdal");
    raster_lat = isce_read_write.read_scalar_data(config_params.ts_output_dir + "/cut_lat.gdal");
    i_found, j_found = stacking_utilities.get_nearest_pixel_in_raster(raster_lon, raster_lat, reflon, reflat);
    print("From lon/lat, found Row and Column at %d, %d " % (i_found, j_found));
    print("STOPPING ON PURPOSE: Please write your reference pixel in your config file.");
    return i_found, j_found;


def stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, ts_output_dir, label=""):
    # WE MAKE THE SIGNAL SPREAD FOR THE CUT IMAGES
    cor_files = [i.replace("fully_processed.unwrappedphase", "cut.cor") for i in intf_files];  # get for isce
    netcdfname = ts_output_dir + '/signalspread_cut_ref' + label + '.nc'
    cor_value = np.nan;
    cor_data = readmytupledata.reader_isce(cor_files);
    a = stack_corr.stack_corr(cor_data, cor_value);
    rwr.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', netcdfname)
    rwr.produce_output_plot(netcdfname, 'Signal Spread above cor=' + str(cor_value),
                            ts_output_dir + '/signalspread_cut_ref' + label + '.png',
                            'Percentage of coherence', aspect=1 / 4, invert_yaxis=False,
                            dot_points=[[colref], [rowref]]);
    signal_spread_ref = a[rowref, colref];
    print("Signal Spread of the reference pixel = %.2f " % signal_spread_ref);
    if signal_spread_ref < 50:
        print("WARNING: Your reference pixel has very low coherence. Consider picking a different one.");
        print("STOPPING ON PURPOSE.");
        sys.exit(0);
    return;


# --------------- STEP 2: Get Reference Pixel ------------ # 

def get_ref(config_params):
    if config_params.startstage > 2:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 2:  # if we're ending at intf, we don't do this.
        return;

    print("Start Stage 2 - Collecting referenced unwrapped");
    call(['cp', 'stacking.config', config_params.ts_output_dir], shell=False);

    # Very general, takes all files and doesn't discriminate.
    intf_file_tuples = stacking_utilities.get_list_of_intf_all(config_params);
    intf_files = [x[2] for x in intf_file_tuples];

    # If we are starting manually, we find reference pixel by using 100% pixels...
    if config_params.ref_idx == "" and config_params.ref_loc == "":
        rowref, colref = get_100p_pixels_get_ref(intf_files, config_params.ref_idx, config_params.ts_output_dir);
        sys.exit(0);
    # .... or by using a lon/lat to get the reference pixel
    if config_params.ref_idx == "" and config_params.ref_loc != "":
        rowref, colref = from_lonlat_get_rowcol(config_params);
        sys.exit(0);
    # .... or we already know the reference indices:
    if config_params.ref_idx != "":
        rowref = int(config_params.ref_idx.split("/")[0]);
        colref = int(config_params.ref_idx.split("/")[1]);
        print("From the config file, the Rowref and Colref are %d, %d\n" % (rowref, colref));

    stack_corr_for_ref_unwrapped_isce(intf_files, rowref, colref, config_params.ts_output_dir)

    print("End Stage 2 - Collecting referenced unwrapped\n");

    return;


# --------------- STEP 4: Geocoding Velocities ------------ #
def geocode_vels(config_params):
    if config_params.startstage > 4:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 4:  # if we're ending before, we don't do this.
        return;
    print("Start Stage 4 - Geocoding");
    if config_params.SAT == "UAVSAR":
        geocode_directory = config_params.ts_output_dir + "/isce_geocode";
        # Deleting the contents of this folder would be a good automatic step in the future.
        isce_geocode_tools.gmtsar_nc_stack_2_isce_stack(config_params.ts_output_dir + "/TS.nc", geocode_directory,
                                                        bands=2);  # write the TS data into isce binaries
        W, E, S, N = isce_geocode_tools.geocode_UAVSAR_stack(config_params,
                                                             geocode_directory);  # do this once or more than once
        isce_geocode_tools.create_isce_stack_unw_geo(geocode_directory, W, E, S, N);
        isce_geocode_tools.create_isce_stack_rdr_geo(geocode_directory, W, E, S, N);
        isce_geocode_tools.inspect_isce(geocode_directory);
    print("End Stage 4 - Geocoding");    
    return;
