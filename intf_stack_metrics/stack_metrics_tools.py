# August 2020
# Here we have a number of tools for post-analysis a stack of interferograms
# Some functions are plots.
# Some functions make masks.
# Hopefully reducing the number of times I have to write these functions over and over again. 

import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import netcdf_read_write


def produce_min_max(filename, xyz=False):
    if xyz == False:
        x, y, z = netcdf_read_write.read_netcdf4_variables(filename, 'lon', 'lat', 'z');
    else:
        x, y, z = netcdf_read_write.read_grd_xyz(filename);
    print("File:", filename);
    print("Max: ", np.nanmax(z));
    print("Min: ", np.nanmin(z));
    print("Shape: ", np.shape(z));
    return;


def make_outlier_mask_for_stack(filelist, maskfile, outlier_cutoff=1e4):
    # Make a mask that is ones and nans
    # Given a co-registered stack
    # If a pixel is above the outlier cutoff in any image of the stack, make a nanmask that masks that pixel.
    x, y, z = netcdf_read_write.read_grd_xyz(filelist[1])  # just to get the shape of the outputs
    crazy_mask = np.ones(np.shape(z));
    for ifile in filelist:
        print(ifile);
        x, y, ztemp = netcdf_read_write.read_grd_xyz(ifile);
        for i in range(len(y)):
            for j in range(len(x)):
                if abs(ztemp[i][j]) > outlier_cutoff:
                    crazy_mask[i][j] = np.nan;
    # Put all the crazy pixels into a mask (across all images in the stack).
    netcdf_read_write.produce_output_netcdf(x, y, crazy_mask, "", maskfile);
    return;


def make_residual_plot(file1, file2, plotname, histname, vmin=-20, vmax=5,
                       title1='', title2='', scalelabel='LOS Velocity', units='mm/yr', flip_sign1=False,
                       flip_sign2=False):
    """
    A basic function that takes two co-registered grids and subtracts them, showing residuals in the third panel
    and histogram of residuals in separate plot.
    """
    data1 = netcdf_read_write.read_grd(file1);
    data2 = netcdf_read_write.read_grd(file2);
    if flip_sign1:
        data1 = -1 * data1;
    if flip_sign2:
        data2 = -1 * data2;
    residuals = np.subtract(data1, data2);
    residuals_vector = np.reshape(residuals, (np.shape(residuals)[0] * np.shape(residuals)[1],));

    fig, axarr = plt.subplots(1, 3, sharey='all', figsize=(20, 8), dpi=300);
    axarr[0].imshow(data1, vmin=vmin, vmax=vmax, cmap='rainbow');
    axarr[0].tick_params(labelsize=16);
    axarr[0].set_title(title1, fontsize=20);
    axarr[0].invert_yaxis()

    axarr[1].imshow(data2, vmin=vmin, vmax=vmax, cmap='rainbow');
    axarr[1].tick_params(labelsize=16);
    axarr[1].set_title(title2, fontsize=20);
    axarr[1].invert_yaxis()

    axarr[2].imshow(residuals, vmin=-10, vmax=10, cmap='rainbow');
    axarr[2].tick_params(labelsize=16);
    axarr[2].set_title('Residuals', fontsize=20);
    axarr[2].invert_yaxis()

    # Fancy color bar #1
    cbarax = fig.add_axes([0.85, 0.08, 0.1, 0.9], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax, 0.1));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label(scalelabel + ' (' + units + ')', fontsize=18);
    cb.ax.tick_params(labelsize=16);

    # Fancy color bar for residuals
    cbarax = fig.add_axes([0.68, 0.05, 0.1, 0.9], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=-10, vmax=10);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax, 0.1));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='horizontal');
    cb.set_label('Residual (' + units + ')', fontsize=16);
    cb.ax.tick_params(labelsize=14);
    plt.savefig(plotname);
    plt.close();

    plt.figure(dpi=300, figsize=(8, 6))
    plt.hist(residuals_vector, bins=50, color='orange');
    plt.yscale('log');
    plt.ylabel('Number of Pixels', fontsize=20)
    plt.xlabel('Residuals (' + units + ')', fontsize=20)
    plt.gca().tick_params(axis='both', which='major', labelsize=16);
    plt.savefig(histname);
    plt.close();
    return;


def plot_two_general_grids(file1, file2, plotname,
                           vmin1=-20, vmax1=5, flip_sign1=False, title1='', scalelabel1='Velocity (mm/yr)',
                           vmin2=None, vmax2=None, flip_sign2=False, title2='', scalelabel2='Velocity (mm/yr)',
                           readfile=True, invert_yaxis=True, cmap='rainbow'):
    """
    A little function that plots two grid files in subplots side by side
    (they don't have to have the same registration, so no need to compute residuals)
    If readfile=True: then we read files. Otherwise, those two arguments are actually data
    """
    if readfile:
        data1 = netcdf_read_write.read_grd(file1);
        data2 = netcdf_read_write.read_grd(file2);
    else:
        data1 = file1;
        data2 = file2;
    if flip_sign1:
        data1 = -1 * data1;
    if flip_sign2:
        data2 = -1 * data2;
    if vmin2 is None:
        vmin2 = vmin1;
    if vmax2 is None:
        vmax2 = vmax1;

    # First figure
    fig, axarr = plt.subplots(1, 2, sharey='none', figsize=(15, 8), dpi=300);
    axarr[0].imshow(data1, vmin=vmin1, vmax=vmax1, cmap=cmap);
    axarr[0].tick_params(labelsize=16);
    axarr[0].set_title(title1, fontsize=20);
    if invert_yaxis:
        axarr[0].invert_yaxis()

    # Second figure
    axarr[1].imshow(data2, vmin=vmin2, vmax=vmax2, cmap=cmap);
    axarr[1].tick_params(labelsize=16);
    axarr[1].set_title(title2, fontsize=20);
    if invert_yaxis:
        axarr[1].invert_yaxis()

    # Colorbar #1
    cbarax = fig.add_axes([0.15, 0.06, 0.1, 0.9], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin1, vmax=vmax1);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap=cmap);
    custom_cmap.set_array(np.arange(vmin1, vmax1, 0.1));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='horizontal');
    cb.set_label(scalelabel1, fontsize=18);
    cb.ax.tick_params(labelsize=16);

    # Colorbar #2
    cbarax = fig.add_axes([0.58, 0.06, 0.1, 0.9], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin2, vmax=vmax2);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap=cmap);
    custom_cmap.set_array(np.arange(vmin2, vmax2, 0.1));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='horizontal');
    cb.set_label(scalelabel2, fontsize=18);
    cb.ax.tick_params(labelsize=16);

    plt.savefig(plotname);
    plt.close();
    return;


def histogram_of_grd_file_values(filename, varname='Deviation', plotname='histogram_values.png'):
    """
    simple plot to make a histogram of a grid file
    """
    z = netcdf_read_write.read_grd(filename);
    z_vector = np.reshape(z, (np.shape(z)[0] * np.shape(z)[1],))
    plt.figure(dpi=250, figsize=(8, 7));
    plt.hist(z_vector, bins=50, color='orange');
    plt.yscale('log');
    plt.ylabel('Number of Pixels', fontsize=20)
    plt.xlabel(varname, fontsize=20)
    plt.gca().tick_params(axis='both', which='major', labelsize=16);
    plt.savefig(plotname);
    plt.close();
    return;


def scatterplot_of_grd_values(data1, data2, plotname='scatter.png', xlabel='', ylabel=''):
    """
    One-to-one scatter plot
    """
    xy = np.shape(data1);
    length = xy[0] * xy[1];
    data1 = data1.reshape((length,));
    data2 = data2.reshape((length,));

    plt.figure(figsize=(10, 10), dpi=300);
    plt.plot(data1, data2, '.', markersize=0.5);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    plt.xlim([0, 100])
    plt.grid(True);
    plt.savefig(plotname);
    plt.close();
    return;


def all_gridded_histograms(data_all, date_pairs):
    """
    Here, we make a series of 12-panel plots that histogram the values in each interferogram
    data_all is a list of 2d array data for each intf
    date_pairs is a list of strings like '2016292_2016316' for each intf
    """
    num_plots_x = 4;
    num_plots_y = 3;

    for i in range(len(data_all)):
        if np.mod(i, num_plots_y * num_plots_x) == 0:
            count = i;
            fignum = i / (num_plots_y * num_plots_x);
            f, axarr = plt.subplots(num_plots_y, num_plots_x, figsize=(10, 10));
            for k in range(num_plots_y):
                for m in range(num_plots_x):
                    if count == len(data_all):
                        break;

                    # How many days separate this interferogram?
                    day1 = date_pairs[count].split('_')[0];
                    day2 = date_pairs[count].split('_')[1];
                    dt1 = dt.datetime.strptime(day1, '%Y%j');
                    dt2 = dt.datetime.strptime(day2, '%Y%j');
                    deltat = dt2 - dt1;
                    daysdiff = deltat.days;

                    numelements = np.shape(data_all[count]);
                    mycorrs = np.reshape(data_all[count], (numelements[0] * numelements[1], 1));
                    nonans = mycorrs[~np.isnan(mycorrs)];
                    above_threshold = mycorrs[np.where(mycorrs > 0.2)];
                    above_threshold = int(len(above_threshold) / 1000);

                    axarr[k][m].hist(nonans);

                    axarr[k][m].set_title(str(date_pairs[count]) + '   ' + str(daysdiff) + ' days', fontsize=8);
                    axarr[k][m].set_yscale('log');
                    axarr[k][m].set_xlim([0, 1]);
                    axarr[k][m].set_ylim([1, 1000 * 1000]);
                    axarr[k][m].plot([0.1, 0.1], [0, 1000 * 1000], '--r');
                    axarr[k][m].text(0.75, 200 * 1000, str(above_threshold) + 'K', fontsize=8);

                    count = count + 1;
            plt.savefig("corr_hist_" + str(fignum) + ".eps");
            plt.close();
            print("writing corr_hist_" + str(fignum) + ".eps");
    return;
