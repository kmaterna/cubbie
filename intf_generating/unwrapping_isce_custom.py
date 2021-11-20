"""
Deal with unwrapping in ISCE interferograms
This is a code to take interferograms, cut and mask them, interpolate over bad areas,
unwrap them, and re-mask them again.
The new unwrapped files live in an alt_unwrapped directory inside each Igram dir.
Inspired by Jiang and Lohman, 2020, S1 Salton Sea project.
It is able to re-make the interferogram stack with the right multilooks and filter.
It is meant to be called from a directory parallel to the Igram and SLC directories
in the ISCE Stack processing workflow.
February 26, 2020
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import subprocess, sys, glob
from ..read_write_insar_utilities import isce_read_write
from ..math_tools import mask_and_interpolate


def add_plot(axarr, plotnumber, data, title, colormap, aspect=1/7, is_complex=0, vmin=None, vmax=None):
    x, y = get_axarr_numbers(2, 5, plotnumber);
    print(x, y);
    if is_complex == 1:
        data = np.angle(data);
    if vmin is not None:
        axarr[x][y].imshow(data, aspect=aspect, cmap=colormap, vmin=vmin, vmax=vmax);
    else:
        axarr[x][y].imshow(data, aspect=aspect, cmap=colormap);
    axarr[x][y].set_title(title, fontsize=20);
    return axarr;


def add_rectangle(axarr, plotnumber, xbounds, ybounds):
    x, y = get_axarr_numbers(2, 5, plotnumber);
    ymin = 1-ybounds[0];
    ymax = 1-ybounds[1];  # some of these axes coordinate systems increase downward
    xmin = xbounds[0];
    xmax = xbounds[1];
    axarr[x][y].plot([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], linewidth=2, color='blue',
                     transform=axarr[x][y].transAxes);
    return axarr;


def get_axarr_numbers(rows, cols, idx):
    # Given an incrementally counting idx number and a subplot dimension, where is our plot?
    total_plots = rows*cols;
    col_num = np.mod(idx, cols);
    row_num = int(np.floor(idx/cols));
    return row_num, col_num;


def write_local_iscestack_config(source_file, target_file, date_string, alt_unwrapping=0, rlks=None, alks=None,
                                 filt=None):
    # A small function to take the default ISCE config file and copy it with modifications
    # into an interferogram directory.
    ifile = open(source_file, 'r');
    ofile = open(target_file, 'w');
    stage3 = 0;
    if alt_unwrapping == 1:  # If we're modifying the file for special unwrapping instructions
        for line in ifile:
            if "ifg : " in line:
                temp = line.split('/')[:-2];
                temp = [x+'/' for x in temp];
                newline = ''.join(temp);
                newline = newline+date_string+"/alt_unwrapped/filt_"+date_string+"_manually_masked.int\n";
                ofile.write(newline);
                stage3 = 1;
            elif "coh : " in line and stage3 == 1:
                temp = line.split('/')[:-2];
                temp = [x+'/' for x in temp];
                newline = ''.join(temp);
                newline = newline+date_string+"/alt_unwrapped/filt_"+date_string+"_cut.cor\n";
                ofile.write(newline);
            elif "unwprefix : " in line:
                temp = line.split('/')[:-2];
                temp = [x+'/' for x in temp];
                newline = ''.join(temp);
                newline = newline+date_string+"/alt_unwrapped/filt_"+date_string+"_manually_masked\n";
                ofile.write(newline);
            else:
                ofile.write(line);
    else:
        for line in ifile:
            if "rlks : " in line and rlks is not None:
                ofile.write("rlks : %d\n" % rlks);
            elif "alks : " in line and alks is not None:
                ofile.write("alks : %d\n" % alks);
            elif "strength : " in line and filt is not None:
                ofile.write("strength : %.1f\n" % filt);
            else:
                ofile.write(line);
    ifile.close();
    ofile.close();
    return;


def alt_rlks_alks_workflow(date_string, rlks, alks, filt):
    filedir = "../Igrams/"+date_string+"/";  # *** Might change based on where you're calling from
    orig_config_file = "../configs/config_igram_"+date_string;
    target_file = filedir+"config_igram_"+date_string+"_rlks";
    write_local_iscestack_config(orig_config_file, target_file, date_string, rlks=rlks, alks=alks, filt=filt);
    subprocess.call(['stripmapWrapper.py', '-c', target_file, '-s', 'Function-1', '-e', 'Function-3'], shell=False);
    # for doing everything in the interferogram stage, we'll do 1 to 3
    return;


def alt_isce_unwrapping_workflow(date_string, xbounds, ybounds, coherence_cutoff):
    """
    # Step 1: Read file, and read coherence. Plot.
    # Step 2: Perform appropriate cut. Plot.
    # Step 3: Perform appropriate mask. Plot.
    # Step 4: Perform interpolation. Plot.
    # Step 5: Unwrap. Plot.
    # Step 6: Save a 10-panel plot with all stages of processing.
    # Meant to be called from the TimeSeries directory.
    """

    # CONFIGURATION PARAMETERS
    unw_max = 4*np.pi;  # how high do we let the unwrapped colorscale go?
    plot_aspect = 1/5;  # Based on multilooking, we may need an aspect ratio to make pretty plots
    filedir = "../Igrams/"+date_string+"/";  # *** This may change based on where you're calling from?
    filestem = "filt_"+date_string;
    orig_config_file = "../configs/config_igram_"+date_string;  # *** this may change
    alt_filedir = filedir+"alt_unwrapped/";
    subprocess.call(["mkdir", "-p", alt_filedir], shell=False);

    f, axarr = plt.subplots(2, 5, figsize=(18, 16));

    # Step 1: Read the automatic data
    slc = isce_read_write.read_complex_data(filedir + filestem + ".int");
    cor = isce_read_write.read_scalar_data(filedir + filestem + ".cor");
    # for unwrapped files, band = 2
    orig_unw = isce_read_write.read_scalar_data(filedir + filestem + "_snaphu.unw", band=2);
    orig_comps = isce_read_write.read_scalar_data(filedir + filestem + "_snaphu.unw.conncomp");
    axarr = add_plot(axarr, 0, slc, 'phasefilt', colormap='rainbow', is_complex=1);
    axarr = add_rectangle(axarr, 0, xbounds, ybounds);
    axarr = add_plot(axarr, 1, cor, 'coherence', colormap='gray', is_complex=0);
    axarr = add_rectangle(axarr, 1, xbounds, ybounds);

    # Step 2: Cut data
    slc_cut = mask_and_interpolate.cut_grid(slc, xbounds, ybounds, buffer_rows=3);
    cor_cut = mask_and_interpolate.cut_grid(cor, xbounds, ybounds, buffer_rows=3);
    unw_cut = mask_and_interpolate.cut_grid(orig_unw, xbounds, ybounds, buffer_rows=3);
    comps_cut = mask_and_interpolate.cut_grid(orig_comps, xbounds, ybounds, buffer_rows=3);
    axarr = add_plot(axarr, 2, slc_cut, 'phasefilt_cut', colormap='rainbow', aspect=plot_aspect, is_complex=1,
                     vmin=-np.pi, vmax=np.pi);
    axarr = add_plot(axarr, 3, unw_cut, 'unwrapped', colormap='rainbow', aspect=plot_aspect, is_complex=0, vmin=0,
                     vmax=unw_max);
    axarr = add_plot(axarr, 4, comps_cut, 'Connected Components', colormap='rainbow', aspect=plot_aspect, is_complex=0,
                     vmin=0, vmax=8);
    if np.sum(np.isnan(cor_cut)) != 0:
        print("Error! There are %d nans in the Correlation grid!" % (np.sum(np.isnan(cor_cut))));
        sys.exit(0);

    # Step 3: Perform appropriate mask
    coherence_mask = mask_and_interpolate.make_coherence_mask(cor_cut, coherence_cutoff);
    # an experiment to get more pixels back
    coherence_mask_liberal = mask_and_interpolate.make_coherence_mask(cor_cut, coherence_cutoff - 0.0);
    masked_slc = mask_and_interpolate.apply_coherence_mask(slc_cut, coherence_mask, is_complex=1);
    axarr = add_plot(axarr, 5, masked_slc, 'masked_phasefilt', colormap='rainbow', aspect=plot_aspect, is_complex=1,
                     vmin=-np.pi, vmax=np.pi);

    # Step 4: Perform interpolation
    interp_array = mask_and_interpolate.interpolate_2d(masked_slc);

    # Step 5: Write the interpolated phase out.
    ny, nx = np.shape(slc_cut);
    isce_read_write.write_isce_data(interp_array, nx, ny, dtype='CFLOAT',
                                    filename=alt_filedir + filestem + "_manually_masked.int");
    manually_masked = isce_read_write.read_complex_data(alt_filedir + filestem + "_manually_masked.int");
    axarr = add_plot(axarr, 6, manually_masked, 'Interpolated', colormap='rainbow', aspect=plot_aspect, is_complex=1,
                     vmin=-np.pi, vmax=np.pi);

    # Step 5a: Write the correlation out, which we use for unwrapping.
    isce_read_write.write_isce_data(cor_cut, nx, ny, dtype='FLOAT', filename=alt_filedir + filestem + "_cut.cor");

    # Step 6: UNWRAP
    # PUT THE LOCAL UNWRAP SCRIPT IN ALT-FILEDIR
    # CALL UNWRAPPING (stripmapWrapper start = Function-3, end = Function-3).
    target_file = alt_filedir+"config_igram_"+date_string+"_local";
    write_local_iscestack_config(orig_config_file, target_file, date_string, alt_unwrapping=1);
    subprocess.call(['stripmapWrapper.py', '-c', target_file, '-s', 'Function-3', '-e', 'Function-3'], shell=False);
    post_unwrapping = isce_read_write.read_scalar_data(alt_filedir + filestem + "_manually_masked_snaphu.unw", band=2);
    axarr = add_plot(axarr, 7, post_unwrapping, 'Unwrapped', colormap='rainbow', aspect=plot_aspect, is_complex=0,
                     vmin=0, vmax=unw_max);
    comps = alt_filedir+filestem+"_manually_masked_snaphu.unw.conncomp.vrt"
    comps = isce_read_write.read_scalar_data(comps);
    axarr = add_plot(axarr, 8, comps, 'ConnectedComps', colormap='rainbow', aspect=plot_aspect, is_complex=0,
                     vmin=0, vmax=8);

    # Step 7: Re-apply the mask
    re_masked = mask_and_interpolate.apply_coherence_mask(post_unwrapping, coherence_mask_liberal, is_complex=0,
                                                          is_float32=True);
    axarr = add_plot(axarr, 9, re_masked, 'UnwrappedMasked', colormap='rainbow', aspect=plot_aspect, is_complex=0,
                     vmin=0, vmax=unw_max);

    # MUST WRITE THE FINAL MASKED UNWRAPPED PHASE BACK INTO A FILE.
    isce_read_write.write_isce_data(re_masked, nx, ny, dtype='FLOAT',
                                    filename=alt_filedir + filestem + "_fully_processed.uwrappedphase");

    # Color bar for wrapped phase.
    cbarax = f.add_axes([0.2, 0.35, 0.25, 0.8], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=-np.pi, vmax=np.pi);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(-np.pi, np.pi));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='horizontal');
    cb.set_label('Wrapped Phase (rad)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    # Color bar for unwrapped phase.
    cbarax = f.add_axes([0.6, 0.35, 0.25, 0.8], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=0, vmax=unw_max);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(0, unw_max));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='horizontal');
    cb.set_label('Unrapped Phase (rad)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(filedir+date_string+'_image_development.png')

    return re_masked;


def main_function(rlks, alks, filt, xbounds, ybounds, coherence_cutoff):
    # # For making all interferograms in their unwrapped form.
    xbounds = [float(xbounds.split(',')[0]), float(xbounds.split(',')[1])];
    ybounds = [float(ybounds.split(',')[0]), float(ybounds.split(',')[1])];
    igrams = glob.glob("../Igrams/????????_????????");  # **** this may change depending on where you are
    for i in igrams:
        date_string = i.split('/')[-1];
        alt_rlks_alks_workflow(date_string, rlks=rlks, alks=alks, filt=filt);  # re-makes the igrams
        alt_isce_unwrapping_workflow(date_string, xbounds, ybounds, coherence_cutoff);  # unwraps the igrams
    return;


if __name__ == "__main__":
    rlks = 20
    alks = 18
    filt = 1.5
    xbounds = [0.4, 1.0];  # These are fractional units
    ybounds = [0.15, 0.7];
    coherence_cutoff = 0.6;  # what value of coherence cutoff do we want for masks?

    # # For making all interferograms in their unwrapped form.
    # # Takes a little while (maybe an hour for 50 images)
    igrams = glob.glob("../Igrams/????????_????????");  # **** this may change depending on where you are
    for i in igrams:
        date_string = i.split('/')[-1];
        alt_rlks_alks_workflow(date_string, rlks=rlks, alks=alks, filt=filt);  # re-makes the igrams
        re_masked = alt_isce_unwrapping_workflow(date_string, xbounds, ybounds, coherence_cutoff);  # unwraps the igrams
    # sys.exit(0);

# # THE SIGNAL SPREAD
# cor_value=0.5;
# filepathslist=glob.glob("../Igrams/????????_????????/filt*.cor");  # *** This may change
# cor_data = rdr.reader_isce(filepathslist);
# a = stack_corr.stack_corr(cor_data, cor_value);
# netcdf_read_write.produce_output_netcdf(cor_data.xvalues, cor_data.yvalues, a, 'Percentage', 'signalspread.nc')
# netcdf_read_write.produce_output_plot('signalspread.nc', 'Signal Spread above cor='+str(cor_value),
#                                       'signalspread.png', 'Percentage of coherence', aspect=1/4, invert_yaxis=False)
