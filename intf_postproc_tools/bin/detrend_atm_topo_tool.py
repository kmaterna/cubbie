#!/usr/bin/env python

"""
Post-processing of InSAR phase data, i.e., removing ramps and topography-correlated trend.
Pass in a GRD file of unwrapped phase, and possibly a co-registered DEM.
Option: Solve for best-fitting linear trend globally across the whole scene.
Option: Remove the topography-correlated trend and save the adjusted image to a file.
"""

import numpy as np
import argparse, scipy
from Tectonic_Utils.read_write import netcdf_read_write as rw
from S1_batches.intf_postproc_tools import plots

def arg_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--data_file', type=str,
                   help='''filename for phase information, grd file, REQUIRED''', required=True);
    p.add_argument('-o', '--outname', type=str,
                   help='''Output filename for corrected phase information, grd file, REQUIRED''', required=True);
    p.add_argument('-w', '--produce_plots', action="store_true",
                   help='''Plot the before-and-after images, default is False''');
    p.add_argument('-r', '--detrend_topography', action="store_true",
                   help='''Remove topography-correlated trend, default is False''');
    p.add_argument('-d', '--dem_file', type=str, help='''filename for dem information, grd file''', required=False);
    p.add_argument('-p', '--remove_xy_plane', action="store_true", help='''Planar removal feature, default is False''');
    p.add_argument('-m', '--mask_polygon', type=str, help='''Future: will have a masked polygon feature''');
    p.add_argument('-c', '--coherence', type=str, help='''Future: will have a coherence file, grd file''');
    p.add_argument('-t', '--coherence_cutoff', type=float, help='''Future: will have a coherence mask cutoff''');
    exp_dict = vars(p.parse_args())
    return exp_dict;


def coordinator(exp_dict):
    """ Driver for the whole command-line program. """
    [xdata, ydata, phase_data] = rw.read_any_grd(exp_dict['data_file']);
    corrected_phase_2d = phase_data.copy();
    if exp_dict['detrend_topography']:
        [_, _, demdata] = rw.read_any_grd(exp_dict['dem_file']);
        corrected_phase_2d = correct_for_topo_trend(corrected_phase_2d, demdata, exp_dict);
    if exp_dict['remove_xy_plane']:
        corrected_phase_2d = correct_for_plane(xdata, ydata, corrected_phase_2d, exp_dict);
    rw.produce_output_netcdf(xdata, ydata, corrected_phase_2d, 'unwrapped_phase', exp_dict['outname']);
    return;


# ---------- COMPUTE FUNCTIONS ------------ #
def correct_for_topo_trend(phasedata, demdata, exp_dict):
    """
    :param phasedata: 2d array of unwrapped phase values
    :param demdata: 2d array of topography, same size as zdata
    :param exp_dict: dictionary of parameters, including outfilename
    :returns: 2d array of corrected unwrapped phase values
    """
    print("Removing topography-correlated trend.")
    defensive_checks_topo_phase(phasedata, demdata);
    full_phase_array_1d, full_dem_array_1d = np.ravel(phasedata), np.ravel(demdata);
    phase_array_1d = full_phase_array_1d[np.where(~np.isnan(full_phase_array_1d))];  # non-nan data to 1D arrays
    demarray_1d = full_dem_array_1d[np.where(~np.isnan(full_phase_array_1d))];
    coef = np.polyfit(demarray_1d, phase_array_1d, 1);  # Generate a best-fitting slope between phase and topography
    corrected_array_1d = full_phase_array_1d - coef[0] * full_dem_array_1d;  # Remove slope
    corrected_array_2d = np.reshape(corrected_array_1d, np.shape(phasedata));  # Reshape back into right shape
    print("Best-fitting Slope: %f " % coef[0]);
    if exp_dict["produce_plots"]:
        plots.before_after_images(phasedata, corrected_array_2d,
                                  outfilename=exp_dict['outname'].split('.grd')[0] + '_before_after_topo.png');
        plots.linear_topo_phase_plot(phase_array_1d, demarray_1d, corrected_array_1d,
                                     outfilename=exp_dict['outname'].split('.grd')[0] + '_phase_topo.png');
    return corrected_array_2d;


def correct_for_plane(xdata, ydata, phasedata, exp_dict):
    """
    :param xdata: 1d array
    :param ydata: 1d array
    :param phasedata: 2d array
    :param exp_dict: dictionary of parameter values, including outfilename
    :returns: 2d array of corrected unwrapped phase values
    """
    print("Removing bilinear plane.")
    X, Y = np.meshgrid(xdata, ydata);
    full_x, full_y, full_phase = np.ravel(X), np.ravel(Y), np.ravel(phasedata);
    x_nonans = full_x[np.where(~np.isnan(full_phase))];
    y_nonans = full_y[np.where(~np.isnan(full_phase))];
    phase_nonans = full_phase[np.where(~np.isnan(full_phase))];
    A = np.vstack((x_nonans, y_nonans, np.ones(np.shape(x_nonans)))).T;
    coeffs = scipy.linalg.lstsq(np.array(A), np.array(phase_nonans));  # the actual optimization step
    model_coeffs = coeffs[0];  # model: [z = ax + by + c]
    full_phase_model = model_coeffs[0]*full_x + model_coeffs[1] * full_y + model_coeffs[2] * np.ones(np.shape(full_y));
    full_corrected_phase = np.subtract(full_phase, full_phase_model);
    corrected_phase_2d = np.reshape(full_corrected_phase, np.shape(phasedata));
    if exp_dict["produce_plots"]:
        plots.before_after_images(phasedata, corrected_phase_2d,
                                  outfilename=exp_dict['outname'].split('.grd')[0] + '_before_after_planar.png');
    return corrected_phase_2d;


def defensive_checks_topo_phase(zdata, demdata):
    """
    :param zdata: 2d array of unwrapped phase
    :param demdata: 2d array of dem information
    """
    if np.shape(zdata) != np.shape(demdata):
        raise ValueError("Error! Phase and Topography arrays do not have the same shape.");
    if np.sum(np.isnan(zdata)) == np.shape(zdata)[0] * np.shape(zdata)[1]:
        raise ValueError("Error! Phase contains only nans");
    print('Defensive checks passed for topography-correlated trend');
    return;


if __name__ == "__main__":
    exp_dict = arg_parser();
    coordinator(exp_dict);
