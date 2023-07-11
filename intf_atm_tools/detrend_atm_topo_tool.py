#!/usr/bin/env python

"""
Post-processing of InSAR phase data, i.e., removing ramps and topography-correlated trend.
Pass in a co-registered image and a DEM. Solve for best-fitting linear trend globally across the whole scene.
Remove the topography-correlated trend and save the adjusted image to a file.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from Tectonic_Utils.read_write import netcdf_read_write as rw


def arg_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--data_file', type=str, help='''filename for phase information, grd file''', required=True);
    p.add_argument('-d', '--dem_file', type=str, help='''filename for dem information, grd file''', required=True);
    p.add_argument('-o', '--outname', type=str, help='''Output filename for corrected phase information, grd file''',
                   required=True);
    p.add_argument('-r', '--detrend_topography', type=bool, default=True,
                   help='''Remove topography-correlated trend, default is True''');
    p.add_argument('-m', '--mask_polygon', type=str, help='''Future: will have a masked polygon feature''');
    p.add_argument('-p', '--remove_xy_plane', type=bool, default=False,
                   help='''Future: will have planar removal feature, default is False''');
    p.add_argument('-c', '--coherence', type=str, help='''Future: will have a coherence file, grd file''');
    p.add_argument('-t', '--coherence_cutoff', type=float, help='''Future: will have a coherence mask cutoff''');
    exp_dict = vars(p.parse_args())
    return exp_dict;


def coordinator(exp_dict):
    """ Driver for the whole command-line program. """
    xdata, ydata, phase_data, demdata = data_reader(exp_dict);
    corrected_phase_2d = phase_data.copy();
    if exp_dict['detrend_topography']:
        [phase_1d, dem_1d, corrected_phase_1d, corrected_phase_2d] = correct_for_topo_trend(corrected_phase_2d, demdata)
        output_plots(phase_1d, dem_1d, corrected_phase_1d, phase_data, corrected_phase_2d, exp_dict['outname']);
    if exp_dict['remove_xy_plane']:
        corrected_phase_2d = correct_for_plane(xdata, ydata, corrected_phase_2d);
    rw.produce_output_netcdf(xdata, ydata, corrected_phase_2d, 'unwrapped_phase', exp_dict['outname']);
    return;


def data_reader(exp_dict):
    [_, _, demdata] = rw.read_any_grd(exp_dict['dem_file']);
    [xdata, ydata, phase_data] = rw.read_any_grd(exp_dict['data_file']);
    return xdata, ydata, phase_data, demdata;


# ---------- COMPUTE FUNCTIONS ------------ #
def correct_for_topo_trend(phasedata, demdata):
    """
    :param phasedata: 2d array of unwrapped phase values
    :param demdata: 2d array of topography, same size as zdata
    :returns: 2d array of corrected unwrapped phase values
    """
    print("Removing topography-correlated trend.")
    defensive_checks(phasedata, demdata);
    full_phase_array_1d, full_dem_array_1d = np.ravel(phasedata), np.ravel(demdata);
    phase_array_1d = full_phase_array_1d[np.where(~np.isnan(full_phase_array_1d))];  # non-nan data to 1D arrays
    demarray_1d = full_dem_array_1d[np.where(~np.isnan(full_phase_array_1d))];
    coef = np.polyfit(demarray_1d, phase_array_1d, 1);  # Generate a best-fitting slope between phase and topography
    corrected_array_1d = full_phase_array_1d - coef[0] * full_dem_array_1d;  # Remove slope
    corrected_array_2d = np.reshape(corrected_array_1d, np.shape(phasedata));  # Reshape back into right shape
    print("Best-fitting Slope: %f " % coef[0]);
    return [full_phase_array_1d, full_dem_array_1d, corrected_array_1d, corrected_array_2d];


def correct_for_plane(_xdata, _ydata, phasedata):
    """
    :param _xdata: 1d array
    :param _ydata: 1d array
    :param phasedata: 2d array
    :returns: 2d array of corrected unwrapped phase values
    """
    print("Removing bilinear plane.")
    corrected_phasedata = phasedata.copy();
    return corrected_phasedata;


def defensive_checks(zdata, demdata):
    """
    :param zdata: 2d array of unwrapped phase
    :param demdata: 2d array of dem information
    """
    if np.shape(zdata) != np.shape(demdata):
        raise ValueError("Error! Phase and Topography arrays do not have the same shape.");
    if np.sum(np.isnan(zdata)) == np.shape(zdata)[0] * np.shape(zdata)[1]:
        raise ValueError("Error! Phase contains only nans");
    print('Defensive checks passed');
    return;


# --------------- OUTPUTS ---------------- #
def output_plots(phase_array_1d, dem_array_1d, corrected_array_1d, phase_2d, corrected_array_2d, outgrdname):
    plt.figure();
    plt.plot(dem_array_1d, phase_array_1d, '.', label='Before correction');
    plt.plot(dem_array_1d, corrected_array_1d, '.r', alpha=0.15, label='After correction')
    plt.ylabel('phase');
    plt.xlabel('topo');
    plt.title('Initial and Corrected Phase vs. Topography');
    plt.legend();
    outfilename = outgrdname.split('.grd')[0] + '_phase_topo.png'
    print("Saving figure %s " % outfilename);
    plt.savefig(outfilename);
    plt.close();

    plt.figure();
    plt.subplot(1, 2, 1);
    plt.imshow(phase_2d, cmap='hsv', vmin=np.min(phase_array_1d), vmax=np.max(phase_array_1d));
    plt.title('Before Correction');
    plt.subplot(1, 2, 2);
    plt.imshow(corrected_array_2d, cmap='hsv', vmin=np.min(phase_array_1d), vmax=np.max(phase_array_1d));
    plt.title('After Correction');
    cb = plt.colorbar();
    cb.set_label("Unwrapped Phase (radians)", size=12);
    outfilename = outgrdname.split('.grd')[0] + '_before_after.png'
    print("Saving figure %s " % outfilename);
    plt.savefig(outfilename);
    plt.close();
    return;


if __name__ == "__main__":
    """Example runstring: detrend_atm_topo_tool.py detrend_atm_topo_tool.py unwrap.grd dem.grd detrended_unwrap.grd"""
    exp_dict = arg_parser();
    coordinator(exp_dict);
