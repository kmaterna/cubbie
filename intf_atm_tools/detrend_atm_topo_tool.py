#!/usr/bin/env python

"""
A script to take a co-registered image and a DEM,
Solve for best-fitting linear trend globally across the whole scene.
Remove trend and save the adjusted stack to a file.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from Tectonic_Utils.read_write import netcdf_read_write as rw


def detrend_topo_whole_box(igram_filename, demfile, outfilename):
    """ Driver for the whole command-line program. """
    [_, _, demdata] = rw.read_any_grd(demfile);
    [xdata, ydata, phase_data] = rw.read_any_grd(igram_filename);
    [phase_array_1d, dem_array_1d, corrected_array_1d, corrected_array_2d] = correct_for_topo_trend(phase_data, demdata)
    rw.produce_output_netcdf(xdata, ydata, corrected_array_2d, 'unwrapped_phase', outfilename);
    output_plots(phase_array_1d, dem_array_1d, corrected_array_1d, phase_data, corrected_array_2d, outfilename);
    return;


# ---------- COMPUTE FUNCTION ------------ #
def correct_for_topo_trend(zdata, demdata):
    """
    :param zdata: 2d array of unwrapped phase values
    :param demdata: 2d array of topography, same size as zdata
    :returns: 2d array of corrected unwrapped phase values
    """
    defensive_checks(zdata, demdata);
    full_phase_array_1d, full_dem_array_1d = np.ravel(zdata), np.ravel(demdata);
    phase_array_1d = full_phase_array_1d[np.where(~np.isnan(full_phase_array_1d))];  # non-nan data to 1D arrays
    demarray_1d = full_dem_array_1d[np.where(~np.isnan(full_phase_array_1d))];
    coef = np.polyfit(demarray_1d, phase_array_1d, 1);  # Generate a best-fitting slope between phase and topography
    corrected_array_1d = full_phase_array_1d - coef[0] * full_dem_array_1d;  # Remove slope
    corrected_array_2d = np.reshape(corrected_array_1d, np.shape(zdata));  # Reshape back into right shape
    print("Best-fitting Slope: %f " % coef[0]);
    return [full_phase_array_1d, full_dem_array_1d, corrected_array_1d, corrected_array_2d];


def defensive_checks(zdata, demdata):
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
    plt.savefig(outfilename);
    plt.close();
    return;


if __name__ == "__main__":
    """Example runstring: detrend_atm_topo_tool.py unwrap.grd dem.grd detrended_unwrap.grd"""
    if len(sys.argv) < 4:
        sys.exit(0);
    igram_filename, dem_filename, outfilename = sys.argv[1], sys.argv[2], sys.argv[3];
    detrend_topo_whole_box(igram_filename, dem_filename, outfilename);
