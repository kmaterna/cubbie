"""
August 2020
Calculate the misfit between a geocoded InSAR velocity field and a GPS velocity field
that has been projected into the Line of Sight
Also make a 1-to-1 plot of LOS velocities
"""

import numpy as np
import matplotlib.pyplot as plt
from . import los_projection_tools, file_io


def top_level_driver(gps_los_file, geocoded_insar_struct, plotname, txtname, logname):
    """Geocoded insar: a structure"""
    gps_los_velfield = file_io.input_gps_as_los(gps_los_file)  # input GPS data
    [xarray, yarray, LOS_array] = file_io.inputs_insar_data(geocoded_insar_struct)  # input insar data
    insar_array, gps_array, lons, lats, rms_misfit = compute(gps_los_velfield, xarray, yarray, LOS_array)
    one_to_one_plot(insar_array, gps_array, lons, lats, rms_misfit, plotname, txtname)
    write_output(logname, gps_array, rms_misfit)
    return


def compute(gps_los_velfield, xarray, yarray, LOS_array):
    insar_array, gps_array, lonarray, latarray = los_projection_tools.paired_gps_geocoded_insar(gps_los_velfield,
                                                                                                xarray, yarray,
                                                                                                LOS_array,
                                                                                                window_pixels=10)
    misfit_array = np.subtract(insar_array, gps_array)
    smaller_misfits = np.array([x for x in misfit_array if abs(x) < 15])  # remove the biggest outliers
    rms_misfit = np.sqrt(np.nanmean(smaller_misfits ** 2))
    print("Results: RMS Misfit Between these two fields is %f mm/yr at %d GPS stations \n" % (rms_misfit,
                                                                                              len(insar_array)))
    return insar_array, gps_array, lonarray, latarray, rms_misfit


def one_to_one_plot(insar_array, gps_array, lonarray, latarray, rms_misfit, plotname, txtname):
    plt.figure(figsize=(9, 9), dpi=300)
    plt.plot(gps_array, insar_array, '.', markersize=10)
    bottom_level, top_level = -35, 35
    plt.plot([bottom_level, top_level], [bottom_level, top_level], '--k')
    plt.xlim([bottom_level, top_level])
    plt.ylim([bottom_level, top_level])
    plt.grid(True)
    plt.xlabel('GNSS LOS Velocity (mm/yr)', fontsize=18)
    plt.ylabel('InSAR LOS Velocity (mm/yr)', fontsize=18)
    plt.gca().tick_params(axis='both', labelsize=16)
    plt.title('InSAR vs GNSS Velocities, RMS=%.3fmm/yr' % rms_misfit, fontsize=18)
    plt.savefig(plotname)
    plt.close()

    ofile = open(txtname, 'w')
    ofile.write("# lon lat insar gnss\n")
    for i in range(len(insar_array)):
        ofile.write("%f %f %f %f\n" % (lonarray[i], latarray[i], insar_array[i], gps_array[i]))
    ofile.close()
    return


def write_output(logname, gps_array, rms_misfit):
    ofile = open(logname, 'w')
    ofile.write("Results: RMS Misfit between two fields is %f mm/yr at %d stations \n" % (rms_misfit, len(gps_array)))
    ofile.close()
    return
