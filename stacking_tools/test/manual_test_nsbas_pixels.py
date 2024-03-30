# October 2020, testing pixels

import numpy as np
import matplotlib.pyplot as plt
from . import io_functions
from .. import nsbas, dem_error_correction
from s1_batches.intf_generating import sentinel_utilities


def outputs(x_axis_days, ts):
    plt.figure()
    plt.plot(x_axis_days, ts, '.')
    plt.savefig('Testing_Results/test.png')
    vel = nsbas.compute_velocity_math(ts, x_axis_days)
    print("Velocity is %.12f mm/yr" % vel)
    return


def outputs_dem_error_correction(x_axis_days, ts_raw, ts_corrected):
    plt.figure()
    plt.plot(x_axis_days, ts_raw, '--.r', label='uncorrected')
    plt.plot(x_axis_days, ts_corrected, '--.k', label='corrected')
    plt.xlabel('Days')
    plt.ylabel('Displacement (mm)')
    plt.legend()
    plt.savefig("Testing_Results/DEM_error_results.png")
    vel = nsbas.compute_velocity_math(ts_raw, x_axis_days)
    print(" Original velocity is %.12f mm/yr" % vel)
    vel = nsbas.compute_velocity_math(ts_corrected, x_axis_days)
    print("Corrected velocity is %.12f mm/yr" % vel)
    return


def test_single_pixel_nsbas():
    # # Testing a single pixel, regular NSBAS with weights etc.
    # You can always call this to test a single pixel.
    ifile = 'Testing_Data/testing_pixel_11.txt'
    Igrams = io_functions.read_testing_pixel(ifile, coherence=True)
    Igrams = io_functions.take_coherent_igrams(Igrams, 0.375)
    ts = nsbas.do_nsbas_pixel(Igrams.phase, Igrams.juldays, 56, Igrams.datestrs, coh_value=Igrams.corr)
    outputs(Igrams.x_axis_days, ts)
    return


def test_real_pixel_dem_error():
    # Testing the baseline correction of Fattahi and Amelung, 2013 on an example pixel
    ifile = 'Testing_Data/testing_pixel_3.txt'
    baseline_table = 'Testing_Data/baseline_table.dat'
    baseline_tuple = sentinel_utilities.read_baseline_table(baseline_table)
    Igrams = io_functions.read_testing_pixel(ifile, coherence=False)
    ts = nsbas.do_nsbas_pixel(Igrams.phase, Igrams.juldays, 56, Igrams.datestrs)
    ts_corrected, K_z_error = dem_error_correction.driver(ts, Igrams.datestrs, baseline_tuple)
    print("DEM ERROR: %.10f*rsintheta m " % K_z_error)
    outputs_dem_error_correction(Igrams.x_axis_days, ts, ts_corrected)
    # [stems, times, baselines, missiondays] = sentinel_utilities.read_baseline_table(baseline_table)
    # sentinel_utilities.make_network_plot(Igrams.juldays, stems, times, baselines,
    # "Testing_Results/pixel_baseline_plot.png")
    return


def test_synthetic_pixel_dem_error():
    """
    Figure 1 of Fattahi and Amelung, 2013
    A synthetic test that isolates the DEM error correction.
    """
    baseline_table = 'Testing_Data/Fattahi_Fig1/baseline_table.dat'
    pixel_file = 'Testing_Data/Fattahi_Fig1/pixel_value.txt'
    baseline_tuple = sentinel_utilities.read_baseline_table(baseline_table)
    [datestrs, ts] = np.loadtxt(pixel_file, unpack=True, dtype={'names': ('a', 'b'), 'formats': ('U7', np.float)})
    ts_corrected, K_z_error = dem_error_correction.driver(ts, datestrs, baseline_tuple)
    print("DEM ERROR: %.10f*rsintheta m " % K_z_error)
    outputs_dem_error_correction(range(0, 17 * 8, 8), ts, ts_corrected)
    return


if __name__ == "__main__":
    # test_real_pixel_dem_error()
    test_synthetic_pixel_dem_error()
