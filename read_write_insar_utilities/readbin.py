#!/bin/env/python
import numpy as np
import matplotlib.pyplot as plt
import math, os
import subprocess
import struct
from read_write_insar_utilities import netcdf_read_write
from math_tools import phase_math


def read_binary_roipac_real_imag(filename):
    # Reads a binary file that expects two fields (such as real and imaginary)
    # returns one-dimensional arrays.
    with open(filename, mode='rb') as file:  # b is important -> binary
        fileContent = file.read()
    data = struct.unpack("f" * (len(fileContent) // 4), fileContent);
    print(type(fileContent));
    print("data has %d floats" % (len(data) / 2));

    real = [];
    imag = [];
    for x in range(len(data)):
        if np.mod(x, 2) == 0:
            real.append(data[x]);
        else:
            imag.append(data[x]);
    return [real, imag];


def read_binary_topo(filename, width):
    # Reads a topo file and returns 1d arrays for topography and some other quantity.
    with open(filename, mode='rb') as file:  # b is important -> binary
        fileContent = file.read()
    data = struct.unpack("f" * (len(fileContent) // 4), fileContent);
    num_per_array = int(len(data) / 2);
    print("data has %d floats" % num_per_array);

    topo1 = [];
    topo2 = [];
    for x in range(len(data)):
        if np.mod(x, 2 * width) < width:  # topography is organized a different way.
            topo1.append(data[x]);
        else:
            topo2.append(data[x]);

    print("max of topo1 is %d" % (np.nanmax(topo1)));
    print("min of topo1 is %d" % (np.nanmin(topo1)));
    print("max of topo2 is %d" % (np.nanmax(topo2)));
    print("min of topo2 is %d" % (np.nanmin(topo2)));
    return [topo1, topo2];


# Functions I will write tomorrow. 
def write_binary_roipac_real_imag(filename, real, imag):
    # Takes 1D arrays and writes the 1D binary file.
    total_struct = [];
    for i in range(len(real)):
        if not math.isnan(real[i]):  # Removing bad nans from the stack.
            total_struct.append(real[i]);
        else:
            total_struct.append(0.0);
        if not math.isnan(imag[i]):
            total_struct.append(imag[i]);
        else:
            total_struct.append(0.0);
    data = struct.pack("f" * len(real) * 2, *total_struct);
    print("Packing %d real and imaginary numbers into binary array. " % (len(total_struct) / 2));
    with open(filename, mode='wb') as file:
        file.write(data);
    return;


def write_binary_topo(filename, topo1, topo2, width):
    # Takes 1D array and writes the 1D binary file.
    total_struct = [];
    for i in range(len(topo1)):
        if np.mod(i, width) == 0:
            for j in range(i, i + width):
                total_struct.append(topo1[j]);
            for j in range(i, i + width):
                total_struct.append(topo2[j]);

    data = struct.pack("f" * len(topo2) * 2, *total_struct);
    print("Packing %d topo numbers into binary array. " % (len(total_struct) / 2));
    with open(filename, mode='wb') as file:
        file.write(data);
    return;


def write_rsc_file():
    return;


# Plotting
def output_plots(phase, amp, width, length, plotname):
    # Takes 1D arrays and plots them in graphical format.
    phase = np.reshape(phase, (length, width));
    amp = np.reshape(amp, (length, width));

    f, axarr = plt.subplots(1, 2, figsize=(10, 10));
    axarr[0].imshow(phase, cmap='jet');
    axarr[1].imshow(amp, cmap='jet');
    plt.savefig(plotname);
    plt.close();
    return;


def write_gmtsar2roipac_phase(input_directory, phasefile, phasefilt_file, ampfile, outfilename, outfile_filt):
    """
    A function that reads phase and amplitude grids in GMT format and writes a Binary format file
    with the real and imaginary components of the values.
    """

    # INPUTS
    [_, _, phase] = netcdf_read_write.read_grd_xyz(phasefile);
    [_, _, phasefilt] = netcdf_read_write.read_grd_xyz(phasefilt_file);
    [xdata, ydata, amp] = netcdf_read_write.read_grd_xyz(ampfile);
    width = len(xdata);
    length = len(ydata);

    # # FOR THE PHASE AND AMPLITUDE DATA, reformatting the data and making initial plot.
    phase_1d = np.reshape(phase, (np.size(phase),));
    phasefilt_1d = np.reshape(phasefilt, (np.size(phasefilt),));
    amp_1d = np.reshape(amp, (np.size(amp),));

    amp_1d = [x * 1e12 for x in amp_1d];  # Fixing the different GMTSAR definition of amplitude.

    # # WRITE THEM OUT IN BINARY FORMAT.
    print("converting phase_1d to real,imag.");
    [real, imag] = phase_math.phase_amp2real_imag(phase_1d, amp_1d);
    write_binary_roipac_real_imag(outfilename, real, imag);

    print("converting phase_1d_filt to real,imag.");
    [real, imag] = phase_math.phase_amp2real_imag(phasefilt_1d, amp_1d);
    write_binary_roipac_real_imag(outfile_filt, real, imag);
    return [width, length];


def write_gmtsar2roipac_topo(infile, out_topo):
    """
    A function that reads a GMT topographic grid and writes the corresponding .hgt format
    for use with roipac functions.
    """
    [xdata, ydata, topo] = netcdf_read_write.read_grd_xyz(infile);
    width = len(xdata);
    length = len(ydata);
    topo = np.flipud(topo);  # formatting correct when you flip up/down.

    topo_1d = np.reshape(topo, (np.size(topo),));

    # WRITE THE TOPO OUT IN BINARY FORMAT
    write_binary_topo(out_topo, topo_1d, topo_1d, width);
    # outputs(topo_1d, topo_1d, width, length, 'mendocino_topo_orig.eps');
    return [width, length];


def prep_files(phasefile, phasefilt_file, orig_phasefile, orig_phasefilt_file):
    # For new interferograms, you need to save a copy of phase.grd and phasefilt.grd because they're going to be
    # overwritten when we're doing Marie-Pierre's insar workflow.
    if not os.path.isfile(orig_phasefile):
        subprocess.call(['cp', phasefile, orig_phasefile], shell=False);
    if not os.path.isfile(orig_phasefilt_file):
        subprocess.call(['cp', phasefilt_file, orig_phasefilt_file], shell=False);
    return;


if __name__ == "__main__":
    # THE MAIN PROGRAM
    # filename="tibet_intf/20160510-20160721_sd_16rlks.int"; length=1732; width=1067;
    # filename="tibet_intf/20160510-20160721_sd_16rlks.int"; length=1732; width=1067;
    # [real,imag]=read_binary_roipac_real_imag(filename, width);
    # [phase,amp]=real_imag2phase_amp(real,imag);
    # outputs(phase, amp, width, length, 'output_mpd.eps');

    # filename="2015225_2015249/2015225_2015249_sd.int"; length=1524; width=661;
    # [real,imag]=read_binary_roipac_real_imag(filename, width);
    # [phase,amp]=real_imag2phase_amp(real,imag);
    # outputs(phase, amp, width, length, '2015225_2015249/input_phase_amp.eps');

    filename = "2015225_2015249/out.int";
    length = 1524;
    width = 661;
    [real, imag] = read_binary_roipac_real_imag(filename);
    [phase, amp] = phase_math.real_imag2phase_amp(real, imag);
    output_plots(phase, amp, width, length, '2015225_2015249/working_example.eps');

    filename = "2015225_2015249/out_filtered.int";
    length = 1524;
    width = 661;
    [real, imag] = read_binary_roipac_real_imag(filename);
    [phase, amp] = phase_math.real_imag2phase_amp(real, imag);
    output_plots(phase, amp, width, length, '2015225_2015249/working_example_filt.eps');

# Write the file out.
# outfile_name="tibet_intf/written_by_me.int";
# [real,imag]=phase_amp2real_imag(phase,amp);
# write_binary_roipac_real_imag(outfile_name, real, imag);

# WHAT DOES THE REAL/IMAG FILE THAT I WROTE LOOK LIKE?
# filename="tibet_intf/written_by_me.int"; length=1732; width=1067;
# [real,imag]=read_binary_roipac_real_imag(filename, width);
# [phase,amp]=real_imag2phase_amp(real,imag);
# outputs(phase, amp, width, length, 'written_by_me.eps');
# YES IT LOOKS LIKE THE ORIGINAL ONE!!!

# # # READ THE TOPOGRAPHY
# filename="raw_topo/radar_16rlks.hgt"; length=1733; width=1067;
# [topo1, topo2]=read_binary_topo(filename, width);
# outputs(topo1, topo2, width, length,'testimage_topo.eps');

# # # # WRITE THE TOPOGRAPHY
# outfilename="raw_topo/topo.hgt";
# write_binary_topo(outfilename,topo1,topo2,width);
