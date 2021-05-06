"""
Run Marie Pierre's code on a general stack of GMTSAR interferograms. Last worked through 2021.
You need:
  -- A series of interferogram directories with phase, phasefilt, amp, corr
  -- topo_ra.grd that covers the whole range of the file in same inc as intf (special scripts from Xiaohua for doing so)
  --    topo_radar.hgt.rsc (just for length metadata; can be copied from example .rsc, same dir as topo_radar.hgt)
  -- An example .rsc file with the width and length fields filled in appropriately (MANUAL)
  -- flatten_topo.f compiled and on your system path
Then, this script does the following:
  -- Writes the data for each intf into roipac format
  -- Calls fortran code
  -- Makes output plots
"""


import numpy as np
import subprocess
import glob, os
from read_write_insar_utilities import readbin
from math_tools import phase_math
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.read_write.netcdf_read_write import read_any_grd


def main_function(intf_directory, flattentopo_directory, topo_ra_file, example_rsc):
    """
    :param intf_directory: location where all your interferograms are stored
    :param flattentopo_directory: location where all your output corrected igrams are stored
    :param topo_ra_file: name of topo_ra.grd, registered with same pixels as intf files
    :param example_rsc: name of rsc file with proper length and width set
    """

    # GLOBAL PARAMETERS
    nfit = 0
    ivar = 1
    alt_ref = 100   # changing this during experiments
    thresh_amp = 0.2   # changing this during experiments
    bin_demfile = flattentopo_directory+"/topo_radar.hgt"          # binary topo file

    # INPUTS
    intf_list = glob.glob(intf_directory + "/???????_???????");
    print(len(intf_list), " interferograms found for flatten_topo");

    # Turning dem into roipac format
    readbin.write_gmtsar2roipac_topo(topo_ra_file, bin_demfile);
    [width, length] = readbin.get_file_shape(topo_ra_file);

    # COMPUTE
    for data_dir in intf_list:
        print(data_dir);
        intf_name = data_dir.split('/')[1];   # is this general or specific to one-level-deep directories?
        outdir = flattentopo_directory + intf_name + '/';
        if os.path.isfile(outdir+'/phase.grd'):
            print("skipping %s " % outdir);
            continue;
        subprocess.call(["mkdir", "-p", outdir], shell=False);

        infile = outdir + "/intf_sd.int";
        infile_filtered = outdir + "/intf_filt.int";
        stratfile = outdir + "/strat.unw"
        outfile = outdir + "/out.int"
        outfile_filtered = outdir + "/out_filtered.int"

        # GMTSAR files.
        orig_phasefile = data_dir + "/phase.grd";
        orig_phasefilt_file = data_dir + "/phasefilt.grd";
        orig_ampfile = data_dir + "/amp.grd";
        orig_corrfile = data_dir + "/corr.grd";
        orig_maskfile = data_dir + "/mask.grd";
        out_phasefile = outdir + "/phase.grd";
        out_phasefilt_file = outdir + "/phasefilt.grd";
        out_ampfile = outdir + "/amp.grd";
        out_corrfile = outdir + "/corr.grd";
        out_maskfile = outdir + "/mask.grd";

        # MAKE BINARY INTERFEROGRAMS
        readbin.write_gmtsar2roipac_phase(orig_phasefile, orig_phasefilt_file, orig_ampfile, infile, infile_filtered);

        # # # RUN THE FORTRAN
        print("\nRunning the fortran code to remove atmospheric artifacts from interferogram.")
        subprocess.call(['cp', example_rsc, outdir + '/intf_sd.int.rsc'], shell=False);
        print(
            "flattentopo " + infile + " " + infile_filtered + " " + bin_demfile + " " + outfile + " " + outfile_filtered
            + " " + str(nfit) + " " + str(ivar) + " " + str(alt_ref) + " " + str(thresh_amp) + " " + stratfile + "\n");
        subprocess.call(
            ["flattentopo", infile, infile_filtered, bin_demfile, outfile, outfile_filtered, str(nfit), str(ivar),
             str(alt_ref), str(thresh_amp), stratfile], shell=False);
        subprocess.call(['mv', 'ncycle_topo', outdir + '/ncycle_topo'], shell=False);
        subprocess.call(['mv', 'ncycle_topo_az', outdir + '/ncycle_topo_az'], shell=False);
        subprocess.call(['cp', orig_ampfile, out_ampfile], shell=False);
        subprocess.call(['cp', orig_corrfile, out_corrfile], shell=False);
        subprocess.call(['cp', orig_maskfile, out_maskfile], shell=False);

        # Output handling. First reading 1D arrays
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile);
        [phase_out, _] = phase_math.real_imag2phase_amp(real, imag);
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile_filtered);
        [phasefilt_out, _] = phase_math.real_imag2phase_amp(real, imag);  # 1d arrays

        # Reshape grids into two-dimensional arrays
        phase_out_grd = np.reshape(phase_out, (length, width));
        phasefilt_out_grd = np.reshape(phasefilt_out, (length, width));

        # Write GRD files of output quantities
        [xdata_p, ydata_p] = read_any_grd(orig_phasefile)[0:2];
        [xdata_pf, ydata_pf, phasefilt_early] = read_any_grd(orig_phasefilt_file);
        netcdf_read_write.produce_output_netcdf(xdata_p, ydata_p, phase_out_grd, 'radians', out_phasefile);
        netcdf_read_write.produce_output_netcdf(xdata_pf, ydata_pf, phasefilt_out_grd, 'radians', out_phasefilt_file);

        # Making plot
        readbin.output_plots(phasefilt_early, phasefilt_out, width, length,
                             outdir + "/" + intf_name + "_corrected.eps");

    return;
