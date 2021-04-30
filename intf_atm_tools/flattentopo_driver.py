"""
Run Marie Pierre's code on a general stack of GMTSAR interferograms. Last worked through 2021.
You need:
  -- A series of interferogram directories with phase, phasefilt, amp, corr
  -- topo_ra.grd that covers the whole range of the file (for merged files, I've used dem2topo_ra.csh separately)
  -- ex: dem2topo_ra.csh master.PRM dem.grd xmin/xmax/ymin/ymax,
  --    for this example, range set MANUALLY during merged workflow
  --    for this example, master.PRM must be copied MANUALLY into directory
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
import glob
from read_write_insar_utilities import readbin
from math_tools import phase_math
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.read_write.netcdf_read_write import read_any_grd


def main_function(intf_directory, topo_ra_file, example_rsc):
    """
    :param intf_directory: location where all your interferograms are stored
    :param topo_ra_file: name of topo_ra.grd
    :param example_rsc: name of rsc file with proper length and width set
    """

    # GLOBAL PARAMETERS
    nfit = 0
    ivar = 1
    alt_ref = 100
    thresh_amp = 0.2
    reg_demfile = intf_directory+"topo_ra_registered.grd"   # topo_ra registered to data grid
    bin_demfile = intf_directory+"topo_radar.hgt"          # binary topo file

    # INPUTS
    intf_list = glob.glob(intf_directory + "/???????_???????");
    print(len(intf_list), " interferograms found for flatten_topo");

    # Turning topo_ra into same range/inc as rest of data, and copying dem into roipac format
    subprocess.call('gmt grdsample '+topo_ra_file+' -G'+reg_demfile+' `gmt grdinfo -I ' +
                    intf_list[0]+'/phase.grd` `gmt grdinfo -I- '+intf_list[0]+'/phase.grd`', shell=True);
    readbin.write_gmtsar2roipac_topo(reg_demfile, bin_demfile);
    [width, length] = readbin.get_file_shape(reg_demfile);

    # COMPUTE
    for data_dir in intf_list:
        print(data_dir);
        intf_name = data_dir.split('/')[1];
        infile = data_dir + "/intf_sd.int";
        infile_filtered = data_dir + "/intf_filt.int";
        stratfile = data_dir + "/strat.unw"
        outfile = data_dir + "/out.int"
        outfile_filtered = data_dir + "/out_filtered.int"

        # Pretty standard GMTSAR files.
        phasefile = data_dir + "/phase.grd";
        phasefilt_file = data_dir + "/phasefilt.grd";
        ampfile = data_dir + "/amp.grd";
        orig_phasefile = data_dir + "/orig_phase.grd";
        orig_phasefilt_file = data_dir + "/orig_phasefilt.grd";

        readbin.prep_files(phasefile, phasefilt_file, orig_phasefile, orig_phasefilt_file);  # save orig files

        # MAKE BINARY INTERFEROGRAMS
        readbin.write_gmtsar2roipac_phase(orig_phasefile, orig_phasefilt_file, ampfile, infile, infile_filtered);

        # # # RUN THE FORTRAN
        print("\nRunning the fortran code to remove atmospheric artifacts from interferogram.")
        subprocess.call(['cp', example_rsc, data_dir + '/intf_sd.int.rsc'], shell=False);
        print(
            "flattentopo " + infile + " " + infile_filtered + " " + bin_demfile + " " + outfile + " " + outfile_filtered
            + " " + str(nfit) + " " + str(ivar) + " " + str(alt_ref) + " " + str(thresh_amp) + " " + stratfile + "\n");
        subprocess.call(
            ["flattentopo", infile, infile_filtered, bin_demfile, outfile, outfile_filtered, str(nfit), str(ivar),
             str(alt_ref), str(thresh_amp), stratfile], shell=False);
        subprocess.call(['mv', 'ncycle_topo', data_dir + '/ncycle_topo'], shell=False);
        subprocess.call(['mv', 'ncycle_topo_az', data_dir + '/ncycle_topo_az'], shell=False);

        # Output handling. First reading 1D arrays
        [real, imag] = readbin.read_binary_roipac_real_imag(infile);
        [phase_early, _] = phase_math.real_imag2phase_amp(real, imag);
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile);
        [phase_out, _] = phase_math.real_imag2phase_amp(real, imag);
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile_filtered);
        [phasefilt_out, _] = phase_math.real_imag2phase_amp(real, imag);  # 1d arrays

        # Reshape grids into two-dimensional arrays
        phase_out_grd = np.reshape(phase_out, (length, width));
        phasefilt_out_grd = np.reshape(phasefilt_out, (length, width));

        # Write GRD files of output quantities
        [xdata_p, ydata_p] = read_any_grd(orig_phasefile)[0:2];
        [xdata_pf, ydata_pf] = read_any_grd(orig_phasefilt_file)[0:2];
        netcdf_read_write.produce_output_netcdf(xdata_p, ydata_p, phase_out_grd, 'radians', phasefile);
        netcdf_read_write.produce_output_netcdf(xdata_pf, ydata_pf, phasefilt_out_grd, 'radians', phasefilt_file);

        # Making plot
        readbin.output_plots(phase_early, phasefilt_out, width, length, data_dir + "/" + intf_name + "_corrected.eps");

    return;
