# Use to run Marie Pierre's code on lots of interferograms. 
# Need: phase, phasefilt, amp, topo_ra.grd
# Write the data for each intf into roipac format
# Call fortran code
# Make output plots


import numpy as np
import subprocess
import glob
from read_write_insar_utilities import readbin
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.read_write.netcdf_read_write import read_netcdf3


def main_function():
    # GLOBAL PARAMETERS
    nfit = 0
    ivar = 1
    alt_ref = 100
    thresh_amp = 0.2
    subsampled_dem_grd = "topo/topo_ra_subsampled_june.grd"  # subsampled in the same way as the intf grd files.
    demfile = "topo/topo_radar.hgt"
    example_rsc = "rsc_files/example_sd.int.rsc"

    # INPUTS
    intf_list = glob.glob("intf_all/???????_???????");
    print(intf_list);

    # [width, length] = readbin.write_gmtsar2roipac_topo(subsampled_dem_grd,demfile);
    # copying the demfile into roipac format just in case we've switched the master.

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

        readbin.prep_files(phasefile, phasefilt_file, orig_phasefile,
                           orig_phasefilt_file);  # copy files into safekeeping if necessary

        # MAKE BINARY INTERFEROGRAMS
        [width, length] = readbin.write_gmtsar2roipac_phase(data_dir, orig_phasefile, orig_phasefilt_file, ampfile,
                                                            infile, infile_filtered);
        # width=661; length=1524

        # # # RUN THE FORTRAN
        print("\nRunning the fortran code to remove atmospheric artifacts from interferogram.")
        subprocess.call(['cp', example_rsc, data_dir + '/intf_sd.int.rsc'], shell=False);
        print(
            "flatten_topo " + infile + " " + infile_filtered + " " + demfile + " " + outfile + " " + outfile_filtered
            + " " + str(nfit) + " " + str(ivar) + " " + str(alt_ref) + " " + str(thresh_amp) + " " + stratfile + "\n");
        subprocess.call(
            ["flatten_topo", infile, infile_filtered, demfile, outfile, outfile_filtered, str(nfit), str(ivar),
             str(alt_ref), str(thresh_amp), stratfile], shell=False);
        subprocess.call(['mv', 'ncycle_topo', data_dir + '/ncycle_topo'], shell=False);
        subprocess.call(['mv', 'ncycle_topo_az', data_dir + '/ncycle_topo_az'], shell=False);

        # Output handling. First reading 1D arrays
        [real, imag] = readbin.read_binary_roipac_real_imag(infile);
        [phase_early, amp_early] = readbin.real_imag2phase_amp(real, imag);
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile);
        [phase_out, amp_late] = readbin.real_imag2phase_amp(real, imag);
        [real, imag] = readbin.read_binary_roipac_real_imag(outfile_filtered);
        [phasefilt_out, amp_late] = readbin.real_imag2phase_amp(real, imag);  # 1d arrays

        # # Write the GRD files , fixing an issue with pixel node registration.
        [xdata_p, ydata_p] = read_netcdf3(orig_phasefile)[0:2];
        [xdata_pf, ydata_pf] = read_netcdf3(orig_phasefilt_file)[0:2];
        phase_out_grd = np.reshape(phase_out, (length, width));
        phasefilt_out_grd = np.reshape(phasefilt_out, (length, width));

        netcdf_read_write.produce_output_netcdf(xdata_p, ydata_p, phase_out_grd, 'radians', phasefile);
        netcdf_read_write.flip_if_necessary(phasefile);
        origrange = subprocess.check_output('gmt grdinfo -I- ' + orig_phasefile,
                                            shell=True);  # fixing pixel node registration issue.
        subprocess.call('gmt grdsample ' + phasefile + ' ' + origrange.split()[0] + ' -G' + phasefile + ' -r',
                        shell=True);

        netcdf_read_write.produce_output_netcdf(xdata_pf, ydata_pf, phasefilt_out_grd, 'radians', phasefilt_file);
        netcdf_read_write.flip_if_necessary(phasefilt_file);
        origrange = subprocess.check_output('gmt grdinfo -I- ' + orig_phasefilt_file, shell=True);
        subprocess.call('gmt grdsample ' + phasefilt_file + ' ' + origrange.split()[0] + ' -G' + phasefilt_file + ' -r',
                        shell=True);

        # Making plot
        readbin.output_plots(phase_early, phasefilt_out, width, length, data_dir + "/" + intf_name + "_corrected.eps");

    return;
