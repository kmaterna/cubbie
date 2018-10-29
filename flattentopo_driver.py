# Use to run Marie Pierre's code on lots of interferograms. 
# Need: phase, phasefilt, amp, topo_ra.grd
# Write the data for each intf into roipac format
# Call fortran code
# Make output plots


import numpy as np 
from subprocess import call
import glob, sys
import readbin


# GLOBAL PARAMETERS
nfit=0
ivar=1
alt_ref=100
thresh_amp=0.2
width=661
length=1524
demfile="mend_topo/topo_radar.hgt"
example_rsc="rsc_files/example_sd.int.rsc"

# INPUTS
intf_list=glob.glob("intf_all/???????_???????");
print(intf_list);

# COMPUTE
for data_dir in intf_list:
	intf_name=data_dir.split('/')[1];
	infile=data_dir+"/intf_sd.int";
	infile_filtered=data_dir+"/intf_filt.int";
	stratfile=data_dir+"/strat.unw"
	outfile=data_dir+"/out.int"
	outfile_filtered=data_dir+"/out_filtered.int"

	# MAKE BINARY INTERFEROGRAMS
	[width, length] = readbin.write_gmtsar2roipac_phase(data_dir,infile,infile_filtered);

	# # RUN THE FORTRAN
	print("\nRunning the fortran code to remove atmospheric artifacts from interferogram.")
	call(['cp',example_rsc,data_dir+'/intf_sd.int.rsc'],shell=False);
	print("./flattentopo "+infile+" "+infile_filtered+" "+demfile+" "+outfile+" "+outfile_filtered+" "+str(nfit)+" "+str(ivar)+" "+str(alt_ref)+" "+str(thresh_amp)+" "+stratfile+"\n");
	call(["./flattentopo",infile,infile_filtered,demfile,outfile,outfile_filtered,str(nfit),str(ivar),str(alt_ref),str(thresh_amp),stratfile],shell=False);
	call(['mv','ncycle_topo',data_dir+'/ncycle_topo'],shell=False);
	call(['mv','ncycle_topo_az',data_dir+'/ncycle_topo_az'],shell=False);

	# Output plots
	[real,imag]=readbin.read_binary_roipac_real_imag(infile, width);
	[phase_early,amp_early]=readbin.real_imag2phase_amp(real,imag);
	[real,imag]=readbin.read_binary_roipac_real_imag(outfile, width);
	[phase_late,amp_late]=readbin.real_imag2phase_amp(real,imag);

	readbin.outputs(phase_early, phase_late, width, length, data_dir+"/"+intf_name+"_corrected.eps");


