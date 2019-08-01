#!/bin/bash

# Reverse-geocode some GPS LOS velocities into RA coordinates
# Run from the processing directory
# Set the other directories up top. 

working_dir="simple_stack"
vel_ll=$working_dir/gps_ll_enu_los.txt
vel_llxyz=$working_dir/gps_ll_los.xyz
outfile=$working_dir/gps_ra_los.xyz

awk {'print $1, $2, $5'} $vel_ll > $vel_llxyz


proj_ll2ra_ascii.csh topo/trans.dat $vel_llxyz $outfile

rm $vel_llxyz