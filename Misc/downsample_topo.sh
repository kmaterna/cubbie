#!/bin/bash

# June 12, 2019
# I realized that topo_ra.grd is not the same downsamlping as the intf_all/phasefilt.grd files. 
# This might be a problem down the line. 
# We also want to re-do the topo_radar.hgt for the new super-master. 


# Step 1: Downsample topo_ra.grd according to phasefilt.grd
ds_range=`gmt grdinfo intf_all/2019072_2019108/phasefilt.grd -I-`
echo $ds_range
ds_interval=`gmt grdinfo intf_all/2019072_2019108/phasefilt.grd -I`
echo $ds_interval

echo "gmt grdsample topo/topo_ra.grd -Gtopo/temp.grd $ds_interval $ds_range"
gmt grdsample topo/topo_ra.grd -Gtopo/temp.grd $ds_interval $ds_range
# Make sure the -T flag is set or not set to allow the right file size

echo "Converting into a NetCDF3 file."
nccopy -k classic topo/temp.grd topo/topo_ra_subsampled_june.grd

gmt grdinfo topo/topo_ra_subsampled_june.grd


# IN PYTHON:
# For re-doing the topo_radar.hgt given the new supermaster
# Necessary for Marie Pierre's code 
# This is also now put in the first part of the flattentopo_driver.py
#in_topo="topo/topo_ra_subsampled_june.grd"
#out_topo="topo/topo_radar.hgt"
#[width, length]=readbin.write_gmtsar2roipac_topo(in_topo, out_topo);