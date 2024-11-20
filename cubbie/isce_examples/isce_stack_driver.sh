#!/bin/bash
# Running ISCE stack processing
# using isceenv

stackSentinel.py -s ../data/DATA -o ../../../../kmaterna/Documents/S1_orbits/ -a ../../../../kmaterna/Documents/S1_orbits/ -d ../../event_2018_D173/dem/demLat_N32_N34_Lon_W117_W114.dem.wgs84 -b '32.60 33.4 -116.0 -115.0' -c 2 -C geometry -r 5 -z 2 -f 0.3
# -o orbits
# -s SLC
# -d dem file
# -b SNWE bbox. Must be contained within the images' overlapping bbox. 
# -a aux dir
# -C coregistration (geometry coregistration or NESD coregistration; trying geom for Salton Sea)
# -c number of connections
# -m reference_date [default: first date]
# -z AZIMUTHLOOKS [default: 3]
# -r RANGELOOKS [default: 9]
# -f FILTSTRENGTH [default: 0.5]
# -W workflow
# This will print the driver files to be used. 
# The result is a lot of config files and run files.

# Change into the run directory, and then: 
# Run file 01: Unpack reference image. 
# Run file 02: Unpack 7 secondary images. 
# Run file 03: Average parallel/perp Baselines for each swath. Makes n-1 folders, from reference to each secondary acquisition. Very fast step and easy to deal with. 
# Run file 04: geo2rdr. Makes this step 7 times, one for each secondary image. Computer takes off for a few minutes. 
# Run file 05: resample. Takes 20 minutes. 
# Run file 06: extract valid region. Takes almost no time. 
# Run file 07: merge reference secondary slc. This is maybe where multilooking begins. Actually, it took almost no time... didn't do anything? 
# Run file 08: generate burst igram. Complex multiplication seems pretty easy actually. Leaves the igrams divided by burst and by swath. Significantly faster than coregistering them. Probably involves multilooking and filtering? 
# Run file 09: merge burst igram. Starts with fine interferogram, which probably isn't multilooked. 
# Run file 10: filter coherence. Filtering takes a little while. If you don't do this step, maybe you wouldn't get coherence under this workflow? 
# Run file 11: unwrap (optional) 

cd run_files
chmod u+x *run_*

./run_01_unpack_topo_reference
./run_02_unpack_secondary_slc
./run_03_average_baseline
./run_04_fullBurst_geo2rdr
./run_05_fullBurst_resample
./run_06_extract_stack_valid_region
./run_07_merge_reference_secondary_slc
./run_08_generate_burst_igram
./run_09_merge_burst_igram
./run_10_filter_coherence
./run_11_unwrap


