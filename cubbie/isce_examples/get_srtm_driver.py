#!/bin/bash
# use ISCE to download DEM from SRTM

# EASY WAY:
dem.py -a stitch -b 32 34 -117 -114 -r -s 1 -c
# -r means report 
# -s 1 means srtm 1 



# Way that didn't work: isce stack processor couldn't find the .vrt somehow. 
# Create dummy files: reference.xml, secondary.xml, and topsApp.xml
# Point to the right files, two examples from the Westmorland stack. 
# Then call up to --dostep verifyDEM from here. 
# topsApp.py --start=startup --end=verifyDEM. 
# DON'T DO THIS ONE. 

