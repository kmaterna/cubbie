#!/bin/csh -f
#       $Id$
#
#
#    Xiaohua(Eric) XU, July 7, 2016
#
# Script for merging 3 subswaths TOPS interferograms and then unwrap and geocode. 
# Call this from your main processing directory. 
# Example: 
# quick_geocode.csh stacking/nsbas_apr20/combined merged 20150514.grd 20150514_ll
#
    set grd_directory = $1
    set trans_dat_directory = $2
    set datestr_grd = $3
    set datestr_ll = $4

    echo ""
    echo "GEOCODE-START"

    cp "$trans_dat_directory/trans.dat" $grd_directory
    cd $grd_directory

    proj_ra2ll.csh trans.dat $datestr_grd "$datestr_ll.grd"
    set BT = 26
    set BL = -4
    echo "VARIABLES:"
    echo $BT
    echo $BL
    gmt makecpt -T$BL/$BT/0.5 -Z > unwrap.cpt
    grd2kml.csh $datestr_ll unwrap.cpt
    
    echo "GEOCODE END"