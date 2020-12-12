#!/bin/csh -f
#       $Id$

# generate dem2topo_ra for stacks of tops interferograms
# used for time series analysis
# Kathryn Materna, February 2018. 
# Only change is separating this from intf_tops.csh, and changing 
# scl2amp range dec factor = 2, 
# which matches the range dec factor I've found in dem2topo_ra.csh for sentinel images. 
# Also changing swath directories

# read parameters from config file
#
  set stage = `grep proc_stage $1 | awk '{print $3}'`
  set master = `grep master_image $1 | awk '{print $3}'`
#
# if filter wavelength is not set then use a default of 200m
#
  set filter = `grep filter_wavelength $1 | awk '{print $3}'`
  if ( "x$filter" == "x" ) then
  set filter = 200
  echo " "
  echo "WARNING filter wavelength was not set in config.txt file"
  echo "        please specify wavelength (e.g., filter_wavelength = 200)"
  echo "        remove filter1 = gauss_alos_200m"
  endif
  set dec = `grep dec_factor $1 | awk '{print $3}'`
  set topo_phase = `grep topo_phase $1 | awk '{print $3}'`
  set shift_topo = `grep shift_topo $1 | awk '{print $3}'`
  set threshold_snaphu = `grep threshold_snaphu $1 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $1 | awk '{print $3}'`
  set region_cut = `grep region_cut $1 | awk '{print $3}'`
  set switch_land = `grep switch_land $1 | awk '{print $3}'`
  set defomax = `grep defomax $1 | awk '{print $3}'`
  set swath = `grep swath $1 | awk '{print $3}'` 

##################################
#  start from make topo_ra  #
##################################

#
# clean up
#
  #cleanup.csh topo
#
# make topo_ra
#
if ($topo_phase == 1) then
  echo " "
  echo "DEM2TOPOPHASE.CSH - START"
  echo "USER SHOULD PROVIDE DEM FILE"
  cd F$swath
  cd topo
  cp ../raw/$master.PRM ./master.PRM
  ln -s ../raw/$master.LED .
  if (-f dem.grd) then
    dem2topo_ra.csh master.PRM dem.grd
  else
    echo "no DEM file found: " dem.grd
    exit 1
  endif
  cd ../../
  echo "DEM2TOPOPHASE.CSH - END"

#
#  shift topo_ra
#  
  if ($shift_topo == 1) then
    echo " "
    echo "OFFSET_TOPO - START"
    cd F$swath
    cd topo
    ln -s ../raw/$master.SLC .
    slc2amp.csh master.PRM 2 amp-$master.grd
    offset_topo amp-$master.grd topo_ra.grd 0 0 7 topo_shift.grd
    cd ..
    echo  "OFFSET_TOPO - END"
  else if ($shift_topo == 0) then
    echo "NO TOPOPHASE SHIFT "
  else
    echo "Wrong paramter: shift_topo "$shift_topo
    exit 1
  endif
else if ($topo_phase == 0) then
  echo "NO TOPOPHASE IS SUBSTRACTED"
else
  echo "Wrong paramter: topo_phase "$topo_phase
  exit 1
endif

##################################################


