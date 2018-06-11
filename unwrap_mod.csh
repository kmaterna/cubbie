#!/bin/csh -f
#       $Id$

# generate interferograms for tops stacks
# used for time series analysis

# Xiaohua(Eric) Xu, Jan 20 2016
#

  if ($#argv != 2) then
    echo ""
    echo "Usage: unwrap_tops.csh intf.in batch_tops.config"
    echo "  unwrap a set of tops images in intf.in, dem required in ./topo"
    echo "  supermaster's name required in batch_tops.config"
    echo ""
    echo "  format of data.in:"
    echo "    master_image_stem:slave_image_stem"
    echo ""
    echo "  example of intf.in"
    echo "    S1A20150628_ALL_F1:S1A20150720_ALL_F1"
    echo "    S1A20150720_ALL_F1:S1A20150809_ALL_F1"
    echo ""
    echo "  outputs:"
    echo "    to ./intf_all"
    echo ""
    exit 1
  endif



#
# read parameters from config file
#

  set stage = `grep proc_stage $2 | awk '{print $3}'`
  set master = `grep master_image $2 | awk '{print $3}'`
#
# if filter wavelength is not set then use a default of 200m
#
  set filter = `grep filter_wavelength $2 | awk '{print $3}'`
  if ( "x$filter" == "x" ) then
  set filter = 200
  echo " "
  echo "WARNING filter wavelength was not set in config.txt file"
  echo "        please specify wavelength (e.g., filter_wavelength = 200)"
  echo "        remove filter1 = gauss_alos_200m"
  endif
  set dec = `grep dec_factor $2 | awk '{print $3}'`
  set topo_phase = `grep topo_phase $2 | awk '{print $3}'`
  set shift_topo = `grep shift_topo $2 | awk '{print $3}'`
  set threshold_snaphu = `grep threshold_snaphu $2 | awk '{print $3}'`
  set threshold_geocode = `grep threshold_geocode $2 | awk '{print $3}'`
  set region_cut = `grep region_cut $2 | awk '{print $3}'`
  set switch_land = `grep switch_land $2 | awk '{print $3}'`
  set defomax = `grep defomax $2 | awk '{print $3}'`
  set range_dec = `grep range_dec $2 | awk '{print $3}'`
  set azimuth_dec = `grep azimuth_dec $2 | awk '{print $3}'`




  foreach line (`awk '{print $0}' $1`)
    set ref = `echo $line | awk -F: '{print $1}'`
    set rep = `echo $line | awk -F: '{print $2}'`
    set ref_id  = `grep SC_clock_start ./raw/$ref.PRM | awk '{printf("%d",int($3))}' `
    set rep_id  = `grep SC_clock_start ./raw/$rep.PRM | awk '{printf("%d",int($3))}' `

    #
    # unwrapping
    #

    echo $ref_id"_"$rep_id

    cd intf_all
    cd $ref_id"_"$rep_id

    if ($region_cut == "") then
      set region_cut = `gmt grdinfo phasefilt.grd -I- | cut -c3-20`
      echo $region_cut
    endif
    
    if ($threshold_snaphu != 0 ) then
      echo "will unwrap"
      if ($switch_land == 1) then
        cd ../../topo
        if (! -f landmask_ra.grd) then
          landmask.csh $region_cut
        endif
        cd ../intf_all
        cd $ref_id"_"$rep_id
        ln -s ../../topo/landmask_ra.grd .
      endif

      echo ""
      echo "SNAPHU.CSH - START"
      echo "threshold_snaphu: $threshold_snaphu"
      snaphu_interp.csh $threshold_snaphu $defomax $region_cut  # if you've done writing of grid files in python, use snaphu_interp_mod.csh
      echo "SNAPHU.CSH - END"

    else
      echo ""
      echo "SKIP UNWRAP PHASE"
    endif  

    cd ../..

  end
 # for some reason I need an extra line here. 
 # Ending the file on "end" forces the loop to exit after 1 iteration