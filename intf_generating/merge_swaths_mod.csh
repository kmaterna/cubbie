#!/bin/csh -f
#       $Id$
#
#
#    Xiaohua(Eric) XU, July 7, 2016
#
# Script for merging 3 subswaths TOPS interferograms for a stack of interferograms.
#
#

  if ($#argv != 2) then
    echo ""
    echo "Usage: merge_batch.csh inputfile config_file"
    echo ""
    echo "Note: Inputfiles should be as following:"
    echo ""
    echo "      IF1_Swath1_Path:master.PRM:repeat.PRM,IF1_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
    echo "      IF2_Swath1_Path:master.PRM:repeat.PRM,IF2_Swath2_Path:master.PRM:repeat.PRM,IF1_Swath3_Path:master.PRM:repeat.PRM"
    echo "      (Use the repeat PRM which contains the shift information.)"
    echo "      e.g. ../F1/intf_all/2015092_2015128/:S1A20150403_ALL_F1.PRM:S1A20150509_ALL_F1.PRM,../F2/intf_all/2015092_2015128/:S1A20150403_ALL_F2.PRM:S1A20150509_ALL_F2.PRM,../F3/intf_all/2015092_2015128/:S1A20150403_ALL_F3.PRM:S1A20150509_ALL_F3.PRM"
    echo ""
    echo "      Make sure under each path, the processed phasefilt.grd, corr.grd and mask.grd exist."
    echo "      Also make sure the dem.grd is linked. "
    echo "      If trans.dat exits, recomputation of projection matrix will not proceed."
    echo "      The master image of firet line should be the super_master."
    echo "      KZM: CALL FROM UPPER LEVEL PROCESSING DIRECTORY, not merged."

    echo ""
    echo "      config_file is the same one used for processing."
    echo ""
    echo "Example: merge_batch.csh filelist batch.config"
    echo ""
    exit 1
  endif

  cd merged/

  if (! -f dem.grd) then
    echo "dem.grd is required ..."
    exit 1
  endif

  set input_file = $1
  awk 'NR==1{print $0}' $input_file | awk -F, '{for (i=1;i<=NF;i++) print "../"$i}' | awk -F: '{print $1$2}'> tmpm.filelist

  set now_dir = `pwd`
  set merged_top_dir = `pwd`

  foreach line (`awk '{print $0}' $input_file`)
    set dir_name = `echo $line | awk -F, '{print $1}' | awk -F: '{print $1}' | awk -F"/" '{print $(NF-1)}'`
    mkdir $dir_name
    cd $dir_name
    echo $line | awk -F, '{for (i=1;i<=NF;i++) print "../"$i}' > tmp.filelist
    paste ../tmpm.filelist tmp.filelist | awk '{print $1","$2}' > tmp
    rm tmp.filelist

    foreach f_name (`awk '{print $0}' < tmp`)
        set mm = `echo $f_name | awk -F, '{print $1}'`
        set pth = `echo $f_name | awk -F, '{print $2}' | awk -F: '{print $1}'`
        set f1 = `echo $f_name | awk -F, '{print $2}' | awk -F: '{print $2}'`
        set f2 = `echo $f_name | awk -F, '{print $2}' | awk -F: '{print $3}'`
        cp $mm ./supermaster.PRM
        set rshift = `grep rshift $pth$f1 | tail -1 | awk '{print $3}'`
        update_PRM supermaster.PRM rshift $rshift
        set fs1 = `grep first_sample supermaster.PRM | awk '{print $3}'`
        set fs2 = `grep first_sample $pth$f1 | awk '{print $3}'`
        if ($fs2 > $fs1) then
          update_PRM supermaster.PRM first_sample $fs2
        endif
        cp supermaster.PRM $pth
        echo $pth":supermaster.PRM:"$f2 >> tmp.filelist
    end

    if (-f ../trans.dat) ln -s ../trans.dat .
    if (-f ../raln.grd) ln -s ../raln.grd .
    if (-f ../ralt.grd) ln -s ../ralt.grd .
    if (-f ../landmask_ra.grd ) ln -s ../landmask_ra.grd .
    ln -s ../dem.grd .
    ln -s ../$2 .
    rm tmp


    # Copying the guts of merge_unwrap_geocode_tops.csh.
    # EQUIVALENT TO: merge_unwrap_geocode_tops.csh tmp.filelist $2

    if (-f tmp_phaselist) rm tmp_phaselist
    if (-f tmp_corrlist) rm tmp_corrlist
    if (-f tmp_masklist) rm tmp_masklist

    if (! -f dem.grd ) then
      echo "Please link dem.grd to current folder"
      exit 1
    endif

    set region_cut = `grep region_cut $2 | awk '{print $3}'`

    # Creating inputfiles for merging
    foreach line (`awk '{print $0}' tmp.filelist`)
      set now_dir = `pwd`
      set pth = `echo $line | awk -F: '{print $1}'`
      set prm = `echo $line | awk -F: '{print $2}'`
      set prm2 = `echo $line | awk -F: '{print $3}'`
      cd $pth
      set rshift = `grep rshift $prm2 | tail -1 | awk '{print $3}'`
      set fs1 = `grep first_sample $prm | awk '{print $3}'`
      set fs2 = `grep first_sample $prm2 | awk '{print $3}'`
      cp $prm tmp.PRM
      if ($fs2 > $fs1) then
        update_PRM tmp.PRM first_sample $fs2
      endif
      update_PRM tmp.PRM rshift $rshift
      cd $now_dir

      echo $pth"tmp.PRM:"$pth"phasefilt.grd" >> tmp_phaselist
      echo $pth"tmp.PRM:"$pth"corr.grd" >> tmp_corrlist
      echo $pth"tmp.PRM:"$pth"mask.grd" >> tmp_masklist
    end

    set pth = `awk -F: 'NR==1 {print $1}' tmp.filelist`
    set stem = `awk -F: 'NR==1 {print $2}' tmp.filelist | awk -F"." '{print $1}'`
    #echo $pth $stem

    echo ""
    echo "Merging START"
    merge_swath tmp_phaselist phasefilt.grd $stem
    merge_swath tmp_corrlist corr.grd
    merge_swath tmp_masklist mask.grd
    echo "Merging END"
    echo ""

    set scale = "-JX6.5i"
    gmt makecpt -Crainbow -T-3.15/3.15/0.1 -Z -N > phase.cpt
    gmt grdimage phasefilt.grd $scale -Bxaf+lRange -Byaf+lAzimuth -BWSen -Cphase.cpt -X1.3i -Y3i -P -K > phasefilt.ps
    gmt psscale -Rphasefilt.grd -J -DJTC+w5i/0.2i+h -Cphase.cpt -B1.57+l"Phase" -By+lrad -O >> phasefilt.ps
    gmt psconvert -Tf -P -Z phasefilt.ps
    gmt makecpt -T0./.8/0.1 -Cgray -Z -N > corr.cpt
    gmt grdimage corr.grd $scale -Bxaf+lRange -Byaf+lAzimuth -BWSen -Ccorr.cpt -X1.3i -Y3i -P -K > corr.ps
    gmt psscale -Rcorr.grd -J -DJTC+w5i/0.2i+h -Ccorr.cpt -B1.57+l"Phase" -By+lrad -O >> corr.ps
    gmt psconvert -Tf -P -Z corr.ps

    rm tmp_phaselist
    rm tmp_corrlist
    rm tmp_masklist
    rm phase.cpt
    rm corr.cpt
    rm tmp_phaselist tmp_corrlist tmp_masklist *.eps *.bb


    if (! -f ../trans.dat && -f trans.dat) then
      mv trans.dat ../
      ln -s ../trans.dat .
    endif
    if (! -f ../landmask_ra.grd && -f landmask_ra.grd ) then
      mv landmask_ra.grd  ../
      ln -s ../landmask_ra.grd .
    endif
    if (! -f ../raln.grd && -f raln.grd) then
      mv raln.grd ../
      ln -s ../raln.grd .
    endif
    if (! -f ../ralt.grd && -f raln.grd) then
      mv ralt.grd ../
      ln -s ../ralt.grd .
    endif

    cd $merged_top_dir

  end
  cd ../
