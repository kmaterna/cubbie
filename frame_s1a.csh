#!/bin/csh -f
#       $Id$
#
# Xiaohua(Eric) Xu, Mar 20 2017
#

  if ($#argv != 4) then
    echo ""
    echo "Usage: frame_s1a.csh filelist pins.ll mode dir_name"
    echo "  organize one track of S1A TOPS data, redefine frames, precise/restituted orbit is required"
    echo ""
    echo "filest:"
    echo "    pth_filename1"
    echo "    pth_filename2"
    echo "    ......"
    echo ""
    echo "pins.ll:"
    echo "    lon1 lat1"
    echo "    lon2 lat2"
    echo "    ......"
    echo ""
    echo "Note: "
    echo "    files listed in filelist should be the .SAFE directory with absolute path."
    echo "    mode = 1 will tell how many records are gonna be generated. mode = 2 will do the organizing."
    echo ""
    exit 1
  endif



  set ii = 0
  set mode = $3
  set direname = $4   # maybe it wasn't so great for someone to name a variable dirname? 
  set called=($_)
  if ( "$called" != "" ) then ### called by source
    set script_fn = `greadlink -f $called[1]`
  else 
    set script_fn = `greadlink -f $0`
  endif
  set code_path = `dirname $script_fn`  # the place where the scripts all live. 
  echo "code path:"
  echo $code_path

  if (-f frame.list) rm -f frame.list
  if (-f tmprecord)  rm -f tmprecord

  # divide the list of files into sets, and create frames based on the given pins
  foreach line (`awk '{print $0}' $1`)
    set file1 = `echo $line | awk -F"," '{print $1}'`  # a full path to a .SAFE directory
    set date1 = `echo $file1 | awk '{print substr($1,length($1)-54,8)}'`   # something like 20171201
    set SAT1 = `echo $file1 | awk '{print substr($1,length($1)-71,3)}'`    # something like S1A
    echo $file1
    echo $date1

    if ($ii == 0) then
      set file0 = `echo $file1`
      set date0 = `echo $date1`
      set SAT0 = `echo $SAT1`
      echo $file1 > tmprecord
      set ii = 1
    else
      if ($date1 == $date0 && $SAT1 == $SAT0) then
        echo $file1 >> tmprecord
      else

        echo "" | awk '{printf("%s ","Combing")}' 
        set jj = 1
        set t2 = 9999999999
        if (-f tmprecord_new) rm -f tmprecord_new

        foreach line2 (`awk '{print $0}' tmprecord`)
          echo $line2 | awk '{printf("%s ",$1)}'
          set tt = `echo $line2 | awk '{print substr($1,length($1)-54,15)}'`
       #   set t1 = `gdate -jf "%Y%m%dT%H%M%S" $tt +%s`
         set ss1 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'` 
          set t1 = `gdate --date="$ss1" +%s`  
         if ($t2 == "9999999999") then
             echo $line2 >> tmprecord_new
         else
           if ($t1 <= $t2) then
             echo $line2 >> tmprecord_new
           endif
         endif
   
          set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
        set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
#          set t2 = `gdate -jf "%Y%m%dT%H%M%S" $tt +%s`
          set t2 = `gdate --date="$ss2" +%s`
        end
        echo "" | awk '{printf("%s\n","...")}'

        set n1 = ` gdate --date="$date0 - 1 day" +%Y%m%d `
        set n2 = ` gdate --date="$date0 + 1 day" +%Y%m%d `
        if ($SAT0 == "S1A") then
            python $code_path/get_s1_orbits.py $date0 s1a
            #get_s1a_orbit.csh $date0
            set orb = ` find . -name "$SAT0*$n1*$n2*.EOF" ` # the name of the file we just copied
        else
            python $code_path/get_s1_orbits.py $date0 s1b
            #get_s1b_orbit.csh $date0
            set orb = ` find . -name "$SAT0*$n1*$n2*.EOF" ` # the name of the file we just copied
        endif


        set pin1 = `head -1 $2` 
        set f1 = `head -1 tmprecord_new`
        echo $f1
        make_s1a_tops $f1/annotation/*iw1*vv*xml $f1/measurement/*iw1*vv*tiff tmp2 0  # seems to go okay? Makes a PRM file with 53 lines. 
        ext_orb_s1a tmp2.PRM $orb tmp2   # I think this works okay as well. 
        set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
        # refinie the calculation
        echo "Starting shift_atime_PRM"
        shift_atime_PRM.csh tmp2.PRM $tmpazi  #
        set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`
        
        set pin2 = `tail -1 $2`
        set f2 = `tail -1 tmprecord_new`
        make_s1a_tops $f2/annotation/*iw1*vv*xml $f2/measurement/*iw1*vv*tiff tmp2 0
        ext_orb_s1a tmp2.PRM $orb tmp2
        set tmpazi = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
        shift_atime_PRM.csh tmp2.PRM $tmpazi   # 
        set azi2 = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`

        set nl = `grep num_lines tmp2.PRM | awk '{print $3}'`  # This was 12192 when I printed it. 
        echo $nl
        echo "We are past the shift_PRM step. "  # we get to here, and I think it's successful. 

        # It looks like this block of code works!
        echo $azi1
        echo $azi2
        if ($azi1 > 0 && $azi2 < $nl ) then  
          awk '{print $1","$2}' $2 > tmpllt
          set pin0 = `awk NR==1'{print $0}' tmpllt`
          foreach line2 (`awk '{print $0}' tmpllt`)
            if ($line2 != $pin0) then
              echo $pin0 | awk -F"," '{print $1,$2}' > tmp1llt
              echo $line2 | awk -F"," '{print $1,$2}' >> tmp1llt
              if ($mode != 1) then
                echo "we are inside mode2"  #  
                create_frame_tops.csh tmprecord_new $orb tmp1llt 1
                mkdir -p $direname
                mv *.SAFE $direname
              else
                echo ""
                echo "Frames on date $date0 will be re-organized..."
                echo ""
              endif
              set pin0 = `echo $line2`
            endif
          end 
         
          echo $date0 >> dates.merge.$direname
        else
          echo $date0 >> dates.skip.$direname
          if ($jj == 0) then
            echo ""
            echo "SKIP $date0, as it stopped observation in the middle ..."
            echo ""
          else
            echo ""
            echo "SKIP $date0, as it does not have enough scenes ..."
            echo ""
          endif
        endif

        echo $file1 > tmprecord
        set file0 = `echo $file1`
        set date0 = `echo $date1`
        set SAT0 = `echo $SAT1`
      endif

    endif
  end 

  # proces the last set of files
  echo "" | awk '{printf("%s ","Combing")}' 
  set jj = 1
  set t2 = 9999999999
  if (-f tmprecord_new) rm -f tmprecord_new

  foreach line2 (`awk '{print $0}' tmprecord`)
    echo $line2 | awk '{printf("%s ",$1)}'
          set tt = `echo $line2 | awk '{print substr($1,length($1)-54,15)}'`
       #   set t1 = `gdate -jf "%Y%m%dT%H%M%S" $tt +%s`
         set ss1 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
          set t1 = `gdate --date="$ss1" +%s`
#          if ($t1 > $t2) set jj = 0
         if ($t2 == "9999999999") then
             echo $line2 >> tmprecord_new
         else
           if ($t1 <= $t2) then
             echo $line2 >> tmprecord_new
           endif
         endif
         
        set tt = `echo $line2 | awk '{print substr($1,length($1)-38,15)}'`
        set ss2 = `echo $tt|awk '{print substr($1,1,4)"/"substr($1,5,2)"/"substr($1,7,2)" "substr($1,10,2)":"substr($1,12,2)":"substr($1,14,2)}'`
#          set t2 = `gdate -jf "%Y%m%dT%H%M%S" $tt +%s`
          set t2 = `gdate --date="$ss2" +%s`
  end
  echo "" | awk '{printf("%s\n","...")}'

  # get the orbit file names and download

  set n1 = ` gdate --date="$date0 - 1 day" +%Y%m%d `
  set n2 = ` gdate --date="$date0 + 1 day" +%Y%m%d `
  if ($SAT0 == "S1A") then
    python $code_path/get_s1_orbits.py $date0 s1a
    #get_s1a_orbit.csh $date0
    set orb = ` find . -name "$SAT0*$n1*$n2*.EOF" ` # the name of the file we just copied
  else
    python $code_path/get_s1_orbits.py $date0 s1b
    #get_s1b_orbit.csh $date0
    set orb = ` find . -name "$SAT0*$n1*$n2*.EOF" ` # the name of the file we just copied
  endif

  set pin1 = `head -1 $2` 
  set f1 = `head -1 tmprecord_new`
  make_s1a_tops $f1/annotation/*iw1*vv*xml $f1/measurement/*iw1*vv*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orb tmp2
  set tmpazi = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
  shift_atime_PRM.csh tmp2.PRM $tmpazi
  set azi1 = `echo $pin1 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`

  set pin2 = `tail -1 $2` 
  set f2 = `tail -1 tmprecord_new`
  make_s1a_tops $f2/annotation/*iw1*vv*xml $f2/measurement/*iw1*vv*tiff tmp2 0
  ext_orb_s1a tmp2.PRM $orb tmp2
  set tmpazi = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5)}'`
  shift_atime_PRM.csh tmp2.PRM $tmpazi
  set azi2 = `echo $pin2 | awk '{print $1,$2,0}' | SAT_llt2rat tmp2.PRM 1 | awk '{printf("%d",$2+0.5 + '$tmpazi')}'`
  set nl = `grep num_lines tmp2.PRM | awk '{print $3}'`

  if ($azi1 >= 0 && $azi2 < $nl) then  
    awk '{print $1","$2","$3}' $2 > tmpllt
    set pin0 = `awk NR==1'{print $0}' tmpllt`
    foreach line2 (`awk '{print $0}' tmpllt`)
      if ($line2 != $pin0) then
        echo $pin0 | awk -F"," '{print $1,$2}' > tmp1llt
        echo $line2 | awk -F"," '{print $1,$2}' >> tmp1llt
        if ($mode != 1) then
          create_frame_tops.csh tmprecord_new $orb tmp1llt 1
         #set newfile = `cat newfile`
         #echo "Created new file " $newfile
          #echo ""
          if (! -d $direname) mkdir $direname
          #if (! -d $direname/data) mkdir $direname/data
          mv *.SAFE $direname
        else
          echo ""
          echo "Frames on date $date0 will be re-organized..."
          echo ""
        endif
        set pin0 = `echo $line2`
      endif
    end
   echo $date0 >>dates.merge.$direname
  else 
    echo $date0 >> dates.skip.$direname
    if ($jj == 0) then
      echo ""
      echo "SKIP $date0, as it stopped observation in the middle ..."
      echo ""
    else
      echo ""
      echo "SKIP $date0, as it does not have enough scenes ..."
      echo ""
    endif
  endif   


  rm tmp*
  #rm *.EOF


