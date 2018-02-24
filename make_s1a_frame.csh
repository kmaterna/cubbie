#!/bin/csh -f

if ($#argv != 2) then
   echo ""
   echo "Usage: make_s1a_frame.csh data.list frames.ll"
   echo "data.list: /Users/kmaterna/.../S1A_IW_SLC__1SDV_20171201T142315_20171201T142342_019510_0211C3_2A60.SAFE one per line. "
   echo "frames.ll: lon1 lat1 0 \n lon2 lat2 0"
   echo ""
   exit 1
endif

set frm_ll = $2
set nlines = `cat $frm_ll|wc -l`
set nframe = `echo $nlines |awk '{print $1-1}'`
echo "Number of Frames to make: " $nframe

set iframe = 1

while ($iframe <= $nframe)
   echo "Working on Frame_"$iframe
   set i1 = $iframe
   set i2 = `echo $iframe|awk '{print $1+1}'`
   awk 'NR== "'$i1'" {print $0}' $frm_ll > frame$iframe.ll
   awk 'NR== "'$i2'" {print $0}' $frm_ll >>frame$iframe.ll 
   frame_s1a.csh $1 frame$iframe.ll 2 FRAME_$iframe
   @ iframe = $iframe + 1
end
