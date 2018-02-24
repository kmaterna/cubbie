#!/bin/tcsh -f

if ($#argv != 2) then
   echo ""
   echo "Usage: make_s1a_frame.csh data.list frames.ll"
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
