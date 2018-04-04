#!/bin/csh

alias rm 'rm -f'
unset noclobber
if ( -f ~/.quiet ) then
    set V = ""
else
        set V = "-V"
endif


set ifile = vel.grd
set ofile = vel_ll.grd




set maker = $0:t
set today = `date`
set remarked = `echo by $USER on $today with $maker`
echo remarked is $remarked

set boundR = `gmt grdinfo $ifile -C | awk '{print ($3-$2)/4}'`
set boundA = `gmt grdinfo $ifile -C | awk '{print ($5-$4)/4}'`
gmt makecpt -T-150/150/5 -Cjet -Z > defcolors.cpt
gmt grdimage $ifile -JX6.5i -Cdefcolors.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > vel.ps
gmt psscale -D3.3/-1.5/5/0.2h -Cdefcolors.cpt -B1.57:"mm/yr": -O >> vel.ps


echo "starting proj_ra2ll.csh"  

proj_ra2ll.csh trans.dat $ifile $ofile

echo "Finished with proj_ra2ll.csh"

gmt grdedit -D//"mm/yr"/1///"$PWD:t velocity"/"$remarked" $ofile

grd2kml.csh vel_ll defcolors.cpt
