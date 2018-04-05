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

# Datestamp information
set maker = $0:t
set today = `date`
set remarked = `echo by $USER on $today with $maker`
echo remarked is $remarked

set boundR = `gmt grdinfo $ifile -C | awk '{print ($3-$2)/4}'`
set boundA = `gmt grdinfo $ifile -C | awk '{print ($5-$4)/4}'`
set Zmin = `gmt grdinfo $ifile -M -C | awk '{print $6}'`
set Zmax = `gmt grdinfo $ifile -M -C | awk '{print $7}'`
set Zinterval = `gmt grdinfo $ifile -M -C | awk '{print ($7-$6)/100}'`
set Zlabel = `gmt grdinfo $ifile -M -C | awk '{print int(($7-$6)/6) }'`   # provide labels on the psscale
gmt makecpt -T$Zmin/$Zmax/$Zinterval -Cjet -Z > defcolors.cpt
gmt grdimage $ifile -JX6.5i -Cdefcolors.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > vel.ps
gmt psscale -D3.3/-1.5/7/0.3h -Cdefcolors.cpt -B"$Zlabel":"mm/yr": -O >> vel.ps


echo "starting proj_ra2ll.csh"  

gmt grd2xyz $ifile -s -bo3f > rap
#
#   make grids of longitude and latitude versus range and azimuth unless they already exist
#
if (! -f raln.grd || ! -f ralt.grd ) then
  gmt gmtconvert trans.dat -o0,1,3 -bi5d -bo3f > raln
  gmt gmtconvert trans.dat -o0,1,4 -bi5d -bo3f > ralt
#
gmt surface raln `gmt gmtinfo rap -I16/32 -bi3f` -bi3f -I16/32 -T.50 -Graln.grd $V
gmt surface ralt `gmt gmtinfo rap -I16/32 -bi3f` -bi3f -I16/32 -T.50 -Gralt.grd $V
endif
#
gmt grdtrack rap -nl -Graln.grd -bi3f -bo4f > rapln
gmt grdtrack rapln -nl -Gralt.grd -bi4f -bo5f > raplnlt
#
# get the lon, lat, phase columns and grid
#
gmt gmtconvert raplnlt -bi5f -bo3f -o3,4,2 > llp
#
#
set incs="2.0s/1.5s" # choosing pretty high resolution
set R =  `gmt gmtinfo llp -I$incs -bi3f `
gmt xyz2grd llpb $R -I$incs  -r -fg -G$ofile -bi3f
#
# clean
#
rm rap* llp raln ralt


echo "Finished with proj_ra2ll.csh"

gmt grdedit -D//"mm/yr"/1///"$PWD:t velocity"/"$remarked" $ofile

grd2kml.csh vel_ll defcolors.cpt
