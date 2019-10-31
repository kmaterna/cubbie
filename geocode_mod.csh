#!/bin/csh

alias rm 'rm -f'
unset noclobber
if ( -f ~/.quiet ) then
    set V = ""
else
        set V = "-V"
endif

if ($#argv != 4) then
echo ""
echo "Usage: geocode_mod.csh vel.grd vel_ll.grd vel_ll directory"
echo "  geocode a grid file. Call this from the processing directory, and make sure topo/ has topo/trans.dat"
echo "  vel.grd and vel_ll.grd LIVE IN directory."
echo "  vel_ll is just the name for the kml file."
echo "  Most processing takes place in directory."
echo "  outputs:"
echo "    vel_ll.grd"
exit 1
endif


set ifile = $1
set ofile = $2
set kmlfile = $3
set directory = $4

cd $directory
ln -s ../../topo/trans.dat .


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

# gmt makecpt -T$Zmin/$Zmax/$Zinterval -Cjet -Z > defcolors.cpt
set Zmin = -5
set Zmax = 26
set Zinterval = 1
gmt makecpt -T$Zmin/$Zmax/$Zinterval -Cpolar -Z -Ic > defcolors.cpt

gmt grdimage $ifile -JX6.5i -Cdefcolors.cpt -B"$boundR":Range:/"$boundA":Azimuth:WSen -X1.3i -Y3i -P -K > $kmlfile.ps
gmt psscale -D3.3/-1.5/7/0.3h -Cdefcolors.cpt -B"$Zlabel":"mm/yr": -O >> $kmlfile.ps


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
# IMORTANT: CHANGE YOUR SAMPLING FOR GEOCODED VALUES
set incs="6.0s/4.5s" # choosing pretty high resolution  # THIS IS THE KEY. THIS IS A PARAMETER THAT YOU SHOULD ADJUST. 
# Originally was "2.0s/1.5s"
# If you ask for a higher geocoded resolution than your RADAR data actually has, then you 
# will end up filling your vel_ll with nans. 

set R =  `gmt gmtinfo llp -I$incs -bi3f `
gmt xyz2grd llp $R -I$incs  -r -fg -G$ofile -bi3f

# clean
rm rap* llp raln ralt


echo "Finished with proj_ra2ll.csh"

gmt grdedit -D//"mm/yr"/1///"$PWD:t velocity"/"$remarked" $ofile

cp $kmlfile.ps 2$kmlfile.ps
grd2kml.csh $kmlfile defcolors.cpt

echo "Cleaning up"

rm gmt.history gmt.conf defcolors.cpt
