#!/bin/csh -f 
# Special code by Xiaohua to merge several topo_ra.grd swaths
# Works for 2 or 3 swaths
# Can sometimes mess up the edges by one or two rows


cp ../../F1/topo/master.PRM ./t1.PRM
cp ../../F2/topo/master.PRM ./t2.PRM
# cp ../../F3/topo/master.PRM ./t3.PRM

set r1 = `grep earth_radius t1.PRM | awk '{print $3}'`
set r2 = `grep earth_radius t2.PRM | awk '{print $3}'`
# set r3 = `grep earth_radius t3.PRM | awk '{print $3}'`

echo "computing difference between reference earth radius"
set d1 = `echo $r1 $r2 | awk '{printf("%.6f", $1-$2)}'`
# set d2 = `echo $r2 $r3 | awk '{printf("%.6f", $1-$2)}'`

echo "copying original topo_ra.grd-s and make adjustment"
gmt grdmath ../../F1/topo/topo_ra.grd FLIPUD = t1.grd
gmt grdmath ../../F2/topo/topo_ra.grd FLIPUD = t2.grd
# gmt grdmath ../../F3/topo/topo_ra.grd FLIPUD = t3.grd

echo "t1.PRM:t1.grd" > topolist
echo "t2.PRM:t2.grd" >> topolist
# echo "t3.PRM:t3.grd" >> topolist

merge_swath topolist topo_ra.grd

