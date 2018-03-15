#!/bin/bash

output="ps.ps"
range="-124.5/-122.5/39.5/41.5"
projection="M6.1i"


gmt pscoast -R$range -J$projection -Gwhite -Sgray -Dh -B1.0 -Wthin,black -K > $output

gmt grdimage vel_ll.grd -Cvel_ll.cpt -R$range -J$projection -K -O >> $output

gmt psscale -D0.2i/1.1i/-4c/0.6c -Cvel_ll.cpt -B50.0:"LOS":/:mm/yr: -K -O >> $output