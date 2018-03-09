#!/bin/bash

# Get the average correlation coefficient for each interferogram

ls intf_all > intf_data.txt
results="corr_results.txt"
rm $results


while IFS='' read line; do

	cd intf_all/$line
	gmt grdmath corr_ll.grd MEAN = out.grd
	correlation=`gmt grdinfo out.grd | grep z | awk '{print $3}'`
	rm out.grd
	cd ../../

	echo $line $correlation >> $results

done < intf_data.txt