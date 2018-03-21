#!/bin/bash
# Go into a group of intf_all folders and copy some file of interest into a new folder

# ls . > filelist.txt
# mkdir phasefilt
# while IFS= read -r var
# do 
#   echo $var
#   newname=$var"_phasefilt_mask.grd"
#   echo $newname
#   cp $var/phasefilt_mask.grd $newname
# done < filelist.txt
# mv *.grd phasefilt



# Unwrapped radar coordinates
ls . > filelist.txt
mkdir unwrap_ra
while IFS= read -r var
do
	echo $var
	newname=$var"_unwrap.grd"
	echo $newname
	cp $var/unwrap.grd $newname
done < filelist.txt
mv *.grd unwrap_ra