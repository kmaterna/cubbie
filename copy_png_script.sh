#!/bin/bash
# Go into a group of intf_all folders and copy some file of interest into a new folder

ls . > filelist.txt
mkdir phasefilt
while IFS= read -r var
do 
  echo $var
  newname=$var"_phasefilt_mask.ps"
  echo $newname
  cp $var/phasefilt_mask.ps $newname
done < filelist.txt
mv *.ps phasefilt
