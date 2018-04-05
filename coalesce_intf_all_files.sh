#!/bin/bash
# Go into a group of intf_all folders and copy some file of interest into a new folder
# The new folder will live in intf_all as a separate directory with the name you give it ($1)

if [ "$#" -ne 1 ]; then
	echo "ERROR: Please provide a pattern of filenames you'd like to copy. "
	echo "Example: coalesce_intf_all_files.sh unwrap.grd"
	exit 1;
fi


collect_files=$1
destination=intf_all/$1
mkdir -p $destination
ls -d intf_all/???????_??????? > filelist.txt  # list all the directories we go inside. 
while IFS= read -r var
do
	datestring=`echo $var | cut -c10-26`   # selects the '2015177_2015204' part of the string
	echo $datestring
	newname=$datestring"_"$1
	cp $var/$1 $destination/$newname
done < filelist.txt
rm filelist.txt