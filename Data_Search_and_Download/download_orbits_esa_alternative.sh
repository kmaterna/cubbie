#!/bin/bash
################################################
#      Sentinel-1 Precise Orbit downloader
#      alternative needed after
#      qc.sentinel1.eo.esa.int stopped responding,
#      March 2021
#      source: http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/
#      don't use; not as reliable as ASF.
###############################################

###############################
#  CONFIGURATION PARAMETERS   #
###############################
# Download directory:
# Uncomment this line and set the correct directory
DOWN_DIR="S1_orbits_2021" # directory for store precise orbits files

###############################
#      SCRIPT EXECUTION       #
###############################

# check if download directory is set
if [ -z $DOWN_DIR ]; then
    echo "#######################################################"
    echo "You must set download directory before use this program"
    echo "Edit sentinel1_orbit_downloader.sh and set the value of"
    echo "DOWN_DIR variable"
    echo "#######################################################"
    # cleanup
    exit
fi

# check if DOWN_DIR exists
mkdir -p $DOWN_DIR
cd $DOWN_DIR

## # Download the remote directory
#wget \
#    --accept '*.zip' \
#    --execute robots=off \
#    --recursive \
#    --level=0 \
#    --no-parent \
#    --spider \
#    'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/S1A/2021' 2>&1 | tee -a main.log
#wget \
#    --accept '*.zip' \
#    --execute robots=off \
#    --recursive \
#    --level=0 \
#    --no-parent \
#    --spider \
#    'http://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/S1B/' 2>&1 | tee -a main.log

# # Get the list of files
#grep '^--' main.log > list_of_addresses.txt
#grep '.zip' list_of_addresses.txt > zip_files.txt

# # Wget all the files
while read p; do
  stringarray=($p)
  full_address=${stringarray[2]}   # the http address of the orbit file
  echo $full_address
  name_only=`echo "$full_address" | awk -F/ '{print $NF}'`   # the filename only
  if [ ! -f "$name_only" ]; then
    wget $full_address
  fi
done <zip_files.txt

# Now unzip all the files
# Right now I'm unzipping them manually with right-click and then deleting the tmp folder that is created.
# Unzipping does not take long at all.

## Sometimes the unzipping doesn't automatically put the EOF file into the top directory.
## So we have to do it manually.
#find . -type d -name '*.EOF' > directory_list.txt
#
#while read p; do
#  mv $p $p"_tmp"
#done <directory_list.txt
#
#cp *.EOF_tmp/tmp/*.EOF .
#rm -r *.EOF_tmp/
#rm directory_list.txt