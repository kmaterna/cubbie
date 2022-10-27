#!/bin/bash
################################################
#      Sentinel-1 Precise Orbit downloader
#      alternative needed after
#      qc.sentinel1.eo.esa.int stopped responding,
#      March 2021
#      source: https://s1qc.asf.alaska.edu/aux_poeorb/
#      PREFERRED METHOD FOR DOWNLOADING ORBITS.  March 2021.
###############################################

###############################
#  CONFIGURATION PARAMETERS   #
###############################
# Download directory:
# Uncomment this line and set the correct directory
DOWN_DIR="S1_orbits" # directory for store precise orbits files. Looks like this is hard-coded right now.

###############################
#      SCRIPT EXECUTION       #
###############################

# check if download directory is set
if [ -z $DOWN_DIR ]; then
    echo "#######################################################"
    echo "You must set download directory before use this program"
    echo "Edit download_orbits_asf.sh and set the value of"
    echo "DOWN_DIR variable"
    echo "#######################################################"
    # cleanup
    exit
fi

# check if DOWN_DIR exists
mkdir -p $DOWN_DIR
cd $DOWN_DIR

# Download the list of available files (essentially an index.html)
wget https://s1qc.asf.alaska.edu/aux_poeorb

# # Get the list of files
grep -E 'S1[A-B]_OPER_AUX_POEORB_OPOD_[0-9]{8}T[0-9]{6}_V[0-9]{8}T[0-9]{6}_[0-9]{8}T[0-9]{6}.EOF' -o aux_poeorb > list_of_addresses.txt

# # Wget all the files
while read p; do
  name_only=$p
  echo "$name_only"
  if [ ! -f "$name_only" ]; then
    wget https://s1qc.asf.alaska.edu/aux_poeorb/"$name_only"
  fi
done <list_of_addresses.txt

rm aux_poeorb
