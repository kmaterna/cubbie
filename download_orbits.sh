#!/bin/bash
################################################
#      Sentinel-1 Precise Orbit downloader
#     Joaquin Escayo 2016 j.escayo@csic.es
###############################################
# Version 2.0
# Requisites: bash, wget, sed, sort (preinstalled)
# Recommended use of cron to schedule the execution
# Version History:
# v 1.0 - Initial release
# v 2.0 - Now it also download EAP Phase calibration files. Bugfixes. (18/10/2016)
# v 2.1 - Minor revision, now download of AUX_CAL data is optional, directory must be set to download. (19/10/2016)
# TO-DO:
# 1. Detection of corrupted files (incompleted downloads)
# Modified by K. Materna to work on mac and linux, and to use start/end pages. 

###############################
#  CONFIGURATION PARAMETERS   #
###############################
# Download directory:
# Uncomment this line and set the correct directory
DOWN_DIR="S1_orbits" # directory for store precise orbits files
CAL_DIR="S1_orbits" # directory for store AUX_CAL files
# number of pages to check
PAGE_START=1; # website pages for precise orbits
PAGE_END=32 # Pages for precise orbits
CAL_PAGES=4 # Calibration pages

###############################
#        TEMPORAL FILES       #
###############################
# lista de archivos a descargar
list=$(mktemp /tmp/s1list.XXXXX)
# índice html
index=$(mktemp /tmp/s1index.XXXXX)
# lista de archivos en el servidor
remote_files=$(mktemp /tmp/s1remote.XXXXX)
# lista de archivos a descargar
dw_list=$(mktemp /tmp/s1dw.XXXXX)
# lista de archivos previamente descargados
local_files=$(mktemp /tmp/s1no.XXXXX)
am_i_linux=`uname -a | grep 'Linux'`  # this tells us if we're on a linux machine

###############################
#      SCRIPT EXECUTION       #
###############################

# check if download directory is set
if [ -z $DOWN_DIR ]; then
    echo "#######################################################"
    echo "You must set download directory before use this program"
    echo "Edit sentinel1_orbit_downloader.sh and set the value of"
    echo "DOWN_DIR and (OPTIONAL) CAL_DIR variables"
    echo "#######################################################"
    # cleanup
    rm $list
    rm $index
    rm $remote_files
    rm $dw_list
    rm $local_files
    exit
fi

# check if DOWN_DIR exists
mkdir -p $DOWN_DIR

# check if CAL_DIR exists
if [ ! -d "$CAL_DIR" ] && [ ! -z "$CAL_DIR" ]; then
  mkdir $CAL_DIR
fi

# Orbits download

for i in $(eval echo "{$PAGE_START..$PAGE_END}")
do
	wget --quiet -O - --no-check-certificate https://qc.sentinel1.eo.esa.int/aux_poeorb/?page=$i >> $index
done

# Generando la lista de ficheros a descargar
# Me quedo sólo con el nombre de los ficheros
grep -o '<a .*href=.*>' $index | grep 'EOF' | sed -e 's/<a .*href=['"'"'"]//' -e 's/["'"'"'].*$//' -e '/^$/ d' | sort -u > $remote_files

echo "Remote files are: "
cat $remote_files

# compruebo que se haya generado correctamente la lista de archivos a descargar:
if ! [ -s $remote_files ]
then
    echo "ERROR OCURRED, NO REMOTE FILES FOUND"
    # cleanup
    rm $list
    rm $index
    rm $remote_files
    rm $dw_list
    rm $local_files
    exit
fi

# Elimino de la lista los archivos ya descargados

if [ "$(ls $DOWN_DIR)" ]; then
    ls $DOWN_DIR > $local_files
    awk 'NR==FNR{a[$0]=1;next}!a[$0]' $local_files $remote_files > $dw_list
    # cp $remote_files $dw_list
else
    cp $remote_files $dw_list
fi

# #Generando los enlaces
# if [ ! -z "$am_i_linux" ]; then # WE ARE ON A LINUX MACHINE 
#     sed -i 's/^/https:\/\/qc.sentinel1.eo.esa.int\/aux_poeorb\//' $dw_list
# else
#     sed -i '' 's/^/https:\/\/qc.sentinel1.eo.esa.int\/aux_poeorb\//' $dw_list
# fi
# Commented july 2018 to resolve hosts with proper URL. 

echo "dw_list:"
cat $dw_list

wget --no-check-certificate -P $DOWN_DIR -i $dw_list

# Cleaning of the temporal files
cat /dev/null > $index
cat /dev/null > $list
cat /dev/null > $remote_files
cat /dev/null > $dw_list
cat /dev/null > $local_files

# ------------------ #
# Calibration files  #
# ------------------ #
if [ ! -z $CAL_DIR ]; then
    for i in $(eval echo "{1..$CAL_PAGES}")
    do
	    wget --quiet -O - --no-check-certificate https://qc.sentinel1.eo.esa.int/aux_cal/?page=$i >> $index
    done

    # Generando la lista de ficheros a descargar
    # Me quedo sólo con el nombre de los ficheros
    grep -o '<a .*href=.*>' $index | grep 'SAFE' | sed -e 's/<a .*href=['"'"'"]//' -e 's/["'"'"'].*$//' -e '/^$/ d' | sort -u > $remote_files

    # compruebo que se haya generado correctamente la lista de archivos a descargar:
    if ! [ -s $remote_files ]
    then
        echo "ERROR OCURRED, NO REMOTE FILES FOUND"
        # cleanup
        rm $list
        rm $index
        rm $remote_files
        rm $dw_list
        rm $local_files
        exit
    fi

    # Elimino de la lista los archivos ya descargados

    if [ "$(ls $CAL_DIR)" ]; then
        ls $CAL_DIR > $local_files
        awk 'NR==FNR{a[$0]=1;next}!a[$0]' $local_files $remote_files > $dw_list
    else
        cp $remote_files $dw_list
    fi

    #Generando los enlaces
    if [ ! -z "$am_i_linux" ]; then # WE ARE ON A LINUX MACHINE 
        sed -i 's/^/https:\/\/qc.sentinel1.eo.esa.int\/aux_poeorb\//' $dw_list
    else
        sed -i '' 's/^/https:\/\/qc.sentinel1.eo.esa.int\/aux_poeorb\//' $dw_list
    fi

    wget --quiet --no-check-certificate -P $CAL_DIR -i $dw_list
fi

###############################
#          CLEANUP            #
###############################
rm $list
rm $index
rm $remote_files
rm $dw_list
rm $local_files
