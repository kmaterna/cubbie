#!/bin/bash
# Unzip data as shown in the DATA/ folder.
# Feb. 17, 2021.

if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script unzips the results of data queries in a sub-directory called DATA/"
  echo "one level below where the script is called."
  echo "Usage: ./scihub_unzip_s1_sar.sh -options"
  echo "  Input options:"
  echo " -p polarization (vv, vh, both. vv default)"
  echo "Example: ./scihub_download_s1_sar.sh -p vv"
  echo ""
  exit 1
fi


# Read the run string.
polarization='vv'
while getopts p: opt; do
    case $opt in
      p)  # polarization
        echo "-p was triggered, parameter: $OPTARG" >&2
        polarization=$OPTARG
        ;;  
      \?)
        echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $((OPTIND -1))
# echo "The whole list of values is '${multi[@]}'"  # a debugging line. 


# Defining parameters
am_i_linux=`uname -a | grep 'Linux'`
echo $am_i_linux

cd DATA
ls *.zip | while read f; do
  title=$(echo $f | cut -d'.' -f 1);  # extract the name of the zip/SAFE file.
  echo $title
  if [ $polarization = 'vv' ]
  then
    echo "vv";
    unzip $title.zip
    rm $title.SAFE/measurement/*-slc-vh-*.tiff
    rm $title.zip
  elif [ $polarization = 'vh' ]
  then
    echo "vh";
    unzip $title.zip
    rm $title.SAFE/measurement/*-slc-vv-*.tiff
    rm $title.zip
  elif [ $polarization = 'both' ]
  then
    echo 'both'
    unzip $title.zip
    rm $title.zip
  else
    echo "polarization not recognized"
  fi
done;
cd ../
