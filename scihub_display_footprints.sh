#!/bin/bash
# SCRIPT TO DISPLAY THE FOOTPRINTS OF SENTINEL-1 SEARCH RESULTS
# Feb. 14, 2018

# Displaying the footprints
# Parse inputs
if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script makes a GMT plot of the footprint of data queries"
  echo "Usage: ./scihub_display_footprints.sh -options"
  echo "Example: ./scihub_display_footprints.sh -i my_out.txt -i my_out2.txt"
  echo "Please provide one or more input files."
  echo ""
  exit 1
fi

# Configure the output files
echo "Displaying Footprints of Search Results"
footprints="footprints.txt"
mapfile="footprint_map.ps"
am_i_linux=`uname -a | grep 'Linux'`  # this tells us if we're on a linux machine (sed works slightly differently on mac and linux)

# Read the search results. It could be multiple calls of the -i flag. 
while getopts i: opt; do
    case $opt in
        i) multi+=("$OPTARG");;
    esac
done
shift $((OPTIND -1))
# echo "The whole list of values is '${multi[@]}'"  # a debugging line. 

rm $footprints
for val in "${multi[@]}"; do
    grep 'name=\"footprint' $val >> $footprints
done

am_i_linux=`uname -a | grep 'Linux'`
if [ ! -z "$am_i_linux" ]; then # WE ARE ON A LINUX MACHINE 
    sed -i 's/<str name=\"footprint\">POLYGON //g' $footprints
    sed -i 's/))<\/str>//g' $footprints
    sed -i $'s/((/>\\\n/g' $footprints
    sed -i $'s/,/\\\n/g' $footprints

    # Making a file with two columns: lon, lat
    cp $footprints new_footprints.txt
    sed -i $'s/>//g' new_footprints.txt
    sed -i '/^$/d' new_footprints.txt
else       # WE ARE ON A MAC or OTHER non-LINUX MACHINE
    sed -i '' 's/<str name=\"footprint\">POLYGON //g' $footprints
    sed -i '' 's/))<\/str>//g' $footprints
    sed -i '' $'s/((/>\\\n/g' $footprints
    sed -i '' $'s/,/\\\n/g' $footprints

    # Making a file with two columns: lon, lat
    cp $footprints new_footprints.txt
    sed -i '' $'s/>//g' new_footprints.txt
    sed -i '' '/^$/d' new_footprints.txt
fi

# find the max and min of the range of footprints. 
lonmin=`awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1};} END {printf "%.2f", min}' new_footprints.txt`
lonmax=`awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1};} END {printf "%.2f", max}' new_footprints.txt`
latmin=`awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2};} END {printf "%.2f", min}' new_footprints.txt`
latmax=`awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2};} END {printf "%.2f", max}' new_footprints.txt`


# Make GMT plot 
projection="M6.1i"
range="$lonmin/$lonmax/$latmin/$latmax"
gmt pscoast -R$range -J$projection -Dh -N2 -B1.0 -P -Wblack -Gwhite -Swhite -K > $mapfile
gmt psxy $footprints -R$range -J$projection -Wthick,red -K -O -P >> $mapfile

rm gmt.history new_footprints.txt
#open $mapfile
