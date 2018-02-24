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


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Configure the output files
echo "Displaying Footprints of Search Results"
footprints="footprints.txt"
timing="timing.txt"
mapfile="footprint_map.ps"
timing_file="acquisitions.ps"
bounds_file="bounds_file.txt"
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
rm $timing
for val in "${multi[@]}"; do
    search_file_name=$val
    grep 'name=\"footprint' $val >> $footprints
    grep 'title>S1' $val >> $timing
done
 
# Find descriptive numbers for summarizing the search results
num_results=`cat $footprints | wc -l`  # counting the results 
cut -c25-32 $timing > temp_timing.txt  # finding the timing of the acquisitions. 
rm $timing
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "$line 0.5" >> $timing
done < temp_timing.txt
rm temp_timing.txt;

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


# Make GMT plot of footprints 
# In the best case scenario I would also plot a timeseries of the acquisition dates (nice and fancy)
projection="M6.1i"
range="$lonmin/$lonmax/$latmin/$latmax"
echo "Results displayed: " $num_results
gmt pscoast -R$range -J$projection -Dh -N2 -Bp1.0 -B+t"Displaying $num_results Results" -P -Wblack -Gwhite -Swhite -K > $mapfile
gmt psxy $footprints -R$range -J$projection -Wthick,red -K -O -P >> $mapfile

# Nice little thing: adding the point or rectangle that was used for searching. 
python $DIR/get_search_bounds.py $search_file_name $bounds_file
num_lines=`cat $bounds_file | wc -l` 
if [ $num_lines -eq 1 ]; then
  gmt psxy $bounds_file -R$range -J$projection -Sc0.25 -Gblack -K -O -P >> $mapfile
else
  gmt psxy $bounds_file -R$range -J$projection -Wthick,black -K -O -P >> $mapfile
fi


# Make the timing plot
projection="X10iTi/4i" #Make an xy projetion 6 inches in the horizontal direction and 2 inches in the vertical direction
region="2014-10-01T00:00/2018-05-30T00:00/0.1/1"  # The beginning and end of Sentinel. 

gmt psbasemap -R$region -J$projection -Bpxa6Of2O -Bpya2 -BWeSn+t"Displaying $num_results Acquisitions" -Bsxa1YS -K --FORMAT_DATE_MAP=mm/dd > $timing_file
# -Bpx = primary x-axis; -Bs = secondary. a6O means primary annotate every 6 months; f2O means secondary annotate every 2 months
gmt psxy $timing -R$region -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gblack -K -O >> $timing_file  # plotting the actual data. 

rm gmt.history new_footprints.txt
rm $footprints
rm $bounds_file
#open $mapfile
#open $timing_file
