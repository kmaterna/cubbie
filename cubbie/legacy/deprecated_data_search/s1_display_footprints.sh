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

# I will add some new comments and functions. 

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
# Are we getting a problem with MULTIPOLYGON? 

# Find descriptive numbers for summarizing the search results
# Build the $timing file
num_results=`cat $footprints | wc -l`  # counting the results 
cut -c25-32 $timing > temp_timing.txt  # finding the timing of the acquisitions. 
rm $timing
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "$line 0.5" >> $timing
done < temp_timing.txt
rm temp_timing.txt;

# Build the footprints file. 
am_i_linux=`uname -a | grep 'Linux'`
if [ ! -z "$am_i_linux" ]; then # WE ARE ON A LINUX MACHINE 
    sed -i 's/<str name=\"footprint\">POLYGON //g' $footprints
    sed -i 's/<str name=\"footprint\">MULTIPOLYGON //g' $footprints  # added new
    sed -i 's/))<\/str>//g' $footprints
    sed -i $'s/((/>\\\n/g' $footprints
    sed -i $'s/,/\\\n/g' $footprints
    sed -i 's/(/ /g' $footprints  # addd new
    sed -i 's/)/ /g' $footprints  # added new

    # Making a file with two columns: lon, lat
    cp $footprints new_footprints.txt
    sed -i $'s/>//g' new_footprints.txt
    sed -i '/^$/d' new_footprints.txt

else       # WE ARE ON A MAC or OTHER non-LINUX MACHINE
    sed -i '' 's/<str name=\"footprint\">POLYGON //g' $footprints
    sed -i '' 's/<str name=\"footprint\">MULTIPOLYGON //g' $footprints
    sed -i '' 's/))<\/str>//g' $footprints
    sed -i '' $'s/((/>\\\n/g' $footprints
    sed -i '' $'s/,/\\\n/g' $footprints
    sed -i '' 's/(/ /g' $footprints
    sed -i '' 's/)/ /g' $footprints

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
projection="M4.5i"
range="$lonmin/$lonmax/$latmin/$latmax"
echo $range
echo "Results displayed: " $num_results
gmt pscoast -R$range -J$projection -Dh  -N1/thick,black -N2/thick,black -Bp1.0 -B+t"Displaying $num_results Results" -P -Wthick,black -Gwhite -Swhite -K > $mapfile
# gmt pscoast -R -J -Dh -N1/thick,black -N2/thick,black -Wthick,black -O -K >> $mapfile
gmt psxy $footprints -R -J -Wthick,red -K -O -P >> $mapfile

# Nice little thing: adding the point or rectangle that was used for searching. 
python $DIR/get_search_bounds.py $search_file_name $bounds_file
num_lines=`cat $bounds_file | wc -l` 
if [ $num_lines -eq 1 ]; then
  gmt psxy $bounds_file -R$range -J$projection -Sc0.25 -Gblack -K -O -P >> $mapfile
else
  gmt psxy $bounds_file -R$range -J$projection -Wthick,black -K -O -P >> $mapfile
fi
gmt psxy $footprints -R -J -Wthick,red -O -P >> $mapfile
gmt psconvert -Tg $mapfile


# Make the timing plot
projection="X10iTi/0.9i" #Make an xy projetion 6 inches in the horizontal direction and 2 inches in the vertical direction
region1="2014-10-01T00:00/2018-01-01T00:00/0.1/1"  # The beginning and end of Sentinel.
region2="2018-01-01T00:00/2021-04-01T00:00/0.1/1"  # The beginning and end of Sentinel.
region3="2021-04-01T00:00/2024-07-01T00:00/0.1/1"  # The beginning and end of Sentinel.

gmt psbasemap -R$region1 -J$projection -Bpxa6Of2O -Bpya2 -BWeSn+t"Displaying $num_results Acquisitions" -Bsxa1YS -K -Y14 --FORMAT_DATE_MAP=mm/dd > $timing_file
# -Bpx = primary x-axis; -Bs = secondary. a6O means primary annotate every 6 months; f2O means secondary annotate every 2 months
gmt psxy $timing -R$region1 -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gblack -K -O >> $timing_file  # plotting the actual data. 

gmt psbasemap -R$region2 -J$projection -Bpxa6Of2O -Bpya2 -BWeSn -Bsxa1YS -K -O -Y-5 --FORMAT_DATE_MAP=mm/dd >> $timing_file
gmt psxy $timing -R$region2 -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gblack -O -K >> $timing_file  # plotting the actual data.

gmt psbasemap -R$region3 -J$projection -Bpxa6Of2O -Bpya2 -BWeSn -Bsxa1YS -K -O -Y-5 --FORMAT_DATE_MAP=mm/dd >> $timing_file
gmt psxy $timing -R$region3 -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gblack -O >> $timing_file  # plotting the actual data.
gmt psconvert -Tg $timing_file

rm gmt.history new_footprints.txt
rm $footprints
rm $bounds_file
#open $mapfile
#open $timing_file
