#!/bin/bash

# Compare the results of FRAMES with the initial acquisitions

# Get initial acquisitions
# Find descriptive numbers for summarizing the search results
acq_timing="acq_timing.txt"
num_results=`ls DATA | wc -l`  # counting the results 
ls DATA > $acq_timing
cut -c18-25 $acq_timing > temp_timing.txt  # finding the timing of the acquisitions. 
rm $acq_timing
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "$line 0.5" >> $acq_timing
done < temp_timing.txt
rm temp_timing.txt;


frame_timing="frame_timing.txt"
num_results=`ls FRAMES/FRAME_1 | wc -l`  # counting the results 
ls DATA > $frame_timing
cut -c18-25 $frame_timing > temp_timing.txt  # finding the timing of the acquisitions. 
rm $frame_timing
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "$line 0.6" >> $frame_timing
done < temp_timing.txt
rm temp_timing.txt;




# Make the timing plot
timing_file="frame_vs_acq.ps"
projection="X10iTi/4i" #Make an xy projetion 6 inches in the horizontal direction and 2 inches in the vertical direction
region="2014-10-01T00:00/2018-05-30T00:00/0.1/1"  # The beginning and end of Sentinel. 

gmt psbasemap -R$region -J$projection -Bpxa6Of2O -Bpya2 -BWeSn+t"Displaying $num_results Acquisitions" -Bsxa1YS -K --FORMAT_DATE_MAP=mm/dd > $timing_file
# -Bpx = primary x-axis; -Bs = secondary. a6O means primary annotate every 6 months; f2O means secondary annotate every 2 months
gmt psxy $acq_timing -R$region -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gblack -K -O >> $timing_file  # plotting the actual data. 
gmt psxy $frame_timing -R$region -J$projection --FORMAT_DATE_IN=yyyymmdd -Sc0.15 -Gred -K -O >> $timing_file  # plotting the actual data. 

rm gmt.history