#!/bin/bash
# Download data as shown in the search_results. 
# Feb. 14, 2018.

if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script downloads the results of data queries"
  echo "Usage: ./scihub_download_s1_sar.sh -options"
  echo "Example: ./scihub_download_s1_sar.sh -i search_results1.txt -i search_results2.txt"
  echo "Please provide one or more input files."
  echo ""
  exit 1
fi


# Read the search results. It could be multiple calls of the -i flag. 
while getopts i: opt; do
    case $opt in
        i) multi+=("$OPTARG");;
    esac
done
shift $((OPTIND -1))
# echo "The whole list of values is '${multi[@]}'"  # a debugging line. 


# Defining parameters
id_results=uuid_file.txt
am_i_linux=`uname -a | grep 'Linux'`

# Where will the data live? 
mkdir -p DATA
#mkdir -p MANIFEST


# WILL FIX THIS TO CAT ANY DATA IN ANY FILE
# THIS WILL BE IN A LOOP OVER POTENTIALLY MULTIPLE $RAW_RESULTS files
# Processing the raw results to get unique id names
rm $id_results
for val in "${multi[@]}"; do
    grep -E 'uuid|<title>S1' $val >> $id_results
done


# the -i '' is because of mac computers. Might need to delete the '' on a linux machine. 
if [ ! -z am_i_linux ]; then  # I am on a linux machine
    sed -i 's/<str name=\"uuid\">//g' $id_results
    sed -i 's/<title>//g' $id_results
    sed -i 's/<\/title>//g' $id_results
    sed -i 's/<\/str>//g' $id_results
else  # we are on a mac or non-linux machine! Macs need the -i '' in the sed call.
    sed -i '' 's/<str name=\"uuid\">//g' $id_results
    sed -i '' 's/<title>//g' $id_results
    sed -i '' 's/<\/title>//g' $id_results
    sed -i '' 's/<\/str>//g' $id_results
fi

counter=0
while read p; do
  if [ $counter = 0 ]; then
  	title=$p
    
    # In this version, I actually download from the ASF. The download goes about 5x faster for some reason. 
    wget -c --http-user=kmaterna@berkeley.edu --http-password=Access_d4t4 -O DATA/"$title".zip "https://datapool.asf.alaska.edu/SLC/SA/$title.zip"
  	counter=1
  	continue
  else
    uuid=$p
  fi
  echo $title
  #echo $uuid
  
  # I have the option to download the Manifest or Data from ESA. It's pretty slow though. 
  # MANIFEST only
  # wget --no-check-certificate --user=kmaterna --password=access_data -O MANIFEST/"$title"_manifest.safe "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/Nodes('$title.SAFE')/Nodes('manifest.safe')/\$value"

  # DATA FROM COPERNICUS (full thing- will take a long time)!
  # wget --no-check-certificate --user=kmaterna --password=access_data -O DATA/"$title".SAFE.zip "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value"
  # Takes a few hours for each SAFE.zip. 
  # Each one can be unzipped with unzip. 
  
  cd DATA
  unzip $title.zip
  rm $title.SAFE/measurement/*-slc-vh-*.tiff
  rm $title.zip
  cd ../
  counter=0
  
done <$id_results
