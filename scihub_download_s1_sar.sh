#!/bin/bash
# Download data as shown in the search_results. 
# Feb. 14, 2018.

# Defining parameters
raw_results=search_results.txt
id_results=uuid_file.txt
am_i_linux=`uname -a | grep 'Linux'`

# Where will the data live? 
mkdir -p DATA
#mkdir -p MANIFEST

# Processing the raw results to get unique id names
grep -E 'uuid|<title>S1' $raw_results > $id_results

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
  	counter=1
  	continue
  else
    uuid=$p
  fi
  echo $title
  echo $uuid
  
  # MANIFEST only
  # wget --no-check-certificate --user=kmaterna --password=access_data -O MANIFEST/"$title"_manifest.safe "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/Nodes('$title.SAFE')/Nodes('manifest.safe')/\$value"

  # DATA (full thing- will take a long time)!
  wget --no-check-certificate --user=kmaterna --password=access_data -O DATA/"$title".SAFE.zip "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value"
  # Takes a few hours for each SAFE.zip. 
  # Each one can be unzipped with unzip. 

  unzip $title.zip
  rm $title.zip
  counter=0
  
done <$id_results
