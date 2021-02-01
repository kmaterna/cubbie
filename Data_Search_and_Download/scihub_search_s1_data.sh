#!/bin/bash
# SCRIPT TO QUERY THE COPERNICUS DATABASE FOR INSAR SCENES
# K. Materna, Feb. 2018

# Parse inputs
if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script queries the copernicus scihub for Sentinel-1 SLC TOPS-mode images"
  echo "Usage: ./scihub_search_s1_data.sh -options"
  echo "Example: ./scihub_search_s1_data.sh -s 2015-08-01 -e NOW -r -123.0/-123.3/40.0/40.2 -d Descending"
  echo "  Input options:"
  echo " -s start_time [yyyy-mm-dd or NOW]"
  echo " -e end_time [yyyy-mm-dd or NOW]"
  echo " -b bounding_box [lonW/lonE/latS/latN]"
  echo " -r region (4 sided) [lon1/lat1/lon2/lat2/lon3/lat3/lon4/lat4]"
  echo " -c coordinate [lon/lat]"
  echo " -o orbit_number [0-175]"
  echo " -d direction [Ascending/Descending]"
  echo " -z output_file"
  echo " -u username"
  echo " -p password"
  echo ""
  echo "  outputs:"
  echo "    search_results.txt or specified output_file (contains xml-style information on 0-100 results)"
  echo ""
  echo " Note: To use Scihub, you will generally need a Scihub Copernicus login."
  echo "    You can type them directly into the callstring with -u and -p, or "
  echo "    You can put your credentials into a file called ~/.wget_cred with format:"
  echo "    scihub_username: your_username"
  echo "    scihub_password: your_password"
  echo ""
  echo " Note: The scihub system will only print a maximum of 100 results to a file."
  echo "    If you want to see more than 100 results, this script will cat them, but only up to five hundred right now."
  echo "    More details at https://scihub.copernicus.eu/twiki/do/view/SciHubUserGuide/5APIsAndBatchScripting#URI_components"
  echo ""
  exit 1
fi

# Initial values
output_file=search_results.txt

# Parse arguments (:after options means it expects an argument)
while getopts :s:e:r:c:o:d:b:z:u:p: opt; do
  case $opt in
  	s)  # the start time
      echo "-s was triggered, Parameter: $OPTARG" >&2
      starttime=$OPTARG
      ;;
    e)  # the end time
      echo "-e was triggered, parameter: $OPTARG" >&2
      endtime=$OPTARG
      ;;
    b)  # the box
      echo "-b was triggered, parameter: $OPTARG" >&2
      box=$OPTARG
      components=$(echo $box | tr "/" "\n")  # the '/' symbol is the delimiter
      temp=( $components )
      lonW=${temp[0]}
      lonE=${temp[1]}
      latS=${temp[2]}
      latN=${temp[3]}
      if (( $(echo "$lonW >= $lonE" | bc -l) )); then
      	echo "Error! Longitudes not in increasing order!"
      	exit 1
      fi
      if (( $(echo "$latS >= $latN" | bc -l) )); then
      	echo "Error! Latitudes not in increasing order!"
      	exit 1
      fi
      if (( $(echo "$latN >= 90.0" | bc -l) )); then
      	echo "Error! latN not in bounds!"
      	exit 1
      fi      
      if (( $(echo "$latS >= 90.0" | bc -l) )); then
      	echo "Error! latS not in bounds!"
      	exit 1
      fi
      ;;
    r)  # the region (4 arbitrary points)
      echo "-r was triggered, parameter: $OPTARG" >&2
      region=$OPTARG
      components=$(echo $region | tr "/" "\n")  # the '/' symbol is the delimiter
      temp=( $components )
      lon1=${temp[0]}
      lat1=${temp[1]}
      lon2=${temp[2]}
      lat2=${temp[3]}
      lon3=${temp[4]}
      lat3=${temp[5]}
      lon4=${temp[6]}
      lat4=${temp[7]}
      if (( $(echo "$lat1 >= 90.0" | bc -l) )); then
        echo "Error! lat1 not in bounds!"
        exit 1
      fi      
      if (( $(echo "$lat1 >= 90.0" | bc -l) )); then
        echo "Error! lat1 not in bounds!"
        exit 1
      fi
      ;;
    c)  # the coordinate 
      echo "-p was triggered, parameter: $OPTARG" >&2
      point=$OPTARG
      components=$(echo $point | tr "/" "\n")  # the '/' symbol is the delimiter
      temp=( $components )
      lon=${temp[0]}
      lat=${temp[1]}
      ;;
    o)  # the relative orbit
      echo "-o was triggered, parameter: $OPTARG" >&2
      orbit=$OPTARG
      if [ $orbit -gt 175 ]; then
      	echo "Error: orbit must be between 0 and 175"
      	exit 1
      fi
      if [ $orbit -lt 0 ]; then
      	echo "Error: orbit must be between 0 and 175"
      	exit 1
      fi
      ;;
    d)  # the flight direction
      echo "-d was triggered, parameter: $OPTARG" >&2
      direction=$OPTARG
      if [ "$direction" != "Descending" -a "$direction" != "Ascending" ]; then
      	echo "Error: direction must be 'Ascending' or 'Descending' "
      	exit 1
      fi
      ;;      
    z)  # the name of the output file
      echo "-z was triggered, parameter: $OPTARG" >&2
      output_file=$OPTARG
      ;; 
    u)  # username
      echo "-u was triggered, parameter: $OPTARG" >&2
      username=$OPTARG
      ;;
    p)  # password
      echo "-p was triggered, parameter: $OPTARG" >&2
      password=$OPTARG
      ;;        
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Extracting scihub username and password from ~/.wget_cred in case they weren't provided
if [ -z "$username" ]; then
  echo "No Scihub username provided in callstring; finding your Scihub username in ~/.wget_cred"
  username=`grep 'scihub_username' ~/.wget_cred | awk {'print $2'}`
  echo "Scihub username found: "$username
fi
if [ -z "$password" ]; then
  echo "No Scihub password provided in callstring; finding your Scihub password in ~/.wget_cred"
  password=`grep 'scihub_password' ~/.wget_cred | awk {'print $2'}`
  echo "Scihub password found."
fi

# Get the username and password from the user
if [ -z "$username" ]; then
  echo "No Scihub username has been found. Please read the documentation of this script for details. Exiting..."
  exit 1
fi
if [ -z "$password" ]; then
  echo "No Scihub password has been found. Please read the documentation of this script for details. Exiting..."
  exit 1
fi



search_query="https://scihub.copernicus.eu/dhus/search?q=platformname:Sentinel-1 AND producttype:SLC"
# The basic stuff that doesn't change^^

# CONSTRUCTING THE SEARCH QUERY
# If we are searching based on starttime and endtime
if [ ! -z "$starttime" ]; then
	if [ -z "$endtime" ]; then
		search_query+=" AND beginposition:[$starttimeT00:00:00.000Z TO NOW]"
	else
		search_query+=" AND beginposition:[$starttime"
		search_query+="T00:00:00.000Z TO $endtime"
		search_query+="T00:00:00.000Z]"
	fi
fi

# Searching based on a bounding box LonW/LonE/LatS/LatN
if [ ! -z "$box" ]; then
	search_query+=" AND footprint:\"intersects(POLYGON(("
	search_query+="$lonW $latN,"
	search_query+="$lonE $latN,"
	search_query+="$lonE $latS,"
	search_query+="$lonW $latS,"
	search_query+="$lonW $latN)))\""
fi

# Searching based on a four-sided polygon Lon1/Lat1/Lon2/Lat2/Lon3/Lat3/Lon4/Lat4
if [ ! -z "$region" ]; then
  search_query+=" AND footprint:\"intersects(POLYGON(("
  search_query+="$lon1 $lat1,"
  search_query+="$lon2 $lat2,"
  search_query+="$lon3 $lat3,"
  search_query+="$lon4 $lat4,"
  search_query+="$lon1 $lat1)))\""
fi

# Searching based on a point
if [ ! -z "$point" ]; then
	search_query+=" AND footprint:\"intersects($lat,$lon)\""
fi

# Searching based on a relative orbit
if [ ! -z "$orbit" ]; then
	search_query+=" AND relativeorbitnumber:$orbit"
fi

# Searching based on a flight direction
if [ ! -z "$direction" ]; then
	search_query+=" AND orbitdirection:$direction"
fi

# TOPS Mode instead of strip-map
search_query+=" AND sensoroperationalmode:IW"


# how many rows to display and where to start? 
# Max rows = 100 (slightly annoying rule from the Copernicus server)
search_query0=$search_query"&start=0&rows=100"
echo $search_query0

echo "Input options:" $@ > $output_file
echo "wget --no-check-certificate --user="$username" --password= "$search_query0 >> $output_file
# Execute the search using wget
wget --no-check-certificate --auth-no-challenge --user=$username --password=$password "$search_query0" -O ->> $output_file
num_results=`grep 'title>S1' $output_file | wc -l`

# Execute again if we think there's more search results to be found (100-200 and 200-300 range). 
if [ $num_results -eq "100" ]; then
  echo "We have 100+ results... automatically searching for results #100-200"
  search_query1=$search_query"&start=100&rows=100"
  wget --no-check-certificate --auth-no-challenge --user=$username --password=$password "$search_query1" -O ->> $output_file
fi
num_results=`grep 'title>S1' $output_file | wc -l`
if [ $num_results -eq "200" ]; then
  echo "We have 200+ results... automatically searching for results #200-300"
  search_query2=$search_query"&start=200&rows=100"
  wget --no-check-certificate --auth-no-challenge --user=$username --password=$password "$search_query2" -O ->> $output_file
fi
num_results=`grep 'title>S1' $output_file | wc -l`
if [ $num_results -eq "300" ]; then
  echo "We have 300+ results... automatically searching for results #300-400"
  search_query3=$search_query"&start=300&rows=100"
  wget --no-check-certificate --auth-no-challenge --user=$username --password=$password "$search_query3" -O ->> $output_file
fi
num_results=`grep 'title>S1' $output_file | wc -l`
if [ $num_results -eq "400" ]; then
  echo "We have 400+ results... automatically searching for results #400-500"
  search_query4=$search_query"&start=400&rows=100"
  wget --no-check-certificate --auth-no-challenge --user=$username --password=$password "$search_query4" -O ->> $output_file
fi



# Displaying a summary of the results 
grep 'title>S1' $output_file   # displaying the results
echo "number of total results is:"
grep 'total results' $output_file

echo "number of displayed results is: "
grep 'title>S1' $output_file | wc -l  # counting the results 




