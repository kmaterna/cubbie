#!/bin/bash
# Download data as shown in the search_results. 
# Feb. 14, 2018.

if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script downloads the results of data queries and puts it into a directory called DATA/"
  echo "one level below where the script is called."
  echo "Usage: ./scihub_download_s1_sar.sh -options"  
  echo "  Input options:"
  echo " -i input_file (can have multiple calls to -i)"
  echo " -u username"
  echo " -p password"
  echo "Note: To use ASF, you will need a NASA Earthdata login."
  echo "    You can type them directly into the callstring with -u and -p, or "
  echo "    You can put your credentials into a file called ~/.wget_cred with format:"
  echo "    ASF_username: your_username"
  echo "    ASF_password: your_password" 
  echo "Example: ./scihub_download_s1_sar.sh -i search_results1.txt -i search_results2.txt -u username -p password"
  echo "Please provide one or more input files."
  echo ""
  exit 1
fi


# Read the search results. It could be multiple calls of the -i flag. Read in the username and password. 
while getopts i:u:p: opt; do
    case $opt in
      i) multi+=("$OPTARG")
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
    esac
done
shift $((OPTIND -1))
# echo "The whole list of values is '${multi[@]}'"  # a debugging line. 


# Extracting ASF username and password from ~/.wget_cred in case they weren't provided
if [ -z "$username" ]; then
  echo "No ASF username provided in callstring; finding your ASF username in ~/.wget_cred"
  username=`grep 'ASF_username' ~/.wget_cred | awk {'print $2'}`
  echo "ASF username found: "$username
fi
if [ -z "$password" ]; then
  echo "No ASF password provided in callstring; finding your ASF password in ~/.wget_cred"
  password=`grep 'ASF_password' ~/.wget_cred | awk {'print $2'}`
  echo "ASF password found."
fi

# Get the username and password from the user
if [ -z "$username" ]; then
  echo "No ASF username has been found. Please read the documentation of this script for details. Exiting..."
  exit 1
fi
if [ -z "$password" ]; then
  echo "No ASF password has been found. Please read the documentation of this script for details. Exiting..."
  exit 1
fi


# Defining parameters
id_results=uuid_file.txt
am_i_linux=`uname -a | grep 'Linux'`
echo $am_i_linux

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
if [ ! -z "$am_i_linux" ]; then  # I am on a linux machine
	echo "I'm a linux"
    sed -i 's/<str name=\"uuid\">//g' $id_results
    sed -i 's/<title>//g' $id_results
    sed -i 's/<\/title>//g' $id_results
    sed -i 's/<\/str>//g' $id_results
else  # we are on a mac or non-linux machine! Macs need the -i '' in the sed call.
    echo "I'm a mac"
    sed -i '' 's/<str name=\"uuid\">//g' $id_results
    sed -i '' 's/<title>//g' $id_results
    sed -i '' 's/<\/title>//g' $id_results
    sed -i '' 's/<\/str>//g' $id_results
fi

counter=0
while read p; do
  if [ $counter = 0 ]; then
  	title=$p

    if [ ! -d DATA/$title.SAFE ]; then
    # In this version, I actually download from the ASF. The download goes about 5x faster than Copernicus for users in North America. 
    # You need your NASA/Earthdata/ASF credentials, either from the runstring or from your ~/.wget_cred file
      echo "Downloading DATA/"$title
      wget --http-user=$username --http-password=$password -c -O DATA/"$title".zip "https://datapool.asf.alaska.edu/SLC/SA/$title.zip"

      # cd DATA
      # unzip $title.zip
      # rm $title.SAFE/measurement/*-slc-vh-*.tiff
      # rm $title.zip
      # cd ../
    else
      echo "Already in the data directory: Skipping "$title
    fi
    counter=1
  	continue
  else
    uuid=$p
  fi
  #echo $title
  #echo $uuid
  
  counter=0
  
done <$id_results


  # I have the option to download the Manifest or Data from ESA. It's pretty slow though. 
  # MANIFEST only
  # wget --no-check-certificate --user= --password= -O MANIFEST/"$title"_manifest.safe "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/Nodes('$title.SAFE')/Nodes('manifest.safe')/\$value"

  # DATA FROM COPERNICUS (full thing- will take a long time)!
  # wget --no-check-certificate --user= --password= -O DATA/"$title".SAFE.zip "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value"
  # Takes a few hours for each SAFE.zip. 
  # Each one can be unzipped with unzip. 
  
