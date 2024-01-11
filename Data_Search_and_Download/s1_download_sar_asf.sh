#!/bin/bash
# Download data as shown in the search_results. 
# Feb. 14, 2018.

if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script downloads the results of data queries and puts it into a directory called DATA/"
  echo "one level below where the script is called."
  echo "Usage: ./s1_download_sar_asf.sh -options"
  echo "  Input options:"
  echo " -i input_file"
  echo " -u username"
  echo " -p password"
  echo " -z unzip_flag. Will unzip results if set [default 0]"
  echo "Note: To use ASF, you will need a NASA Earthdata login."
  echo "    You can type them directly into the callstring with -u and -p, or "
  echo "    You can put your credentials into a file called ~/.wget_cred with format:"
  echo "    ASF_username: your_username"
  echo "    ASF_password: your_password" 
  echo "Example: ./s1_download_sar_asf.sh -i search_results1.txt -u username -p password"
  echo "Please provide an input file."
  echo ""
  exit 1
fi


# Read the search results. It could be multiple calls of the -i flag. Read in the username and password.
unzip_flag=0
while getopts i:u:p: opt; do
    case $opt in
      i)
        echo "-i was triggered, parameter: $OPTARG" >&2
        input_file=$OPTARG
        ;;
      u)  # username
        echo "-u was triggered, parameter: $OPTARG" >&2
        username=$OPTARG
        ;;
      p)  # password
        echo "-p was triggered, parameter: $OPTARG" >&2
        password=$OPTARG
        ;;
      z)  # unzip_flag
        echo "-z was triggered, parameter: $OPTARG" >&2
        unzip_flag=$OPTARG
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


# Where will the data live? 
mkdir -p DATA
id_results=ids.txt
awk '(NR>1)' $input_file | awk '{print $1}' > $id_results  # extract the .SAFE names from the query results

counter=0
while read p; do
  if [ $counter = 0 ]; then
  	title=$p

    if [ ! -d DATA/$title ]; then
    # In this version, I actually download from the ASF. The download goes about 5x faster than Copernicus for users in North America. 
    # You need your NASA/Earthdata/ASF credentials, either from the runstring or from your ~/.wget_cred file
      echo "Downloading DATA/"$title
      wget --http-user=$username --http-password=$password -c -O DATA/"$title".zip "https://datapool.asf.alaska.edu/SLC/SA/$title.zip"

      if [ $unzip_flag != 0 ]; then
         cd DATA
         unzip $title.zip
         rm $title.SAFE/measurement/*-slc-vh-*.tiff
         rm $title.zip
         cd ../
      fi
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
rm $id_results
