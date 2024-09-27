#!/bin/bash

# SCRIPT TO CHECK DOWNLOAD STATUS
# This will tell us if certain search results did not successfully get downloaded. Will have to fix these some other way.  
# K. Materna, March. 2018

# Parse inputs
if [[ "$#" -eq 0 ]]; then
  echo ""
  echo "This script checks if the download in the given file actually completed"
  echo "Usage: ./scihub_check_dl_status.sh -options"
  echo "Example: ./scihub_check_dl_status.sh -i search_results.txt -d DATA"
  echo "-i: the search_results file that contains the files we want to download"
  echo "-d: the relative path to the place where the .SAFE directories live."
  echo ""
  echo "  outputs:"
  echo "    Will send output to \"undownloaded.txt\"."
  echo ""
  echo ""
  exit 1
fi


while getopts :i:d: opt; do
  case $opt in
  	i)  # the file with search results
      echo "-i was triggered, parameter: $OPTARG" >&2
      search_file=$OPTARG
      ;;
    d)  # the directory where data lives
      echo "-d was triggered, parameter: $OPTARG" >&2
      relpath=$OPTARG
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

# Setting up Parameters and temp files 
safe_results="safe_results.txt"
undownloaded="undownloaded.txt"
search_directory=`pwd`/$relpath
grep "\.SAFE" $search_file > $safe_results


# List the .SAFE files that we expect
am_i_linux=`uname -a | grep 'Linux'`
echo $am_i_linux
# the -i '' is because of mac computers. Might need to delete the '' on a linux machine. 
if [ ! -z "$am_i_linux" ]; then  # I am on a linux machine
    sed -i 's/<str name=\"filename\">//g' $safe_results
    sed -i 's/<\/str>//g' $safe_results
else  # we are on a mac or non-linux machine! Macs need the -i '' in the sed call.
    sed -i '' 's/<str name=\"filename\">//g' $safe_results
    sed -i '' 's/<\/str>//g' $safe_results
fi

# Decide which files we don't have already
echo "# Files we failed to download:" > $undownloaded
echo "The files we should download, i.e. are not in the DATA directory: "
while IFS= read -r var
do
  if [ ! -d $relpath/"$var" ]; then  # if the files that we searched are NOT in the DATA directory after downloading, then we have to chase them down some other way. 
   	echo "$var" >> $undownloaded
    echo "do not have file: "$relpath/"$var"
    # echo $relpath/"$var"
  fi
done < "$safe_results"

# Print outputs
num_of_search_file=`cat $safe_results | wc -l`
num_of_current_file=`ls $relpath | grep 'SAFE' | wc -l`
num_of_missing_file=`grep 'SAFE' $undownloaded | wc -l`
echo "RESULTS"
echo ">> Input $search_file returned "$num_of_search_file" SAFE results"
echo ">> Data directory $relpath contains "$num_of_current_file" SAFE results"
echo ">> We must download "$num_of_missing_file" SAFE results"

rm $safe_results


