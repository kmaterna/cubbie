# S1_batches

This set of scripts is built on top of GMTSAR to perform Sentinel-1 batch processing. It will search the Sentinel-1 archive, download data, organize files, select the super-master, produce interferograms with GMTSAR, and perform other options. 

## Description

Dislcaimer: this library is under development and is not finished.  Nonetheless, please feel free to read, fork, use, and/or contribute. 

### Capabilities: 
* The scripts in Data_Search_and_Download will search through SciHub for Sentinel-1 images and return the results in maps of footprints and plots of timing. 
* The scripts in the main directory perform batch processing for Sentinel-1. The driver is sentinel_driver.py and the controls are manipulated via a batch.config file
* The scripts in stacking_tools have an implementation of SBAS for stacked data
* Example config files are kept in the configs_and_setup directory

### Usage of Data Search Features: 
An example usage of the scripts to search the Sentinel-1 database is given here. The directory containing these scripts must be on your path. 
```bash
scihub_search_s1_data.sh -s 2015-01-01 -e 2020-07-01 -p -124/40.3 -d Descending -o 13
scihub_display_footprints.sh -i search_results.txt 
scihub_download_s1_sar.sh -i search_results.txt
```

## Example: 

A map of a search query:
![Footprint](https://github.com/kmaterna/S1_batches/blob/master/Data_Search_and_Download/MTJ_footprints.png)

A time-plot of a search query:
![Timing](https://github.com/kmaterna/S1_batches/blob/master/Data_Search_and_Download/MTJ_timing.png)


