# S1_batches

This set of codes I use for help with Sentinel-1 data analysis. Its functionality is basically down to searching the Sentinel-1 archive and downloading data.  I used to use it for stack processing with GMTSAR and SBAS. 

## Description

Disclaimer: Just use MintPy instead!  

I still use this library for basic things like searching the Sentinel-1 archive, reading GMTSAR and ISCE files into Python, and maybe doing a simple average, but for my future time series needs, I intend to switch to MintPy.

### Capabilities: 
* The scripts in ```Data_Search_and_Download/``` will search through SentinelHub for Sentinel-1 images and make maps of footprints and timing. 
* You will eventually be able to download all your search results too.

### Notes:
* [SciHub is deprecated as of October 2023](https://dataspace.copernicus.eu/news/2023-9-28-accessing-sentinel-mission-data-new-copernicus-data-space-ecosystem-apis).  January 2024: I switched the internals to the Sentinel-Hub OData API.  
* You will need a Copernicus login, [a Copernicus access token](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Overview/Authentication.html), and an ASF EarthData login to fully use these search and download features. 

### Usage of Data Search Features: 
An example usage of the scripts to search the Sentinel-1 database is given here. The directory containing these scripts (```Data_Search_and_Download/```) must be on your path. Type the name of any program by itself (with no arguments) to see a help menu.  If you want to download the .SAFE folders (work in progress), you should have an ASF login. 
```bash
s1_search_Odata.py
s1_search_Odata.py -s 2014-01-01 -e 2024-02-01 -c="-122.6/38.1" 
```

## Example: 

Map of the above search query:
![Footprint](https://github.com/kmaterna/S1_batches/blob/master/examples/footprints.png)

Time-plot of the above search query:
![Timing](https://github.com/kmaterna/S1_batches/blob/master/examples/timing.png)

