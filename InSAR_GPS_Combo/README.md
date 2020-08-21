# S1_batches/InSAR_GPS_Combo

When you're ready to compare your InSAR with GNSS velocities, this set of functions can help.  To read GNSS velocities, it depends on a GNSS repo (currently called Mendocino_Geodesy), so make sure to find that and put it on your path. 

## Description

Dislcaimer: this library is under development and is frequently updated.  Nonetheless, please feel free to read, fork, use, and/or contribute. 

## Example
Here we perform a simple calculation to interpolate GNSS velocities across a study region and then project that calculation simply into LOS veiwing geometries. We assume a single incidence angle and flight angle for the whole field, which is a simplification of the actual picture (the incidence angle changes from near range to far range). 

GNSS Velocities projected into ascending and descending viewing geometries, and interpolated:
![Footprint](https://github.com/kmaterna/S1_batches/blob/master/InSAR_GPS_Combo/testing_and_results/LOS.png)
