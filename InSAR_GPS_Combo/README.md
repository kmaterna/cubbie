# S1_batches/InSAR_GPS_Combo

When you're ready to compare your InSAR with GNSS velocities, this set of functions can help.  To read GNSS velocities, it depends on a GNSS repo (GNSS_Timeseries_Veiwers), so make sure to find that and put it on your path. 

## Description

Dislcaimer: this library is under development and is frequently updated.  Nonetheless, please feel free to read, fork, use, and/or contribute. 

## Example
Here we perform a simple calculation to interpolate GNSS velocities across a study region and then project that calculation simply into LOS veiwing geometries. We assume a single incidence angle and flight angle for the whole field, which is a simplification of the actual picture (the incidence angle changes from near range to far range). 

GNSS Velocities projected into ascending and descending viewing geometries, and interpolated:
![Footprint](https://github.com/kmaterna/S1_batches/blob/master/InSAR_GPS_Combo/testing_and_results/LOS.png)

## Example
In another example, we can project a GAGE-provided 3D velocity field (e.g., from ftp://data-out.unavco.org/pub/products/velocity/) into LOS using grids of variable look vectors. This code depends on the S1_batches repository and the GNSS_Timeseries_Viewers repository for read functions. Configuration and driver is shown below. 
```python
config_params = {
	"gps_file":"../../GPS_POS_DATA/Velocity_Files/NAM08_pbovelfile_feb2018.txt",
	"look_vector_files":[
		"../Data/LookVector/SIOX_D071_east_look.grd",
		"../Data/LookVector/SIOX_D071_north_look.grd",
		"../Data/LookVector/SIOX_D071_up_look.grd"],
	"reference_gps":"P617",
	"coordbox_gps":[-119.2, -115.8, 33, 36.5],
	"outfile":"Track71_gps_ll_enu_los.txt"
};
calc_gps_LOS_var_incidence.top_level_driver(config_params);
```
