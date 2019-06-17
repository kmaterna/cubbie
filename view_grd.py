

import sys
import netcdf_read_write
import test_outputs

myfile = sys.argv[1]
print(myfile);
test_outputs.plot_grid_file(myfile,"test");

# [xdata_p,ydata_p]=netcdf_read_write.read_grd_xy(orig_phasefile);
