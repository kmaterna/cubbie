

import sys
import glob
import numpy as np 
import netcdf_read_write as rwr
import test_outputs
import readmytupledata
import stack_corr

# myfile = sys.argv[1]
# print(myfile);
# test_outputs.plot_grid_file(myfile,"test");



filelist=glob.glob("stacking/unwrapped/*unwrap.grd");
# datacube = readmytupledata.reader(filelist);
# a=stack_corr.stack_corr(datacube, np.nan);
# rwr.produce_output_netcdf(datacube.xvalues, datacube.yvalues, a, 'Percentage', 'signalspread_please_test.nc')
# rwr.produce_output_plot('signalspread_please_test.nc', 'Signal Spread', 'signalspread_please_test.png', 
# 	'Percentage of coherence (out of 17 images)' )
rwr.produce_output_plot(filelist[0],'unwrapped_phase','single_intf.png','phase');
rwr.produce_output_contourf(filelist[0],'unwrapped_phase','single_intf_contourf.png','phase');
