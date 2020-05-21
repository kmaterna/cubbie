# Testing the compatibility of the read/write functions
# In ISCE and GMTSAR
# March 2020

import isce_read_write
import netcdf_read_write


def test_read_write(filename):
	# A TEST OF READ/WRITE FUNCTIONS
	# Step 1: Read ISCE interferogram as phase
	# Step 2: Write it as ISCE format data
	# Step 3: Read ISCE interferogram again
	# Step 4: Write as .grd
	# Step 5: make plot. 
	# RESULT: PRETTY GOOD! There is a flipup between the ISCE/GMTSAR conventions, 
	# but it might never really be a problem. 

	# Step 1: read phase
	slc = isce_read_write.read_phase_data(filename);
	print("Shape of slc is ", np.shape(slc) );
	isce_read_write.plot_complex_data(filename, aspect=1/10, outname="original.png");

	# Step 2: write ISCE data file
	ny, nx = np.shape(slc);
	# dtype = 'FLOAT';
	isce_written="isce_written_phase.phase"
	isce_read_write.write_isce_data(slc, nx, ny, dtype="FLOAT", filename=isce_written);

	# Step 3: read phase again. 
	phase = isce_read_write.read_scalar_data(isce_written+".vrt");
	print("Shape of phase is ", np.shape(phase) );
	isce_read_write.plot_scalar_data("isce_written_phase.phase",colormap='rainbow', aspect=1/10, outname="isce_written_phase.png");

	# Step 4: write that phase as grd. 
	netcdfname = "netcdf_written_phase.grd"
	xdata=np.arange(0,nx);
	ydata=np.arange(0,ny);
	phase=np.flipud(phase);  # THIS SEEMS TO BE NECESSARY TO SWITCH BETWEEN CONVENTIONS. GRD PLOTS ARE UPSIDE DOWN FROM ISCE. 
	netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, "radians", netcdfname);

	# Step 5: look at what's inside; 
	netcdf_read_write.produce_output_plot(netcdfname, "phase", "grdstyle_phase.png", "radians", aspect=1/10);

	return;


if __name__=="__main__":
	test_read_write("temp.int");