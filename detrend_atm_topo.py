# A script to take a stack of images plus a DEM
# Solve for best-fitting linear trend
# (Globally or locally)
# Remove trend and save the adjusted stack in out_dir



def main_function(staging_directory, outdir):
	configure(staging_directory, outdir);
	return;


def configure(staging_directory, outdir):
	print("Detrending atm/topo for all files in %s and storing the result in %s " % (staging_directory, outdir) )
	return;

def inputs():
	return;
def compute():
	return;
def outputs():
	return;



if __name__=="__main__":
	main_function('intf_all/unwrap.grd', 'intf_all/atm_topo_corrected.grd');