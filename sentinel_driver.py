
"""
	Kathryn's S1A SBAS driver.py
	For this toolbox of scripts to perform GMTSAR SBAS processing on S1A/S1B images, 
	your starting configuration should be:

	--DATA (has a pile of .SAFE directories, each 3.5 or 7.0 GB)
	--batch.config   (very important file!)
	--topo
	    >> dem.grd
	    >>  If you've already run README_proc once before, you will also have other files related to topographic phase:
		[S1A20150403_ALL_F1.LED, gmt.history, topo_ra.grd, trans.dat dem.grd, master.PRM, topo_ra.ps].
		These files take a long time to generate, and only need to be done once I think. 
		At that point you can keep these files and just start at phase 2 in the preproc_batch_tops file. 

	Will Make:
	--raw_orig
	    >> *.xml 
	    >> *.tiff
	    >> *.EOF
	    >> yyyymmdd_manifest.safe
	    >> s1a-aux-cal.xml  (for the satellite)

	In order to decide on your supermaster, you will first produce: 
	--README_prep.txt
	--data.in
	--raw/baseline_table.dat
	--raw/baseline.ps 

	For actual processing with your supermaster, you will then produce:
	--README_proc.txt
	--intf.in 
	--Network_Geometry.eps
	--raw/*
	--intf_all/*
	--SBAS/*

        You will control how the program operates through the processing stages specified in the batch.config file. 
"""

import sentinel_main_functions

if __name__=="__main__":

	# Step 0
	config_params = sentinel_main_functions.read_config();
	sentinel_main_functions.manifest2raw_orig_eof(config_params);

	# # Step 1: choose master and preprocess (Step 2 = aligning; combined for Sentinel)
	sentinel_main_functions.preprocess(config_params);

	# # Step 3
	sentinel_main_functions.topo2ra(config_params);
	# Will separate this later! 

	# # Step 4
	sentinel_main_functions.make_interferograms(config_params);




