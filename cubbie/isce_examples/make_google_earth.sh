#!/bin/bash
# Create some KMZ files so that an intern to perform data quality measurements. 
# For track D173 in Salton Sea, bbox is '32.7 33.5 -115.9 -115.0'
# For Nevada: '40.65 41.2 -118.5 -117.7'
# For Iceland: '63.80 64.10 -22.73 -22.03' with -x 0.00005 -y 0.00005
# For Alaska: '55.0 56.4 -162.0 -158.0'
# For Alaska 2: '55.6 56.0 -155.8 -155.2'
# For Imperial: '32.5 32.9 -115.7 -115.3' with -x 0.0001 -y 0.0001
# 6/30/2022

latfile="../../geom_reference/lat.rdr"
lonfile="../../geom_reference/lon.rdr"
losrdr="../../geom_reference/los.rdr"

geocode_to_kmz() {   # arguments: filename, cmin, cmax, cmap
  # Step 1: geocode the interferogram file.  Will automatically name it geo_$name
  geocodeGdal.py -l $latfile -L $lonfile -f $1 --bbox '63.80 64.10 -22.73 -22.03' -x 0.00005 -y 0.00005
  # geocode.py $1 --lat-file $latfile --lon-file $lonfile
  # -x and -y are the interval for pixel size, in degrees.  Default is 0.001. 0.0005 improves Google Earth experience. 
  # Be aware: The output data in geo$name can be complex, but the vrt always says FLOAT32. We try to pass only floats for safety. 

  gdal_edit.py -a_nodata -9999 geo_$1   # teach the file that nodatavalue = -9999  (useful for displaying as tiff, but not for GRDs)

  # Step 2: Send geocoded intf into geotiff format, specifying the color scale
  isce2geotiff.py -i geo_$1 -o geo_$1.tif -c $2 $3 -m $4
  #   -b BAND       Band number to use if input image is multiband. Default: 0
  #   -c CLIM CLIM  Color limits for the graphics
  #   -m CMAP       Matplotlib colormap to use
  # This returns only real(x) when passed a complex number. Could it be because of the vrt? 

  # Step 3: Place intf into a valid KMZ file.
  gdal_translate -of KMLSUPEROVERLAY geo_$1.tif geo_$1.kmz
}

geocode_to_kmz_custom_table() {   # arguments: filename, cmin, cmax, cmap, for masked stuff the tables already exist
  geocodeGdal.py -l $latfile -L $lonfile -f $1 --bbox '63.80 64.10 -22.73 -22.03' -x 0.00005 -y 0.00005  # topsStack version
  # geocode.py $1 --lat-file $latfile --lon-file $lonfile   # general ISCE version
  gdal_edit.py -a_nodata -9999 geo_$1   # teach the file that nodatavalue = -9999  (useful for displaying as tiff, but not for GRDs)
  isce2geotiff.py -i geo_$1 -o geo_$1.tif -c $2 $3 -t $4    #   -t table    Color table to use, customized to have low values masked
  gdal_translate -of KMLSUPEROVERLAY geo_$1.tif geo_$1.kmz
}


rm geo_*
unw_to_phase.py --unw_file filt_fine.unw --outfile filt_fine_single.unw  # user-written script
isceintf_to_phase.py --isce_file filt_fine.int --outfile phase.int  # user-written script
mask_isce_file.py -i filt_fine_single.unw -o filt_fine_single_masked.unw -c filt_fine.cor -t 0.30 -v -9999 # implement coherence mask on phase
mask_isce_file.py -i phase.int -o phase_masked.int -c filt_fine.cor -t 0.30  -v -9999 # implement coherence mask on phase
gdal_edit.py -a_nodata -9999 filt_fine_single_masked.unw   # teach the file that nodatavalue = -9999
gdal_edit.py -a_nodata -9999 phase_masked.int   # teach the file that nodatavalue = -9999

geocode_to_kmz phase.int -3.14 3.14 gist_rainbow
geocode_to_kmz_custom_table phase_masked.int -3.14 3.14 gist_rainbow.cpt
geocode_to_kmz filt_fine_single.unw -15 15 RdYlBu
geocode_to_kmz_custom_table filt_fine_single_masked.unw -15 15 RdYlBu.cpt
geocode_to_kmz filt_fine.cor 0 1 gray


# # Step 4: Geocode los.rdr file: EXECUTE from within ../../geom_reference with the same parameters as the interferograms.
cd ../../geom_reference/
rm geo_*
geocodeGdal.py -l lat.rdr -L lon.rdr -f los.rdr -b '63.80 64.10 -22.73 -22.03' -x 0.00005 -y 0.00005

