#!/bin/bash
# Driver for Descending Track 173 over the Salton Sea

latfile='../../geom_reference/lat.rdr'
lonfile='../../geom_reference/lon.rdr'

geocode_to_kmz() {   # arguments: filename, cmin, cmax, cmap
  # Step 1: geocode the interferogram file.  Will automatically name it geo_$name
  geocodeGdal.py -l $latfile -L $lonfile -f "$1" --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
  # -x and -y are the interval for pixel size, in degrees.  Default is 0.001. 0.0005 improves Google Earth experience.
  # Be aware: The output data in geo$name can be complex, but the vrt always says FLOAT32. We try to pass only floats for safety.

  # Step 2: Send geocoded intf into geotiff format, specifying the color scale
  isce2geotiff.py -i "geo_$1" -o "geo_$1.tif" -c "$2" "$3" -m "$4"
  #   -b BAND       Band number to use if input image is multiband. Default: 0
  #   -c CLIM CLIM  Color limits for the graphics
  #   -m CMAP       Matplotlib colormap to use
  # This returns only real(x) when passed a complex number. Could it be because of the vrt?

  # Step 3: Place intf into a valid KMZ file.
  gdal_translate -of KMLSUPEROVERLAY "geo_$1.tif" "geo_$1.kmz"
}

cd merged/interferograms/20230409_20230421/ || exit
isceintf_to_phase.py --isce_file fine.int --outfile phase.int
isceintf_to_phase.py --isce_file filt_fine.int --outfile filt_phase.int
mask_isce_file.py -i filt_phase.int -o filt_phase_masked.int -c filt_fine.cor -t 0.37  # implement coherence mask on phase
geocodeGdal.py -l $latfile -L $lonfile -f phase.int --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
geocodeGdal.py -l $latfile -L $lonfile -f filt_phase.int --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
geocodeGdal.py -l $latfile -L $lonfile -f filt_phase_masked.int --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
geocodeGdal.py -l $latfile -L $lonfile -f filt_fine.unw --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
geocodeGdal.py -l $latfile -L $lonfile -f filt_fine.cor --bbox '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
geocode_to_kmz filt_phase.int -3.14 3.14 rainbow
geocode_to_kmz phase.int -3.14 3.14 rainbow


# # Step 4: Geocode los.rdr file: EXECUTE from within ../../geom_reference with the same parameters as the interferograms.
cd ../../geom_reference/ || exit
geocodeGdal.py -l lat.rdr -L lon.rdr -f los.rdr -b '32.5 33.4 -116.15 -115.15' -x 0.0001 -y 0.0001
