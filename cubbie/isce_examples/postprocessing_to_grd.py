#!/usr/bin/env python

import s1_batches.read_write_insar_utilities.isce_read_write as rw
from s1_batches.math_tools import grid_tools
from Tectonic_Utils.read_write import netcdf_read_write
import numpy as np
import matplotlib.pyplot as plt

paramdict = {"imperial_range": (-115.60, -115.35, 32.75, 32.95),
             "filelist": ('geo_filt_fine.unw', 'geo_phase.int', 'geo_filt_phase.int', 'geo_filt_fine.cor',
                          'geo_filt_phase_masked.int'),
             "imperial_intfs": ("20230328_20230409", "20230409_20230421")}


def clip_isce_files_into_grds(directory, filelist=(), cut_range=(),
                              los_inc_outfile=None, los_az_outfile=None):
    """
    Clip a few isce files and convert them into grds. Also convert LOS information too.

    :param directory: string
    :param filelist: list of strings
    :param cut_range: (lonW, lonE, latS, latN)
    :param los_inc_outfile: string, filename
    :param los_az_outfile: string, filename
    """
    if "geo_filt_fine.unw" in filelist:
        x_orig, y_orig, array1 = rw.read_isce_unw_geo(directory+'geo_filt_fine.unw')
        x, y, array1 = grid_tools.clip_array_by_bbox(x_orig, y_orig, array1, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, array1, zunits='unw_phase',
                                                netcdfname=directory+'geo_filt_fine_unw.grd')
    if "geo_phase.int" in filelist:
        x_orig, y_orig, array1 = rw.read_scalar_data(directory+'geo_phase.int')
        x, y, array1 = grid_tools.clip_array_by_bbox(x_orig, y_orig, array1, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, array1, zunits='wr_phase', netcdfname=directory+'geo_phase.grd')
    if "geo_filt_phase.int" in filelist:
        x_orig, y_orig, array1 = rw.read_scalar_data(directory+'geo_filt_phase.int')
        x, y, array1 = grid_tools.clip_array_by_bbox(x_orig, y_orig, array1, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, array1, zunits='wr_phase',
                                                netcdfname=directory+'geo_filt_phase.grd')
    if "geo_filt_phase_masked.int" in filelist:
        x_orig, y_orig, array1 = rw.read_scalar_data(directory+'geo_filt_phase_masked.int')
        x, y, array1 = grid_tools.clip_array_by_bbox(x_orig, y_orig, array1, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, array1, zunits='wr_phase',
                                                netcdfname=directory+'geo_filt_phase_masked.grd')
    if "geo_filt_fine.cor" in filelist:
        x_orig, y_orig, array1 = rw.read_scalar_data(directory+'geo_filt_fine.cor')
        x, y, array1 = grid_tools.clip_array_by_bbox(x_orig, y_orig, array1, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, array1, zunits='cor', netcdfname=directory+'geo_cor.grd')
    if "geo_los.rdr" in filelist:
        # Read the LOS information from ISCE
        x_orig, y_orig, azimuth_array = rw.read_scalar_data(directory+"geo_los.rdr", band=2)
        _, _, inc_array = rw.read_scalar_data(directory+"geo_los.rdr", band=1)
        x, y, inc_array_cut = grid_tools.clip_array_by_bbox(x_orig, y_orig, inc_array, cut_range)
        x, y, azimuth_array_cut = grid_tools.clip_array_by_bbox(x_orig, y_orig, azimuth_array, cut_range)
        netcdf_read_write.produce_output_netcdf(x, y, inc_array_cut, zunits='inc', netcdfname=los_inc_outfile)
        netcdf_read_write.produce_output_netcdf(x, y, azimuth_array_cut, zunits='azimuth', netcdfname=los_az_outfile)
    return


if __name__ == "__main__":
    # Postprocess the Imperial images
    for intf in paramdict["imperial_intfs"]:
        clip_isce_files_into_grds(directory="merged/interferograms/"+intf+"/", filelist=paramdict['filelist'],
                                  cut_range=paramdict['imperial_range'])
    clip_isce_files_into_grds(directory="merged/geom_reference/", filelist='geo_los.rdr',
                              cut_range=paramdict['imperial_range'],
                              los_inc_outfile='merged/geom_reference/imperial_los_inc.grd',
                              los_az_outfile='merged/geom_reference/imperial_los_az.grd')
