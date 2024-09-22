#!/usr/bin/env python
# Specific driver for Brawley Project, 2020-2021

import numpy as np
import datetime as dt
from Tectonic_Utils.read_write import netcdf_read_write
from .. import stacking_utilities


def reshape_TS_into_standard(outdir, earlyfile, cofile, latefile, outfile):
    """
    This is not particular general. It's using hard-coded information about time axis
    and the timing of the earthquake.
    It's for track 26509.
    This function pastes together a pre-seismic, co-seismic, and post-seismic set of time series or jumps.
    On the same xy grid.
    """
    tolerance = 300  # Purposely killing all pixels above this value.
    print("Reshaping UAVSAR file into single TS File")
    [_tdata1, xdata1, ydata1, zdata1] = netcdf_read_write.read_3D_netcdf(earlyfile)
    [_xdata2, _ydata2, zdata2] = netcdf_read_write.read_netcdf3(cofile)
    [_tdata3, _xdata3, _ydata3, zdata3] = netcdf_read_write.read_3D_netcdf(latefile)
    print(np.shape(zdata1), np.shape(zdata2), np.shape(zdata3))
    ynum = np.shape(zdata1)[1]
    xnum = np.shape(zdata1)[2]
    znum = np.shape(zdata1)[0] + np.shape(zdata3)[0]
    total_data = np.zeros([znum, ynum, xnum])
    print(np.shape(total_data))
    for i in range(np.shape(zdata1)[0]):
        temp = zdata1[i, :, :]
        temp[abs(temp) > tolerance] = np.nan  # killing outliers
        total_data[i, :, :] = temp  # doing this six times.
        print("Early data Slice %d" % i)
    temp = zdata2
    temp[abs(temp) > tolerance] = np.nan
    total_data[7, :, :] = np.add(total_data[6, :, :], zdata2)  # The coseismic chunk
    print("Coseismic data slice 7")
    for i in range(1, np.shape(zdata3)[0]):
        temp = zdata3[i, :, :]
        temp[abs(temp) > tolerance] = np.nan  # killing outliers
        total_data[i + np.shape(zdata1)[0], :, :] = np.add(temp, total_data[7, :, :])
        print("Postseismic data slice %d " % (i + np.shape(zdata1)[0]))
    dtarray = []
    dtarray.append(dt.datetime.strptime("2009-04-24", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2009-09-21", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2010-04-12", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2010-07-01", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2010-12-01", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2011-05-18", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2011-11-10", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2012-09-26", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2013-05-24", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2014-06-11", "%Y-%m-%d"))  # Hard-coded
    dtarray.append(dt.datetime.strptime("2017-11-01", "%Y-%m-%d"))  # Hard-coded
    zunits = "mm"
    print(np.shape(total_data))
    netcdf_read_write.produce_output_timeseries(xdata1, ydata1, total_data, dtarray, zunits, outfile)
    stacking_utilities.plot_full_timeseries(outfile, dtarray, outdir + "TS_cumulative.png", vmin=-50, vmax=200,
                                            aspect=1 / 8)
    stacking_utilities.plot_incremental_timeseries(outfile, dtarray, outdir + "TS_incremental.png", vmin=-50, vmax=100,
                                                   aspect=1 / 8)
    return


if __name__ == "__main__":
    outdir = "nsbas_apr20_total_18_14_1.5/"
    earlyfile = "nsbas_apr20_early_18_14_1.5/TS.nc"
    cofile = "nsbas_apr20_coseismic_18_14_1.5/coseismic.grd"
    latefile = "nsbas_apr20_late_18_14_1.5/TS.nc"
    outfile = outdir + "TS.nc"
    reshape_TS_into_standard(outdir, earlyfile, cofile, latefile, outfile)
