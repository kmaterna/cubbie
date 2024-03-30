import numpy as np
import collections
import re
from datetime import datetime
from s1_batches.read_write_insar_utilities import isce_read_write
from Tectonic_Utils.read_write import netcdf_read_write as rwr
from . import stacking_utilities


"""
This should be re-written to be a list of objects instead of an object of lists. 
When I have time. 
"""

data = collections.namedtuple('data', ['filepaths', 'date_pairs_julian', 'date_deltas',
                                       'xvalues', 'yvalues', 'zvalues', 'date_pairs_dt', 'ts_dates'])


def reader(filepathslist):
    """
    This function takes in a list of filepaths to GMTSAR grd files, taking in a cuboid of data.
    It splits and returns this data in a named tuple.
    """
    filepaths = []
    date_pairs_julian, date_deltas, date_pairs = [], [], []
    xdata, ydata, zvalues = [], [], []
    for i in range(len(filepathslist)):
        print(filepathslist[i])
        # Establish timing and filepath information
        filepaths.append(filepathslist[i])
        datesplit = re.findall(r"\d\d\d\d\d\d\d_\d\d\d\d\d\d\d", filepathslist[i])[0]  # example: 2010040_2014052
        # adding 1 to both dates because 000 = January 1
        date_new = datesplit.replace(datesplit[0:7], str(int(datesplit[0:7]) + 1))  # replacing first date
        date_new = date_new.replace(date_new[8:15], str(int(date_new[8:15]) + 1))  # replacing second date
        date_pairs_julian.append(date_new[0:15])  # example: 2015158_2018178
        acq1 = datetime.strptime(date_new[0:7], '%Y%j')
        acq2 = datetime.strptime(date_new[8:15], '%Y%j')
        date_pairs.append([acq1, acq2])
        delta = abs(acq1 - acq2)  # timedelta object
        date_deltas.append(delta.days / 365.24)  # in years. 

        # Read in the data
        xdata, ydata, zdata = rwr.read_netcdf4(filepathslist[i])  # does this work on netcdf3 as well?
        zvalues.append(zdata)
        if i == round(len(filepathslist) / 2):
            print('halfway done reading files...')

    # The sorted list of dates used in this interferogram network
    ts_dates = stacking_utilities.get_unique_dts_from_intf_dates(np.array(date_pairs))

    mydata = data(filepaths=np.array(filepaths), date_pairs_julian=np.array(date_pairs_julian),
                  date_deltas=np.array(date_deltas), xvalues=np.array(xdata), yvalues=np.array(ydata),
                  zvalues=np.array(zvalues), date_pairs_dt=np.array(date_pairs), ts_dates=ts_dates)
    return mydata


def reader_from_ts(filepathslist):
    """ 
    This function makes a tuple of grids in timesteps
    It can read in radar coords or geocoded coords, depending on the use of xvar, yvar
    """
    filepaths, zvalues, ts_dates = [], [], []
    xvalues, yvalues = [], []
    for i in range(len(filepathslist)):
        print(filepathslist[i])
        # Establish timing and filepath information
        filepaths.append(filepathslist[i])
        datestr = re.findall(r"\d\d\d\d\d\d\d\d", filepathslist[i])[0]
        ts_dates.append(datetime.strptime(datestr, "%Y%m%d"))
        # Read in the data, either netcdf3 or netcdf4
        [xvalues, yvalues, zdata] = rwr.read_netcdf4(filepathslist[i])
        zvalues.append(zdata)
        if i == round(len(filepathslist) / 2):
            print('halfway done reading files...')
    mydata = data(filepaths=np.array(filepaths), date_pairs_julian=None, date_deltas=None,
                  xvalues=np.array(xvalues), yvalues=np.array(yvalues), zvalues=np.array(zvalues),
                  date_pairs_dt=None, ts_dates=np.array(ts_dates))
    return mydata


def reader_simple_format(file_names):
    """
    An earlier reading function, works fast, useful for things like coherence statistics
    """
    filename = file_names[0]
    [xdata, ydata] = rwr.read_netcdf3(filename)[0:2]
    data_all = []
    for ifile in file_names:  # this happens to be in date order on my mac
        data = rwr.read_netcdf3(ifile)[2]
        data_all.append(data)
    date_pairs = []
    for name in file_names:
        pairname = name.split('/')[-2][0:15]
        date_pairs.append(pairname)  # returning something like '2016292_2016316' for each intf
        print(pairname)
    return [xdata, ydata, data_all, date_pairs]


def reader_isce(filepathslist, band=1):
    """
    This function takes in a list of filepaths that each contain a 2d array of data, taking
    in a cuboid of data. It splits and stores this data in a named tuple which is returned. This can then be used
    to extract key pieces of information. It reads in ISCE format. 
    """

    filepaths = []
    date_pairs_julian, date_deltas, date_pairs = [], [], []
    xvalues, yvalues, zvalues = [], [], []
    for i in range(len(filepathslist)):
        filepaths.append(filepathslist[i])
        # In the case of ISCE, we have the dates in YYYYMMDD_YYYYMMDD format somewhere within the filepath
        # (maybe multiple times). We take the first.
        datesplit = re.findall(r"\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\d\d", filepathslist[i])[0]  # example: 20100402_20140304
        date1 = datetime.strptime(datesplit[0:8], "%Y%m%d")
        date2 = datetime.strptime(datesplit[9:17], "%Y%m%d")
        date_pairs.append([date1, date2])
        # in order to maintain consistency with GMTSAR formats:
        datestr_julian = datetime.strftime(date1, "%Y%j") + "_" + datetime.strftime(date2, "%Y%j")
        date_pairs_julian.append(datestr_julian)  # example: 2015158_2018178
        delta = abs(date1 - date2)
        date_deltas.append(delta.days / 365.24)  # in years.

        _, _, zdata = isce_read_write.read_scalar_data(filepathslist[i], band,
                                                       flush_zeros=False)  # NOTE: For unwrapped files, will be band=2
        # flush_zeros=False preserves the zeros in the input datasets. Added April 9 2020. uncertain results.
        xvalues = range(0, np.shape(zdata)[1])
        yvalues = range(0, np.shape(zdata)[0])
        zvalues.append(zdata)
        if i == round(len(filepathslist) / 2):
            print('halfway done reading files...')

    # The sorted list of dates used in this interferogram network
    ts_dates = stacking_utilities.get_unique_dts_from_intf_dates(np.array(date_pairs))

    mydata = data(filepaths=np.array(filepaths), date_pairs_julian=np.array(date_pairs_julian),
                  date_deltas=np.array(date_deltas), xvalues=np.array(xvalues), yvalues=np.array(yvalues),
                  zvalues=np.array(zvalues), date_pairs_dt=np.array(date_pairs), ts_dates=ts_dates)

    return mydata
