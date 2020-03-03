#!/usr/bin/python
import numpy as np
import collections
from datetime import datetime
import re
import netcdf_read_write as rwr
import isce_read_write

data = collections.namedtuple('data', ['filepaths', 'dates_correct', 'date_deltas',  'xvalues', 'yvalues', 'zvalues'])

def reader(filepathslist):
    """This function takes in a list of filepaths that each contain a 2d array of data, effectively taking
    in a cuboid of data. It splits and stores this data in a named tuple which is returned. This can then be used
    to extract key pieces of information."""
    filepaths  = []
    dates_correct , date_deltas = [], []
    xvalues, yvalues, zvalues = [], [], []
    for i in range(len(filepathslist)):
        print(filepathslist[i])
        filepaths.append(filepathslist[i])
        datesplit = filepathslist[i].split('/')[-1]  # example: 2015157_2018177_unwrap.grd
        date_new = datesplit.replace(datesplit[0:7], str(int(datesplit[0:7]) + 1))
        date_new = date_new.replace(date_new[8:15], str(int(date_new[8:15]) + 1))  # adding 1 to the date because 000 = January 1
        dates_correct.append(date_new[0:15])  # example: 2015158_2018178

        delta = abs(datetime.strptime(dates_correct[i][0:7], '%Y%j') - datetime.strptime(dates_correct[i][8:15], '%Y%j'))  # timedelta object
        date_deltas.append(delta.days/365.24)  # in years. Is that a good idea? 

        xdata, ydata, zdata = rwr.read_grd_xyz(filepathslist[i])
        xvalues=xdata
        yvalues=ydata
        zvalues.append(zdata)
        if i == round(len(filepathslist)/2):
            print('halfway done reading files...')

    mydata = data(filepaths=np.array(filepaths), dates_correct=np.array(dates_correct), 
        date_deltas=np.array(date_deltas), xvalues=np.array(xvalues), yvalues=np.array(yvalues), zvalues=np.array(zvalues))
    return mydata



def reader_isce(filepathslist, band=1):
    """This function takes in a list of filepaths that each contain a 2d array of data, effectively taking
    in a cuboid of data. It splits and stores this data in a named tuple which is returned. This can then be used
    to extract key pieces of information. It reads in ISCE format. """

    filepaths  = []
    dates_correct , date_deltas = [], []
    xvalues, yvalues, zvalues = [], [], []
    for i in range(len(filepathslist)):
        filepaths.append(filepathslist[i])
        # In the case of ISCE, we have the dates in YYYYMMDD_YYYYMMDD format somewhere within the filepath (maybe multiple times). We take the first. 
        datesplit = re.findall(r"\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\d\d", filepathslist[i])[0]; #  example: 20100402_20140304
        date1 = datetime.strptime(datesplit[0:8],"%Y%m%d");
        date2 = datetime.strptime(datesplit[9:17],"%Y%m%d");
        datestr_julian=datetime.strftime(date1,"%Y%j")+"_"+datetime.strftime(date2,"%Y%j");  # in order to maintain consistency with GMTSAR formats
        dates_correct.append(datestr_julian)  # example: 2015158_2018178
        delta = abs(date1-date2)
        date_deltas.append(delta.days/365.24)  # in years. 

        zdata = isce_read_write.read_scalar_data(filepathslist[i], band);  # NOTE: For unwrapped files, this will be band=2
        xvalues=range(0,np.shape(zdata)[1]);  # is this correct? 
        yvalues=range(0,np.shape(zdata)[0]);
        zvalues.append(zdata)
        if i == round(len(filepathslist)/2):
            print('halfway done reading files...')

    mydata = data(filepaths=np.array(filepaths), dates_correct=np.array(dates_correct), 
        date_deltas=np.array(date_deltas), xvalues=np.array(xvalues), yvalues=np.array(yvalues), zvalues=np.array(zvalues))

    return mydata; 