# Sentinel Utilities

from subprocess import call
import os
import sys
import glob
import datetime
import numpy as np


def get_all_xml_names(directory, polarization, swath):
    pathname1=directory+"/*-vv-*-00"+swath+".xml";
    pathname2=directory+"/*-vv-*-00"+str(int(swath)+3)+".xml";
    list_of_images_temp=glob.glob(pathname1)+glob.glob(pathname2);
    list_of_images=[]
    for item in list_of_images_temp:
        list_of_images.append(item[:])
    return list_of_images;

def get_manifest_safe_names(directory):
    mansafe=glob.glob(directory+'/manifest.safe');
    return mansafe;

def get_all_tiff_names(directory, polarization, swath):
    pathname1=directory+"/*-vv-*-00"+swath+".tiff";
    pathname2=directory+"/*-vv-*-00"+str(int(swath)+3)+".tiff";
    list_of_images_temp=glob.glob(pathname1)+glob.glob(pathname2);
    list_of_images=[]
    for item in list_of_images_temp:
        list_of_images.append(item[:])
    return list_of_images;


def get_previous_and_following_day(datestring):
    """ This is a function that takes a date like 20160827 and generates 
    [20160826, 20160828]: the day before and the day after the date in question. """
    year=int(datestring[0:4]);
    month=int(datestring[4:6]);
    day=int(datestring[6:8]);
    mydate=datetime.date(year, month, day);
    tomorrow =mydate + datetime.timedelta(days=1);
    yesterday=mydate - datetime.timedelta(days=1);
    previous_day=pad_string_zeros(yesterday.year)+pad_string_zeros(yesterday.month)+pad_string_zeros(yesterday.day);
    following_day=pad_string_zeros(tomorrow.year)+pad_string_zeros(tomorrow.month)+pad_string_zeros(tomorrow.day);
    return [previous_day, following_day];

def get_date_from_xml(xml_name):
    """
    xml file has name like s1a-iw1-slc-vv-20150121t134413-20150121t134424-004270-005317-001.xml
    We want to return 20150121. 
    """
    xml_name=xml_name.split('/')[-1];
    mydate=xml_name[15:23];
    return mydate;

def pad_string_zeros(num):
    if num<10:
        numstring="0"+str(num);
    else:
        numstring=str(num);
    return numstring;

def get_eof_from_xml(xml_name, eof_dir):
    """ This returns something like S1A_OPER_AUX_POEORB_OPOD_20160930T122957_V20160909T225943_20160911T005943.EOF. 
    """
    mydate=get_date_from_xml(xml_name);
    [previous_day,following_day]=get_previous_and_following_day(mydate);
    eof_name=glob.glob(eof_dir+"/*"+previous_day+"*"+following_day+"*.EOF"); 
    eof_name=eof_name[0];
    return eof_name;


def make_data_in(polarization, swath, master_date="00000000"):
    """
    data.in is a reference table that links the xml file with the correct orbit file.
    """
    list_of_images=get_all_xml_names("raw_orig",polarization,swath);
    outfile=open("data.in",'w');
    if master_date=="00000000":
        for item in list_of_images:
            item=item.split("/")[-1];  # getting rid of the directory
            eof_name=get_eof_from_xml(item,"raw_orig");
            outfile.write(item[:-4]+":"+eof_name.split("/")[-1]+"\n");
    else:
        # write the master date first. 
        for item in list_of_images:
            mydate=get_date_from_xml(item);
            if mydate==master_date:
                item=item.split("/")[-1];  # getting rid of the directory
                eof_name=get_eof_from_xml(item,"raw_orig");
                outfile.write(item[:-4]+":"+eof_name.split("/")[-1]+"\n");
        # then write the other dates. 
        for item in list_of_images:
            mydate = get_date_from_xml(item);
            if mydate != master_date:
                item=item.split("/")[-1];  # getting rid of the directory
                eof_name=get_eof_from_xml(item,"raw_orig");
                outfile.write(item[:-4]+":"+eof_name.split("/")[-1]+"\n");
    outfile.close();
    print "data.in successfully printed."
    return;


# after running the baseline calculation from the first pre_proc_batch, choose a new master that is close to the median baseline and timespan.
def choose_master_image():    
    # load baseline table
    baselineFile = np.genfromtxt('raw/baseline_table.dat',dtype=str)
    time = baselineFile[:,1].astype(float)
    baseline = baselineFile[:,4].astype(float)
 
    #GMTSAR (currently) guarantees that this file has the same order of lines as baseline_table.dat.
    dataDotIn=np.genfromtxt('raw/data.in',dtype='str').tolist()
    
    # calculate shortest distance from median to scenes
    consider_time=True
    if consider_time:
        time_baseline_scale=1 #arbitrary scaling factor, units of (meters/day)
        sceneDistance = np.sqrt(((time-np.median(time))/time_baseline_scale)**2 + (baseline-np.median(baseline))**2)
    else:
        sceneDistance = np.sqrt((baseline-np.median(baseline))**2)
    
    minID=np.argmin(sceneDistance)
    masterID=dataDotIn[minID]    
   
    # put masterId in the first line of data.in
    dataDotIn.pop(dataDotIn.index(masterID))
    dataDotIn.insert(0,masterID)
    
    os.rename('raw/data.in','raw/data.in.old')
    np.savetxt('raw/data.in',dataDotIn,fmt='%s')
    return masterID

