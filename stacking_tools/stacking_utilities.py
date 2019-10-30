# Sentinel Utilities

import subprocess
import os
import sys
import glob
import datetime as dt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
import numpy as np
import collections
import netcdf_read_write



def write_super_master_batch_config(masterid):
    ifile=open('batch.config','r');
    ofile=open('batch.config.new','w');
    for line in ifile:
        if 'master_image' in line:
            ofile.write('master_image = '+masterid+'\n');
        else:        
            ofile.write(line);
    ifile.close();
    ofile.close();
    subprocess.call(['mv','batch.config.new','batch.config'],shell=False);
    print("Writing master_image into batch.config");
    return;


def get_list_of_intfs(config_params):
    # This is where some hand-picking takes place
    select_intf_list=[];
    total_intf_list=glob.glob("F"+config_params.swath+"/intf_all/???????_???????/unwrap.grd");

    select_criterion=3; # 3 years, 2 years, 1 year

    for item in total_intf_list:
        dates = item.split("/")[-2];
        year1 = dates[0:4];
        year2 = dates[8:12];
        day1  = str(int(dates[4:7])+1);
        day2  = str(int(dates[12:15])+1);
        date1 = dt.datetime.strptime(year1+day1,"%Y%j")
        date2 = dt.datetime.strptime(year2+day2,"%Y%j")
        deltat = date2-date1
        if deltat.days > select_criterion*0.9*365: # a year plus or minus a month
            select_intf_list.append(item);

    print("Out of %d possible interferograms, we are using %d" % (len(total_intf_list), len(select_intf_list)) );
    return select_intf_list;


def make_referenced_unwrapped(intf_list, swath, ref_swath, rowref, colref, ref_dir):
    output_dir="F"+swath+"/"+ref_dir;
    print("Imposing reference pixel on %d files; saving output in %s" % (len(intf_list), output_dir) );
    if ref_swath != swath:
        print("Issue we haven't touched yet: reference pixel is not within this swath. ");
    else:
        saving_file = open(output_dir+"/reference_pixels_orig.txt",'w');
        for item in intf_list:
            [xdata,ydata,zdata] = netcdf_read_write.read_grd_xyz(item);
            individual_name=item.split('/')[-1];
            intf_name=item.split('/')[-2];  # ex: 2015178_2018180
            referenced_zdata=np.zeros(np.shape(zdata));
            for i in range(len(ydata)):
                for j in range(len(xdata)):
                    referenced_zdata[i][j]=zdata[i][j]-zdata[rowref][colref];
            outname=output_dir+"/"+intf_name+"_"+individual_name;
            netcdf_read_write.produce_output_netcdf(xdata, ydata, referenced_zdata, 'phase', outname);
            saving_file.write("%s %d %d %f %f\n" % (item, rowref, colref, zdata[rowref][colref], referenced_zdata[rowref][colref]) )
        saving_file.close();
    return;


def get_ref_index(ref_swath, swath, ref_loc, ref_idx):
    if swath != ref_swath:
        rowref, colref = -1, -1;  # this is a separate case. 
    else:
        if ref_idx != []:  # if we already have an index location... 
            rowref=int(ref_idx.split('/')[0])
            colref=int(ref_idx.split('/')[1])
        else:
            rowref=0; colref=0;
            # Here we will run ll2ra in the future. 
    return rowref, colref;




