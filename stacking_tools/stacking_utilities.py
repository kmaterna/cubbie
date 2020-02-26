# Sentinel Utilities

import subprocess
import os, sys, glob
import datetime as dt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
import numpy as np
import collections
import netcdf_read_write
import get_ra_rc_from_ll


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

def get_list_of_intf_all(config_params):
    # This is mechanical: just takes the list of interferograms in intf_all. 
    # The smarter selection takes place in make_selection_of_intfs. 
    total_intf_list=glob.glob("F"+config_params.swath+"/intf_all/???????_???????/unwrap.grd");
    return total_intf_list;

def make_referenced_unwrapped(intf_list, swath, ref_swath, rowref, colref, ref_dir):
    # This works for both F1 and F2. You should run whichever swath has the reference point first. 
    # This will break for F3, because we need the F3-F2 offset and the F2-F1 offset. 
    # Writes an output file that shows which interferograms were unsuccessfully referenced. 
    output_dir="F"+swath+"/"+ref_dir;
    errors_file="F"+swath+"/"+ref_dir+"/errors.txt"
    ofile=open(errors_file,'w');
    error_count=0;
    print("Imposing reference pixel on %d files; saving output in %s" % (len(intf_list), output_dir) );

    for item in intf_list:
        # Step 1: get reference offset
        individual_name=item.split('/')[-1];  # ex: unwrap.grd
        intf_name=item.split('/')[-2];  # ex: 2015178_2018180
        F1_name=item.replace('F'+swath,'F1');
        F2_name=item.replace('F'+swath,'F2');
        F3_name=item.replace('F'+swath,'F3');
        is1 = (os.path.isfile(F1_name));
        is2 = (os.path.isfile(F2_name));
        is3 = (os.path.isfile(F3_name)); 

        if swath=='1':
            # This is easy. We have no issues with 2*n*pi
            zvalue = get_reference_value(swath, ref_swath, rowref, colref, item);
        else:
            if is1==1 and is2==1 and swath=='2':  # 2 case
                zvalue_pixel = get_reference_value(ref_swath, ref_swath, rowref, colref, F1_name);  # the reference pixel in F1
                zvalue_2npi = get_n_2pi(F1_name, F2_name);  # this is an integer
                zvalue = zvalue_pixel+zvalue_2npi*-2*np.pi;

            elif is1==1 and is2==1 and is3==1 and swath=='3': # 3 case. 
                print("in swath 3. Will find n pi");
                F2_referenced = "F2/stacking/ref_unwrapped/"+intf_name+"_unwrap.grd"
                zvalue_pixel = get_reference_value(ref_swath, ref_swath, rowref, colref, F1_name);  # the reference pixel in F1
                zvalue_2npi_12 = get_n_2pi(F1_name, F2_name);  # this is an integer
                zvalue_2npi_23 = get_n_2pi(F2_name, F3_name);
                print("n is %d and %d" % (zvalue_2npi_12, zvalue_2npi_23) );
                zvalue = zvalue_pixel+zvalue_2npi_12*-2*np.pi+zvalue_2npi_23*-2*np.pi;
            else: # we don't have enough data to do the referencing.
                print("skipping making ref_unwrapped for %s " % item);
                ofile.write("%s\n" % (item) );
                error_count=error_count+1; 
                continue;

        outname=output_dir+"/"+intf_name+"_"+individual_name;
        print("Making %s " % outname);
        [xdata,ydata,zdata] = netcdf_read_write.read_grd_xyz(item);
        referenced_zdata = apply_reference_value(xdata, ydata, zdata, zvalue);
        netcdf_read_write.produce_output_netcdf(xdata, ydata, referenced_zdata, 'phase', outname); 
    print("Done making reference unwrapped")
    total_intf_list=glob.glob("F"+swath+"/"+ref_dir+"/*unwrap.grd");
    print("%s contains %d files " % (ref_dir, len(total_intf_list)) );
    print("%s files were missed due to referencing errors" % (error_count) );
    ofile.close();
    return;

def get_n_2pi(file1, file2):
    # file1 must be the smaller number (F1)
    print(file1)
    print(file2)
    [xdata_f1,ydata_f1,zdata_f1] = netcdf_read_write.read_grd_xyz(file1);
    [xdata_f2,ydata_f2,zdata_f2] = netcdf_read_write.read_grd_xyz(file2);  
    n = (np.nanmean(zdata_f1[:,-20])-np.nanmean(zdata_f2[:,20]))/(2*np.pi);
    return np.round(n);


def apply_reference_value(xdata, ydata, zdata, reference_value):
    # The math component of applying a reference to a grid. 
    referenced_zdata=np.zeros(np.shape(zdata));
    for i in range(len(ydata)):
        for j in range(len(xdata)):
            referenced_zdata[i][j]=zdata[i][j]-reference_value;
    return referenced_zdata;


def get_reference_value(swath, ref_swath, rowref, colref, item):
    # Reading a file and putting it in the cache, or ignoring this scene due to lack of reference value. 
    if swath==ref_swath:
        [xdata,ydata,zdata] = netcdf_read_write.read_grd_xyz(item);
        reference_value = zdata[rowref][colref];
    else:
        reference_value=np.nan;
        print("Skipping %s because we can't find it in F%s" % (item, str(ref_swath)) );
    return reference_value;


def get_ref_index(ref_swath, swath, ref_loc, ref_idx):
    # THIS IS A DEGENERATE FUNCTION. 
    # IT DOESN'T DO ANYTHING RIGHT NOW. 
    # IT WILL BE CODED LATER. 
    # I have thus far been hard-coding the reference idx in the config file. 
    if ref_idx != []:  # if we already have an index location... 
        rowref=int(ref_idx.split('/')[0])
        colref=int(ref_idx.split('/')[1])
    else:
        rowref=0; colref=0;
        # Here we will run ll2ra in the future. This is a later project. 
    print("Reference pixel: swath = %s, index = [%d, %d] " % (ref_swath, rowref, colref) );
    return rowref, colref;


def test_geocoding(lon, lat):
    # This tests whether the geocoding returns the same value as it started with. 
    try_swaths = ["1","2","3"];
    swath=-1;
    row=-1;
    col=-1;
    for try_swath in try_swaths:
        print("Trying %f %f in swath %s " % (lon, lat, try_swath) );
        trans_dat="F"+try_swath+"/topo/trans.dat";
        example_grd=glob.glob("F"+try_swath+"/stacking/ref_unwrapped/*_unwrap.grd")[0];
        # Assumes there are some unwrapped referenced grd files hanging around to be used. 
        [ra, az] = get_ra_rc_from_ll.get_ra_from_ll(trans_dat, example_grd, lon, lat);
        if np.isnan(ra) or np.isnan(az):
            print("WARNING: Cannot Find %f %f in swath %s." % (lon, lat, try_swath) );
            continue;
        else:
            [row, col] = get_ra_rc_from_ll.get_nearest_row_col(example_grd, ra, az); 
            swath=try_swath;
            print("SUCCESS: Found %f %f in Swath %s, Row %d, Col %d" % (lon, lat, swath, row, col) ); 
            break;

    [lon_return, lat_return] = get_ra_rc_from_ll.get_ll_from_ra(trans_dat, ra, az);  # this is a test. 

    ofile=open("geocode_exp.txt",'a');
    ofile.write("START: This point started as: %f %f \n" % (lon, lat) );
    ofile.write("MIDDLE: This point converted to: %s %d %d \n" % (swath, row, col) );
    ofile.write("FINISH: This point converted back to: %f %f \n\n" % (lon_return, lat_return) );
    ofile.close();
    
    return;


def get_rows_cols(ts_points_file):
    # Get a group of swaths, rows, and columns, one for each TS point
    # In the past, I had tried to make this write a cache to save time, but right now it's not necessary. 
    # If I wanted, I could refactor this to make just three computations (F1, F2, F3) if I really wanted. 
    # Since ra_from_ll works on arrays and returns nans in the right places. 
    lons, lats, names, swaths, rows, cols = read_ts_points_file(ts_points_file);

    if len(lons)==0:
        print("No points requested. Ending without doing any time series calculations.");
        return swaths, rows, cols, names, lons, lats;
    
    for i in range(len(lons)):
        if swaths[i] != '':
            print("Skipping %s: we already have its row/col"  % (names[i]) );

        if swaths[i] == '': # If we don't have the existing swath/row/col of the point, then we compute it. 
            print("Computing swath/row/col for %f %f " % (lons[i], lats[i]) );
            swath, row, col = get_swath_row_col(lons[i], lats[i]);
            swaths[i]=swath;
            rows[i]=row;
            cols[i]=col;

    print("Requesting time series for the following pixels:");
    print(lons)
    print(lats)
    print(names)
    print(swaths)
    print(rows)
    print(cols)
    return lons, lats, names, swaths, rows, cols;


def get_swath_row_col(lon, lat):
    # For a single lat/lon point, what swath and row and column is nearest to it? 
    try_swaths = ["1","2","3"];
    swath=-1;
    row=-1;
    col=-1;
    for try_swath in try_swaths:
        print("Trying %f %f in swath %s " % (lon, lat, try_swath) );
        trans_dat="F"+try_swath+"/topo/trans.dat";
        example_grd=glob.glob("F"+try_swath+"/stacking/ref_unwrapped/*_unwrap.grd")[0];
        # Assumes there are some unwrapped referenced grd files hanging around to be used. 
        [ra, az] = get_ra_rc_from_ll.get_ra_from_ll(trans_dat, example_grd, lon, lat);
        if np.isnan(ra) or np.isnan(az):
            print("WARNING: Cannot Find %f %f in swath %s." % (lon, lat, try_swath) );
            continue;
        else:
            [row, col] = get_ra_rc_from_ll.get_nearest_row_col(example_grd, ra, az); 
            swath=try_swath;
            print("SUCCESS: Found %f %f in Swath %s, Row %d, Col %d" % (lon, lat, swath, row, col) ); 
            break;
    return swath, row, col;




def read_ts_points_file(ts_points_file):
    # Here we can use several formats simultaneously. Point name is not required. 
    #Format 1:  -117.76 35.88 2 313 654 coso1
    #Format 2:  -117.76 35.90 coso2
    #Format 3:  -117.76 35.92
    lons=[]; lats=[]; names=[]; swaths=[]; rows=[]; cols=[];
    if os.path.isfile(ts_points_file):
        ifile=open(ts_points_file,'r');
        for line in ifile:
            temp=line.split();
            if len(temp)==2:  # we have provided the lat/lon
                lons.append(float(temp[0]));
                lats.append(float(temp[1]));
                names.append(temp[0]+'.'+temp[1]);  # giving a name based on the latlon
                swaths.append('');
                rows.append('');
                cols.append('');                
            if len(temp)==3:  # we have provided the lat/lon/name
                lons.append(float(temp[0]));
                lats.append(float(temp[1]));
                names.append(temp[2]);
                swaths.append('');
                rows.append('');
                cols.append('');                
            if len(temp)==6:  # we have provided the lat/lon/swath/row/col/name
                lons.append(float(temp[0]));
                lats.append(float(temp[1]));
                swaths.append(temp[2]);
                rows.append(int(temp[3]));
                cols.append(int(temp[4]));
                names.append(temp[5]);
        print("Computing time series at %d geographic points " % (len(lons)) );
    else:
        print("No ts_points_file %s. Not computing time series at points. " % ts_points_file);    
    return lons, lats, names, swaths, rows, cols;



def replace_line_ts_points_file(ts_points_file, name, swath, row, col):
    # Careful: this function re-writes a certain line of the file. 
    ifile=open(ts_points_file,'r');
    ofile=open("temp.txt",'w');
    for line in ifile:
        if ' '+name in line:
            temp=line.split();
            ofile.write("%s %s %s %s %s %s\n" % (temp[0], temp[1], temp[2], swath, row, col) );
        else:
            ofile.write(line);
    ifile.close();
    ofile.close();
    subprocess.call(['mv','temp.txt',ts_points_file],shell=False);
    return;




def make_selection_of_intfs(config_params, swath=0):
    # If you want to pass in the swath, then use an explicit argument. 
    # Otherwise we'll use the swath in the config file. 
    swath=str(swath);
    if swath=='0':
        swath=config_params.swath;

    total_intf_list=glob.glob("F"+swath+"/"+config_params.ref_dir+"/???????_???????_unwrap.grd");

    # ------------------------------ # 
    # HERE IS WHERE YOU SELECT WHICH INTERFEROGRAMS YOU WILL BE USING. 
    # IN A GENERAL CASE, WE WILL NOT BE SELECTING ONLY LONG INTERFEROGRAMS
    # WE MIGHT APPLY A MANUAL EXCLUDE, OR A TIME CONSTRAINT. 
    # THIS DEPENDS ON YOUR CONFIG SETTINGS
    # I THINK WE MIGHT WANT TO SELECT ALL INTERFEROGRAMS
    # FEB 2020
    # select_criterion=0.8; # 3+ years, 2+ years, 1+ year

    # for item in total_intf_list:
    #     dates = item.split("/")[-2];
    #     year1 = dates[0:4];
    #     year2 = dates[8:12];
    #     day1  = str(int(dates[4:7])+1);
    #     day2  = str(int(dates[12:15])+1);
    #     date1 = dt.datetime.strptime(year1+day1,"%Y%j")
    #     date2 = dt.datetime.strptime(year2+day2,"%Y%j")
    #     deltat = date2-date1
    #     if deltat.days > select_criterion*0.9*365: # a year plus or minus a month
    #         select_intf_list.append(item);  
    #  print("Out of %d possible interferograms, we are trying to use %d" % (len(total_intf_list), len(select_intf_list)) );

    # THIS IS WHERE YOU WILL MAKE CHANGES
    # ------------------------------ # 

    select_intf_list=total_intf_list;

    # Writing the exact interferograms used in this run. 
    record_file="F"+swath+"/"+config_params.ts_output_dir+"/"+"intf_record.txt";
    ofile=open(record_file,'w');
    ofile.write("List of %d interferograms used in this run:\n" % (len(select_intf_list)) );
    for item in select_intf_list:
        ofile.write("%s\n" % (item) );
    ofile.close();

    return select_intf_list; 


