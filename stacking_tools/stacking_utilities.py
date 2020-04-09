# Sentinel Utilities

import subprocess
import os, sys, glob
import datetime as dt
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
import numpy as np
import re 
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
    if config_params.SAT=="S1":
        total_intf_list=glob.glob("F"+config_params.swath+"/intf_all/???????_???????/unwrap.grd");
    elif config_params.SAT=="UAVSAR":
        # Specific to the case of UAVSAR stacks with alt-unwrapped taking place
        total_intf_list=glob.glob("../Igrams/*/alt_unwrapped/filt*_fully_processed.uwrappedphase");

    return total_intf_list;

def get_xdates_from_intf_tuple(intf_tuple):
    total_dates=[];
    for item in intf_tuple.dates_correct:
        date1=dt.datetime.strptime(item[0:7],"%Y%j");
        date2=dt.datetime.strptime(item[8:15],"%Y%j");
        if date1 not in total_dates:
            total_dates.append(date1);
        if date2 not in total_dates:
            total_dates.append(date2);
    xdates = sorted(total_dates);
    return xdates;

def make_referenced_unwrapped(intf_list, swath, ref_swath, rowref, colref, output_dir):
    # This is a pretty specific function for Sentinel stacks. Most of the logic relates to swaths. 
    # This works for F1, F2, and F3. You should run whichever swath has the reference point first. 
    # Writes an output file that shows which interferograms were unsuccessfully referenced. 
    errors_file=output_dir+"/errors.txt"
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


def get_ref_index(ref_swath, swath, ref_loc, ref_idx, intf_files):
    # For sentinel, we usually put some kind of lat/lon as the reference point
    # Of we have already identified a ref_idx. 
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
    # GMTSAR version
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


def exclude_intfs_manually(total_intf_list, skip_file):
    print("Excluding intfs based on manual_exclude.");
    print("Started with %d total interferograms. " % (len(total_intf_list)) );
    print("Excluding the following interferograms based on SkipFile %s: " % skip_file);
    ifile=open(skip_file,'r');
    manual_removes = [];
    for line in ifile:
        manual_removes.append(line.split()[0]);
    ifile.close();
    print(manual_removes);

    if manual_removes==[]:
        selected_intf_list = total_intf_list;
    else:
        # Checking to see if each interferogram should be included. 
        selected_intf_list=[];
        for igram in total_intf_list:
            include_flag = 1;
            for scene in manual_removes:
                if scene in igram:
                    include_flag=0;
            if include_flag==1:
                selected_intf_list.append(igram);
    print("Returning %d interferograms " % len(selected_intf_list) );
    return selected_intf_list;

def get_intf_dates_gmtsar(total_intf_list):
    d1 = []; d2 = [];
    for item in total_intf_list:
        datesplit = item.split('/')[-1];  # example: 2015157_2018177_unwrap.grd
        date1 = str(int(datesplit[0:7]) + 1);
        date2 = str(int(date_new[8:15]) + 1);  # adding 1 to the date because 000 = January 1
        d1.append(dt.datetime.strptime(date1,"%Y%j"));
        d2.append(dt.datetime.strptime(date2,"%Y%j"));
    return d1, d2;

def get_intf_dates_isce(total_intf_list):
    d1 = []; d2 = [];
    for item in total_intf_list:
        datesplit = re.findall(r"\d\d\d\d\d\d\d\d_\d\d\d\d\d\d\d\d", item)[0]; #  example: 20100402_20140304
        date1 = dt.datetime.strptime(datesplit[0:8],"%Y%m%d");
        date2 = dt.datetime.strptime(datesplit[9:17],"%Y%m%d");
        d1.append(date1);
        d2.append(date2);
    return d1, d2;


def include_intfs_by_time_range(intf_list, d1, d2, start_time, end_time):
    # Here, we look for each interferogram that falls totally within the time range 
    # given in the config file.
    print("Including only interferograms in time range %s to %s ." % (dt.datetime.strftime(start_time,"%Y-%m-%d"),
        dt.datetime.strftime(end_time, "%Y-%m-%d"))); 
    print("Starting with %d interferograms " % len(intf_list) )
    select_intf_list=[];
    for i in range(len(intf_list)):
        if d1[i]>=start_time and d1[i]<=end_time:
            if d2[i]>=start_time and d2[i]<=end_time:
                select_intf_list.append(intf_list[i]);
    print("Returning %d interferograms " % len(select_intf_list));
    return select_intf_list;



def make_selection_of_intfs(config_params, swath=0):

    # Get all ref_unwrapped
    if config_params.SAT=="S1":
        total_intf_list=glob.glob(config_params.ref_dir+"/???????_???????_unwrap.grd");  # the GMTSAR workflow
        d1, d2 = get_intf_dates_gmtsar(total_intf_list);
    elif config_params.SAT=="UAVSAR":
        total_intf_list=glob.glob(config_params.ref_dir+"/????????_????????.refunwrapped");  # The ISCE workflow
        d1, d2 = get_intf_dates_isce(total_intf_list);

    # Use the config file to excluse certain time ranges
    select_intf_list = include_intfs_by_time_range(total_intf_list, d1, d2, config_params.start_time, config_params.end_time);

    # Use the config file to impose MUST COVER COSEISMIC
    # This would be for making an average coseismic interferogram however. 
    # Slightly different workflow. 

    # Employing the Manual Removes
    select_intf_list = exclude_intfs_manually(select_intf_list, config_params.skip_file);


    # ------------------------------ # 
    # HERE IS WHERE YOU SELECT WHICH INTERFEROGRAMS YOU WILL BE USING. 
    # IN A GENERAL CASE, WE WILL NOT BE SELECTING ONLY LONG INTERFEROGRAMS
    # WE MIGHT APPLY A MANUAL EXCLUDE, OR A TIME CONSTRAINT. 
    # THIS DEPENDS ON YOUR CONFIG SETTINGS
    # I THINK WE MIGHT WANT TO SELECT ALL INTERFEROGRAMS
    # FEB 2020
    # select_criterion=0.8; # 3+ years, 2+ years, 1+ year
    # THIS IS WHERE YOU WILL MAKE CHANGES
    # ------------------------------ # 


    # Writing the exact interferograms used in this run. 
    record_file=config_params.ts_output_dir+"/"+"intf_record.txt";
    print("Writing out list of %d interferograms used in this run to %s" % (len(select_intf_list), record_file) );
    ofile=open(record_file,'w');
    ofile.write("List of %d interferograms used in this run:\n" % (len(select_intf_list)) );
    for item in select_intf_list:
        ofile.write("%s\n" % (item) );
    ofile.close();

    return select_intf_list; 

def make_selection_of_coh_files(config_params, intf_files):
    coh_files=[];
    for i in range(len(intf_files)):
        if config_params.SAT=="UAVSAR":
            datestring = intf_files[i].split('/')[-1][0:17];
            coh_file = "../Igrams/"+datestring+"/alt_unwrapped/filt_"+datestring+"_cut.cor";
            coh_files.append(coh_file);
    print("Returning %d files with coherence information " % (len(coh_files)) )
    return coh_files;


def get_axarr_numbers(rows, cols, idx):
    # Given an incrementally counting idx number and a subplot dimension, where is our plot? 
    total_plots = rows*cols;
    col_num = np.mod(idx, cols);
    row_num = int(np.floor(idx/cols));
    return row_num, col_num;


def plot_full_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1):
    # Make a nice time series plot. 
    tdata, xdata, ydata, TS_array = netcdf_read_write.read_3D_netcdf(TS_NC_file);
    num_rows_plots=3;
    num_cols_plots=4;

    f, axarr = plt.subplots(num_rows_plots,num_cols_plots,figsize=(16,10),dpi=300);
    for i in range(len(xdates)):
        rownum, colnum = get_axarr_numbers(num_rows_plots,num_cols_plots,i);
        axarr[rownum][colnum].imshow(TS_array[i,:,:],aspect=aspect,cmap='rainbow',vmin=vmin,vmax=vmax);
        titlestr = dt.datetime.strftime(xdates[i],"%Y-%m-%d");
        axarr[rownum][colnum].get_xaxis().set_visible(False);
        axarr[rownum][colnum].set_title(titlestr,fontsize=20);

    cbarax = f.add_axes([0.75,0.35,0.2,0.3],visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(TS_image_file);
    return;



def plot_incremental_timeseries(TS_NC_file, xdates, TS_image_file, vmin=-50, vmax=200, aspect=1):
    # Make a nice time series plot. 
    # With incremental displacement data. 
    tdata, xdata, ydata, TS_array = netcdf_read_write.read_3D_netcdf(TS_NC_file);
    num_rows_plots=3;
    num_cols_plots=4;

    # Combining the two shortest intervals into one. 
    print(np.shape(TS_array));
    selected = [0,1,3,4,5,6,7,8,9,10];
    TS_array = TS_array[selected,:,:];
    xdates = [xdates[i] for i in range(11) if i in selected];
    print(np.shape(TS_array));
    print(len(xdates));

    f, axarr = plt.subplots(num_rows_plots,num_cols_plots,figsize=(16,10),dpi=300);
    for i in range(1,len(xdates)):
        rownum, colnum = get_axarr_numbers(num_rows_plots,num_cols_plots,i);
        data = np.subtract(TS_array[i,:,:],TS_array[i-1,:,:]);
        axarr[rownum][colnum].imshow(data,aspect=aspect,cmap='rainbow',vmin=vmin,vmax=vmax);
        titlestr = dt.datetime.strftime(xdates[i],"%Y-%m-%d");
        axarr[rownum][colnum].get_xaxis().set_visible(False);
        axarr[rownum][colnum].set_title(titlestr,fontsize=20);

    cbarax = f.add_axes([0.75,0.35,0.2,0.3],visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap,aspect=12,fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(TS_image_file);
    return;


