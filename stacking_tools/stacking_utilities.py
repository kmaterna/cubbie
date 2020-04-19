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
        total_intf_list=glob.glob(config_params.intf_dir+"/???????_???????/unwrap.grd");
    elif config_params.SAT=="UAVSAR":
        # Specific to the case of UAVSAR stacks with alt-unwrapped taking place
        total_intf_list=glob.glob("../Igrams/*/alt_unwrapped/filt*_fully_processed.uwrappedphase");
    print("Identifying all unwrapped intfs in %s: " % config_params.intf_dir);
    print("Found %d interferograms for stacking. " % (len(total_intf_list)) );
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

def get_ref_index_merged(ref_loc, ref_idx, intf_files):
    # Using merged files, get the index of the reference pixel. 
    # If you don't have the reference pixel in the config file, 
    # the program will stop execution so you can write it there.
    print("Identifying reference pixel:"); 
    if ref_idx != "":
        rowref=int(ref_idx.split('/')[0])
        colref=int(ref_idx.split('/')[1])
        print("Returning rowref, colref %d, %d from config file" % (rowref, colref) );
    else:
        lon=float(ref_loc.split('/')[0])
        lat=float(ref_loc.split('/')[1])        
        print("Finding coordinate %.3f, %.3f in grd file" % (lon, lat) )
        trans_dat="merged/trans.dat";
        example_grd=intf_files[0];
        rowref, colref = get_index_merged(lon, lat, trans_dat, example_grd);            
        print("Found reference coordinate at row, col %d, %d " % (rowref, colref) );
        print("Please write this into your config file. ")
        sys.exit(0);
    return rowref, colref;

def get_index_merged(lon, lat, trans_dat, example_grd):
    print("Finding coordinate %.3f, %.3f in grd file" % (lon, lat) )
    [ra, az] = get_ra_rc_from_ll.get_ra_from_ll(trans_dat, example_grd, lon, lat);
    if np.isnan(ra) or np.isnan(az):
        print("WARNING: Cannot Find %f %f in file." % (lon, lat) );
        rowref = np.nan; colref = np.nan; 
    else:
        [rowref, colref] = get_ra_rc_from_ll.get_nearest_row_col(example_grd, ra, az); 
        print("Found Coordinate at row, col %d, %d " % (rowref, colref) );
    return rowref, colref;

def exclude_intfs_manually(total_intf_list, skip_file):
    print("Excluding intfs based on manual_exclude file %s." % skip_file);
    print(" Started with %d total interferograms. " % (len(total_intf_list)) );
    if skip_file=="":
        print(" No manual exclude file provided.\n Returning all %d interferograms. " % (len(total_intf_list)) );
        selected_intf_list=total_intf_list;
    else:
        print(" Excluding the following interferograms based on SkipFile %s: " % skip_file);
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
        print(" Returning %d interferograms " % len(selected_intf_list) );
    return selected_intf_list;

def get_intf_dates_gmtsar(total_intf_list):
    d1 = []; d2 = [];
    for item in total_intf_list:
        datesplit = item.split('/')[-1];  # example: 2015157_2018177_unwrap.grd
        date1 = str(int(datesplit[0:7]) + 1);
        date2 = str(int(datesplit[8:15]) + 1);  # adding 1 to the date because 000 = January 1
        d1.append(dt.datetime.strptime(date1,"%Y%j"));
        d2.append(dt.datetime.strptime(date2,"%Y%j"));
    return d1, d2;

def get_intf_dates_gmtsar_merged(total_intf_list):
    d1 = []; d2 = [];
    for item in total_intf_list:
        datesplit = item.split('/')[-2];  # example: merged/2015133_2015157/unwrap.grd
        date1 = str(int(datesplit[0:7]) + 1);
        date2 = str(int(datesplit[8:15]) + 1);  # adding 1 to the date because 000 = January 1
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

def include_intfs_by_time_range(intf_list, d1, d2, start_time, end_time, coseismic):
    # Here, we look for each interferogram that falls totally within the time range 
    # given in the config file.
    # It also implements a filter for spanning a coseismic interval, if you include one. 
    print("Including only interferograms in time range %s to %s ." % (dt.datetime.strftime(start_time,"%Y-%m-%d"),
        dt.datetime.strftime(end_time, "%Y-%m-%d"))); 
    print(" Starting with %d interferograms " % len(intf_list) )
    select_intf_list=[];
    for i in range(len(intf_list)):
        if d1[i]>=start_time and d1[i]<=end_time:
            if d2[i]>=start_time and d2[i]<=end_time:
                if coseismic=="":
                    select_intf_list.append(intf_list[i]);  # in the case of no coseismic constraint
                else:
                    if d1[i]<coseismic and d2[i]>coseismic:
                        select_intf_list.append(intf_list[i]);  # in the case of a coseismic constraint
    print(" Returning %d interferograms " % len(select_intf_list));
    return select_intf_list;

def make_selection_of_intfs(config_params):
    # For TS and velocities:  
    # Get all ref_unwrapped
    if config_params.SAT=="S1":
        total_intf_list=glob.glob(config_params.intf_dir+"/*/unwrap.grd");  # the GMTSAR workflow
        d1, d2 = get_intf_dates_gmtsar_merged(total_intf_list);
    elif config_params.SAT=="UAVSAR":
        total_intf_list=glob.glob(config_params.ref_dir+"/????????_????????.refunwrapped");  # The ISCE workflow (fix later)
        d1, d2 = get_intf_dates_isce(total_intf_list);

    # Use the config file to excluse certain time ranges and implement coseismic constraints
    select_intf_list = include_intfs_by_time_range(total_intf_list, d1, d2, 
        config_params.start_time, config_params.end_time, config_params.coseismic);

    # Employing the Manual Removes
    select_intf_list = exclude_intfs_manually(select_intf_list, config_params.skip_file);

    # ------------------------------ # 
    # HERE IS WHERE YOU SELECT WHICH INTERFEROGRAMS YOU WILL BE USING.
    # THIS GETS CREATIVE.  
    # IN A GENERAL CASE, WE WILL NOT BE SELECTING ONLY LONG INTERFEROGRAMS
    # WE MIGHT APPLY A MANUAL EXCLUDE, OR A TIME CONSTRAINT. 
    # THIS DEPENDS ON YOUR CONFIG SETTINGS
    # I THINK WE MIGHT WANT TO SELECT ALL INTERFEROGRAMS
    # select_criterion=0.8; # 3+ years, 2+ years, 1+ year
    # THIS IS WHERE YOU WILL MAKE CHANGES
    # ------------------------------ # 
    if config_params.ts_type=="STACK":
        # Then we are stacking only the long interferograms. 
        print("Selecting interferograms for long-intf stacking and velocity formation.");
        long_intfs= [];
        d1, d2 = get_intf_dates_gmtsar_merged(select_intf_list);
        for i in range(len(select_intf_list)):
            delta = d2[i]-d1[i];
            if delta.days > 300:
                long_intfs.append(select_intf_list[i]);
        select_intf_list=long_intfs;

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
    # Selecting coherence information for each intf_file we've already determined. 
    coh_files=[];
    for i in range(len(intf_files)):
        if config_params.SAT=="UAVSAR":
            datestring = intf_files[i].split('/')[-1][0:17];
            coh_file = "../Igrams/"+datestring+"/alt_unwrapped/filt_"+datestring+"_cut.cor";
            coh_files.append(coh_file);
        if config_params.SAT=="S1":
            coh_file = intf_files[i].replace("unwrap.grd","corr.grd");
            coh_files.append(coh_file);
    print("Returning %d files with coherence information " % (len(coh_files)) )
    return coh_files;


# Functions to get TS points in row/col coordinates
def drive_cache_ts_points(ts_points_file):
    # If you want to re-compute things, you need to delete the cache. 
    if os.path.isfile(ts_points_file+".cache"):
        lons, lats, names, rows, cols = read_ts_points_file(ts_points_file+".cache");
    elif os.path.isfile(ts_points_file): # if there's no cache, we will make one. 
        lons, lats, names, rows, cols = read_ts_points_file(ts_points_file);
        lons, lats, names, rows, cols = match_ts_points_row_col(lons, lats, names, rows, cols);
        write_ts_points_file(lons, lats, names, rows, cols, ts_points_file+".cache");
    else:
        print("Error! You are asking for points but there's not ts_points_file %s . No points computed. " % ts_points_file);
        return [], [], [], [], [];
    return lons, lats, names, rows, cols;

def read_ts_points_file(ts_points_file):
    # Here we can use several formats simultaneously. Point name is required. 
    #Format 1:  -117.76 35.88 313 654 coso1
    #Format 2:  -117.76 35.90 coso2
    print("Reading file %s" % ts_points_file);
    lons=[]; lats=[]; names=[]; rows=[]; cols=[];
    ifile=open(ts_points_file,'r');
    for line in ifile:
        temp=line.split();
        if len(temp)==3:  # we have provided the lat/lon/name
            lons.append(float(temp[0]));
            lats.append(float(temp[1]));
            names.append(temp[2]);
            rows.append('');
            cols.append('');                
        if len(temp)==5:  # we have provided the lat/lon/row/col/name
            if np.isnan(float(temp[2])):  # if the cache has nan for those pixels, we skip. 
                continue;
            else:
                lons.append(float(temp[0]));
                lats.append(float(temp[1]));
                rows.append(int(temp[2]));
                cols.append(int(temp[3]));
                names.append(temp[4]);
    print("Computing time series at %d geographic points " % (len(lons)) );
    return lons, lats, names, rows, cols;

def write_ts_points_file(lons, lats, names, rows, cols, ts_points_file):
    ofile=open(ts_points_file,'w');
    for i in range(len(lons)):
        ofile.write("%.5f %.5f %s %s %s\n" % (lons[i], lats[i], str(rows[i]), str(cols[i]), names[i]) );
    ofile.close();
    return;

def match_ts_points_row_col(lons, lats, names, rows, cols):
    # Find each row and col that hasn't been found before. 
    trans_dat="merged/trans.dat";
    example_grd=glob.glob("merged/*/unwrap.grd")[0];  # specific to the merged case. 
    for i in range(len(lons)):
        if rows[i]=='':
            irow, icol = get_index_merged(lons[i], lats[i], trans_dat, example_grd);
            rows[i]=irow;
            cols[i]=icol;
    return lons, lats, names, rows, cols;



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


