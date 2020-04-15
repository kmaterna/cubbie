import collections
import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
from subprocess import call, check_output
import stacking_utilities
import aps 
import detrend_atm_topo
import flattentopo_driver
import phasefilt_plot
import sbas
import nsbas
import nsbas_accessing
import gps_into_LOS
import Super_Simple_Stack as sss
import netcdf_read_write as rwr
import stack_corr


# --------------- STEP 0: Setting up ------------ # 
def set_up_output_directories(config_params):
    if config_params.startstage>0:  
        return;
    if config_params.endstage<0:   
        return;
    print("Stage 0 - Setting up output directories within %s." % config_params.ts_parent_dir);
    call(['mkdir','-p',config_params.ts_parent_dir],shell=False);
    # call(['mkdir','-p',config_params.ref_dir],shell=False);
    call(['mkdir','-p',config_params.ts_output_dir],shell=False);
    call(['cp','stacking.config',config_params.ts_output_dir],shell=False);
    call(['cp',config_params.skip_file, config_params.ts_output_dir],shell=False);
    print("Setting up output directory at %s " % config_params.ts_output_dir);
    return;


# --------------- STEP 1: Make corrections ------------ # 
def make_corrections(config_params):
    if config_params.startstage>1:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<1:   # if we're ending at intf, we don't do this. 
        return;  
    print("Stage 1 - Doing optional atm corrections");  
    # This is where we would implement GACOS, detrending, other atmospheric corrections, or unwrapping errors. 
    # Step 3A: Solve or exclude unwrapping errors
    # Step 3B: Try APS-based atmospheric correction
    # Step 3C: Detrend topo-correlated atmosphere   
    # if config_params.gacos:
    #     #gacos.main_function(prior_staging_directory, post_staging_directory, rowref, colref);
    # if config_params.aps:
    #     # aps.main_function(prior_staging_directory, post_staging_directory, rowref, colref,'');
    # if config_params.detrend_atm_topo:
    #     detrend_atm_topo.main_function(prior_staging_directory, post_staging_directory, rowref, colref);
    return;


# --------------- STEP 2: Make ref_unwrap.grd ------------ # 

def collect_unwrap_ref(config_params):
    if config_params.startstage>2:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<2:   # if we're ending at intf, we don't do this. 
        return;

    print("Stage 2 - Finding Reference Information.");

    # Very general, takes all files and doesn't discriminate. 
    intf_files=stacking_utilities.get_list_of_intf_all(config_params);

    # Here we need to get ref_idx if we don't have it already
    rowref, colref = stacking_utilities.get_ref_index_merged(config_params.ref_loc, config_params.ref_idx, intf_files); 

    # Now we coalesce the files and reference them to the right value/pixel
    # stacking_utilities.make_referenced_unwrapped(intf_files, config_params.swath, config_params.ref_swath, rowref, colref, config_params.ref_dir);
    return;


# --------------- STEP 3: Velocities and Time Series! ------------ # 
def vels_and_ts(config_params):
    if config_params.startstage>3:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<3:   # if we're ending at intf, we don't do this. 
        return;

    # This is where the hand-picking takes place: manual excludes, long intfs only, ramp-removed, atm-removed, etc.
    intfs = stacking_utilities.make_selection_of_intfs(config_params);
    
    # Plumbing stuff
    rowref=int(config_params.ref_idx.split('/')[0]);
    colref=int(config_params.ref_idx.split('/')[1]);
    call(['cp','stacking.config',config_params.ts_output_dir],shell=False);

    # Make signal_spread here. Can be commented if you already have it. 
    # corr_files = [i.replace("unwrap.grd","corr.grd") for i in intfs];
    # stack_corr.drive_signal_spread_calculation(corr_files, 0.1, config_params.ts_output_dir);
 
    if config_params.ts_type=="STACK":
        print("Running velocities by simple stack.")
        # sss.drive_velocity_simple_stack(intfs, config_params.wavelength, config_params.ts_output_dir);
    if config_params.ts_type=="SBAS":
        print("Running velocities and time series by SBAS");
        # sbas.drive_velocity_sbas(intfs, config_params.sbas_smoothing, config_params.wavelength, config_params.ts_output_dir);
    if config_params.ts_type=="NSBAS":
        print("Running velocities and time series by NSBAS");
        # nsbas_accessing.drive_velocity_gmtsar(intfs, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir);
        nsbas_accessing.drive_point_ts_gmtsar(intfs, config_params.ts_points_file, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir);
        # nsbas_accessing.drive_full_TS_gmtsar(intfs, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir); 
    if config_params.ts_type=="WNSBAS":
        print("Running velocities and time series by WNSBAS");
        # coh_files = stacking_utilities.make_selection_of_coh_files(config_params, intfs);
        # nsbas_accessing.drive_velocity_gmtsar(intfs, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir, coh_files=coh_files);
        # nsbas_accessing.drive_point_ts_gmtsar(intfs, config_params.ts_points_file, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir, coh_files=coh_files);
        # nsbas_accessing.drive_full_TS_gmtsar(intfs, config_params.nsbas_min_intfs, config_params.sbas_smoothing, config_params.wavelength, rowref, colref, config_params.ts_output_dir, coh_files=coh_files); 
    return; 


# --------------- STEP 4: Geocoding Velocities ------------ # 
def geocode_vels(config_params):
    if config_params.startstage>4:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<4:   # if we're ending at intf, we don't do this. 
        return; 

    directory = config_params.ts_output_dir
    vel_name = "velo_nsbas"

    outfile=open("geocoding.txt",'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to geocode velocities.\n\n");
    outfile.write("cd F"+str(config_params.swath)+"\n");
    outfile.write("geocode_mod.csh "+vel_name+".grd "+vel_name+"_ll.grd "+vel_name+"_ll "+directory+"\n");
    outfile.close();
    print("Ready to call geocoding.txt.")
    call("chmod +x geocoding.txt",shell=True);
    call("./geocoding.txt",shell=True); 
    
    return; 


    # For later plotting, we want to project available GPS into LOS. 
    # gps_into_LOS.top_level_driver(config_params, rowref, colref);




