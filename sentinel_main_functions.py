import collections
import os,sys,shutil,argparse,time,configparser, glob
import numpy as np
from subprocess import call, check_output
import sentinel_utilities
import analyze_coherence
import choose_reference_pixel
import unwrapping_errors
import aps 
import detrend_atm_topo
import rose_baseline_plot
import flattentopo_driver
import phasefilt_plot
import gps_into_LOS

Params=collections.namedtuple('Params',['config_file','SAT','wavelength','startstage','endstage','master','intf_file','orbit_dir','tbaseline','xbaseline','restart',
    'mode','swath','polarization','frame1','frame2','numproc','ts_type','bypass','sbas_smoothing','nsbas_min_intfs','choose_refpixel',
    'solve_unwrap_errors','detrend_atm_topo','gacos','aps','start_time','end_time','threshold_snaphu','gps_file','flight_angle','look_angle','skip_file','ts_output_dir']);

def read_config():
    ################################################
    # Stage 0: Read and check config parameters
    #
    # read command line arguments and parse config file.
    parser = argparse.ArgumentParser(description='Run GMTSAR batch processing. Default automatically determines master, does alignment, and runs all possible interferograms.')
    parser.add_argument('config',type=str,help='supply name of config file to setup processing options. Required.')
    parser.add_argument('--mpi',action='store_true',help='Use MPI (default: false, uses python multiprocessing library instead).')
    parser.add_argument('--debug',action='store_true',help='Print extra debugging messages (default: false)')
    args = parser.parse_args()

    # read config file
    config=configparser.ConfigParser()
    config.optionxform = str #make the config file case-sensitive
    config.read(args.config)
    
    # Setup gnu parallel multiprocessing tool
    numproc=config.getint('py-config','num_processors')

    # get options from config file
    config_file_orig=args.config;
    SAT=config.get('py-config','satellite')
    wavelength=config.getfloat('py-config','wavelength')
    startstage=config.getint('py-config','startstage')
    endstage=config.getint('py-config','endstage')    
    master=config.get('csh-config','master_image')
    intf_file=config.get('py-config','intf_file')
    orbit_dir=config.get('py-config','orbit_dir')
    tbaseline=config.getint('py-config','max_timespan')
    xbaseline=config.getint('py-config','max_baseline') 
    restart=config.getboolean('py-config','restart')
    mode=config.get('py-config','mode') 
    swath=config.get('py-config','swath') 
    polarization=config.get('py-config','polarization') 
    frame_nearrange1=config.get('py-config','frame_nearrange1')
    frame_nearrange2=config.get('py-config','frame_nearrange2')
    ts_type=config.get('timeseries-config','ts_type')
    choose_refpixel = config.getint('timeseries-config','choose_refpixel');
    solve_unwrap_errors = config.getint('timeseries-config','solve_unwrap_errors');
    gacos = config.getint('timeseries-config','gacos');
    aps = config.getint('timeseries-config','aps');
    detrend_atm_topo = config.getint('timeseries-config','detrend_atm_topo');
    bypass=config.get('timeseries-config','bypass')
    nsbas_min_intfs=config.getint('timeseries-config','nsbas_min_intfs');
    sbas_smoothing = config.getfloat('timeseries-config','sbas_smoothing');
    start_time = config.getint('timeseries-config','start_time');
    end_time = config.getint('timeseries-config','end_time');
    gps_file = config.get('timeseries-config','gps_file');
    flight_angle = config.getfloat('timeseries-config','flight_angle');
    look_angle = config.getfloat('timeseries-config','look_angle');
    skip_file = config.get('timeseries-config','skip_file');
    ts_output_dir = config.get('timeseries-config','ts_output_dir');
    threshold_snaphu=config.getfloat('csh-config','threshold_snaphu');

    
    # print config options
    if args.debug:
        print('Running gmtsar_app.py:')
        if args.mpi:
            print('  Using MPI')
        else:
            print('  Not using MPI')
        print('  config file:',args.config,'contains the following options')    
        print(config.write(sys.stdout))

    print("Running sentinel batch processing, starting with stage %d" % startstage);
    
    # check config options
    if intf_file == '':
        intf_file = 'intf_batch.in'
    
    # only valid option for mode is 'scan'    
    if mode != '' and mode != 'scan':
        print('warning: invalid option mode = %s in config file, using blank'%mode)
        mode = ''

    # if master specified in the config file disagrees with existing data.in, we must re-do the pre-processing.
    if master and startstage > 1 and os.path.isfile('F'+str(swath)+'/raw/data.in'):
        #check master in data.in
        dataDotIn=np.genfromtxt('F'+str(swath)+'/raw/data.in',dtype='str')
        oldmaster=dataDotIn[0]
        if '-'+master[3:11]+'t' not in oldmaster:  # For sentinel, oldmaster is formatted like s1a-iw1-slc-vv-20171201t142317-...; master is formatted like S1A20171213_ALL_F1
            print('Warning: The master specified in the config file disagrees with the old master in data.in.')
            #print('We will re-run starting from pre-processing with the master from the config file.')
            #startstage = 1
            #restart = True
    
    # if data.in is not found, we must do pre-processing.
    if startstage > 1 and not os.path.isfile('F'+str(swath)+'/raw/data.in'):
        print('Warning: Pre-processing has not been run, changing startstage to 1')
        startstage = 1
            
    # enforce startstage <= endstage
    if endstage < startstage:
        print('Warning: endstage is less than startstage. Setting endstage = startstage.')
        endstage = startstage

    #logtime is the timestamp added to all logfiles created during this run
    logtime=time.strftime("%Y_%m_%d-%H_%M_%S")
    config_file='batch.run.'+ logtime +'.cfg'
    with open(config_file, 'w') as configfilehandle:
        config.write(configfilehandle)


    config_params=Params(config_file=config_file_orig, SAT=SAT,wavelength=wavelength,startstage=startstage,endstage=endstage,master=master,intf_file=intf_file,
        orbit_dir=orbit_dir,tbaseline=tbaseline, xbaseline=xbaseline,restart=restart,mode=mode,swath=swath,polarization=polarization,frame1=frame_nearrange1, frame2=frame_nearrange2, 
        numproc=numproc, ts_type=ts_type, bypass=bypass, sbas_smoothing=sbas_smoothing, nsbas_min_intfs=nsbas_min_intfs, choose_refpixel=choose_refpixel, 
        solve_unwrap_errors=solve_unwrap_errors, gacos=gacos, aps=aps, detrend_atm_topo=detrend_atm_topo, start_time=start_time, end_time=end_time, 
        threshold_snaphu=threshold_snaphu, gps_file=gps_file, flight_angle=flight_angle, look_angle=look_angle, skip_file=skip_file, ts_output_dir=ts_output_dir);

    return config_params; 


# --------------- STEP -1: Making frames from bursts (will skip if FRAMES contains data) ------------ # 
def compile_frame_from_bursts(config_params):
    if config_params.startstage>-1:  # don't need to set up if we're starting mid-stream. 
        return;
    if config_params.endstage<-1: # don't need to do this step
        return; 
    file_list = stitch_frames(config_params);  # will assemble frames if necessary. Otherwise will just return the DATA/*.SAFE files        
    return;


def stitch_frames(config_params):
    # This will read the batch.config file, make frames, and put the results into the file_list 
    # (for later porting into raw_orig)
    if config_params.frame1 != '':
        # Here we want a frame to be made. 
        call(['mkdir','-p','FRAMES'],shell=False);
        frame1_def = config_params.frame1.split('/');
        frame2_def = config_params.frame2.split('/');

        # # write the data list to data.list
        outfile=open("make_frame_commands.sh",'w');
        outfile.write("#!/bin/bash\n")
        outfile.write("echo RUNNING make_frame_commands.sh...\n");
        outfile.write("cd FRAMES\n");
        outfile.write("if [ -z \"$(ls . )\" ]; then\n"); # if the directory is empty, then we make more frames. 
        outfile.write("  greadlink -f ../DATA/*.SAFE > data.list\n");  # I needed to change readlink --> greadlink for mac (readlink for linux)
        outfile.write("  echo \"%s %s 0\" > frames.ll\n" % (frame1_def[0], frame1_def[1]) );  # write the near-range edges of the frame to frame.ll
        outfile.write("  echo \"%s %s 0\" >> frames.ll\n" % (frame2_def[0], frame2_def[1]) );
        outfile.write("  make_s1a_frame.csh data.list frames.ll\n");
        outfile.write("  echo \"Assembling new frames!\"\n")
        outfile.write("else\n")
        outfile.write("  echo \"FRAMES contains files already... not assembling new frames\"\n")
        outfile.write("fi\n")
        outfile.write("cd ../\n");
        outfile.close();
        call(['chmod','+x','make_frame_commands.sh'],shell=False);
        call(['./make_frame_commands.sh'],shell=False)
        #call(['rm','make_frame_commands.sh'],shell=False)

        # Copy the scenes where only one scene is exactly covering the pre-defined frame (otherwise will be skipped because there's no combining to do)
        # It turns out that sometimes, the second scene that covers the frame doesn't exist, so there's only one scene for that given date.  
        # Other times, one scene covers the whole frame. 
        # Thankfully, make_s1a_frame.csh already copies the orbit files into the FRAMES directory, even if the .SAFE isn't copied. 

        # Make a list of dates in FRAMES/*.safe and compare with dates in the data directories. 
        call('compare_frames_acquisitions.sh',shell=True);
        print("Please check the frames and acquisitions and see if all your data has been included. ");
    return;



# --------------- STEP 0: Setting up raw_orig with safe, eof, xml, tiff ------------ # 
def manifest2raw_orig_eof(config_params):
	# This will set up the raw_orig directory from the DATA/.SAFE directories
	# Will also go into orbit directory and make copies of the right orbit files into the raw_orig directory. 

    if config_params.startstage>0:  # don't need to set up if we're starting mid-stream. 
        return;
    if config_params.endstage<0: # don't need to do this step
        return;

    file_list = get_SAFE_list(config_params);  # gets list of SAFE from FRAMES or DATA
    print(file_list);

    # When starting from preprocess, the system will often find a new super-master and re-align
    # To make this possible, we start from empty raw and raw-orig directories to avoid conflict.
    swath=config_params.swath
    print("Removing raw/ and raw_orig to proceed fresh. ");
    try:
         shutil.rmtree("F"+str(swath)+"/raw_orig");
    except OSError as e:
         print("Error: %s - %s." % (e.filename, e.strerror));
    try:
         shutil.rmtree("F"+str(swath)+"/raw");
    except OSError as e:
          print("Error: %s - %s." % (e.filename, e.strerror));

    # Unpack the .SAFE directories into raw_orig
    call(["mkdir","-p","F"+str(swath)])
    call(["mkdir","-p","F"+str(swath)+"/raw_orig"],shell=False);
    print("Copying xml files into raw_orig...")
    print("Copying manifest.safe files into raw_orig...")
    print("Copying tiff files into raw_orig...")
    # Copying these files is a lot of space, but it breaks if you only put the links to the files in the space. 
    for onefile in file_list:
        xml_files = sentinel_utilities.get_all_xml_names(onefile+'/annotation',config_params.polarization, config_params.swath);
        call(['cp',xml_files[0],'F'+str(swath)+'/raw_orig'],shell=False);
        manifest_safe_file = sentinel_utilities.get_manifest_safe_names(onefile);
        yyyymmdd=sentinel_utilities.get_date_from_xml(xml_files[0]);
        call(['cp',manifest_safe_file[0],'F'+str(swath)+'/raw_orig/'+yyyymmdd+'_manifest.safe'],shell=False);
        tiff_files = sentinel_utilities.get_all_tiff_names(onefile+'/measurement',config_params.polarization, config_params.swath);
        one_tiff_file=tiff_files[0].split("/")[-1];
        if not os.path.isfile('F'+str(swath)+'/raw_orig/'+one_tiff_file):
            call(['cp',tiff_files[0],'F'+str(swath)+'/raw_orig'],shell=False);  # only copy the tiff files if they don't already exist. 

        # STEP 2: get orbit files into the raw_orig directory
    for onefile in file_list:
        xml_name = glob.glob(onefile+'/annotation/*vv*.xml')[0];
        mydate = sentinel_utilities.get_date_from_xml(xml_name);
        sat    = sentinel_utilities.get_sat_from_xml(xml_name);
        eof_name = sentinel_utilities.get_eof_from_date_sat(mydate, sat, config_params.orbit_dir);
        print("Copying %s to raw_orig..." % eof_name);
        call(['cp',eof_name,'F'+str(swath)+'/raw_orig'],shell=False);
    print("copying s1a-aux-cal.xml to raw_orig...");
    call(['cp',config_params.orbit_dir+'/s1a-aux-cal.xml','F'+str(swath)+'/raw_orig'],shell=False);
    check_raw_orig_sanity(swath);
    return;

def check_raw_orig_sanity(swath):
    number_of_tiffs = int(check_output('ls F'+str(swath)+'/raw_orig/*.tiff | wc -l',shell=True));
    number_of_safes = int(check_output('ls F'+str(swath)+'/raw_orig/*.safe | wc -l',shell=True));
    number_of_EOFs  = int(check_output('ls F'+str(swath)+'/raw_orig/*.EOF | wc -l',shell=True));
    print('number of tiffs is %d ' % number_of_tiffs);
    print('number of safes is %d ' % number_of_safes);
    print('number of EOFs is %d ' % number_of_EOFs);
    if number_of_tiffs != number_of_safes:
        raise sentinel_utilities.Directory_error('Huge error: Your raw_orig directory has the wrong number of tiff/safe files. You should stop!'); 
    if number_of_tiffs != number_of_EOFs:
        raise sentinel_utilities.Directory_error('Huge error: Your raw_orig directory has the wrong number of tiff/EOF files. You should stop!');
    return;


# Return the files we're going to put into raw_orig
def get_SAFE_list(config_params):
    if config_params.frame1 != '':
        file_list = glob.glob("FRAMES/FRAME_1/*.SAFE");  # if we're assembling frames, we use the FRAMES directory. 
    else:
        file_list = glob.glob("DATA/*.SAFE");   # if we're not assembling frames, we use the DATA directory.     
    return file_list;


# --------------- STEP 1: Pre-processing (also aligning for Sentinel) ------------ # 
def preprocess(config_params):

    if config_params.startstage>1:  # don't need to pre-process if we're starting at topo or intf. 
        return;
    if config_params.endstage<1: # don't need to pre-process if we're just doing frames (stage 0)
        return;

    # # MODE 1: Before you know your super-master
    write_xml_prep(config_params.polarization, config_params.swath);   # writes the beginning, common part of README_prep.txt
    sentinel_utilities.make_data_in(config_params.polarization, config_params.swath, config_params.master);  # makes data.in the first time, with no super_master
    write_preproc_mode1(config_params.swath);              # writes the bottom of README_prep
    # call("./README_prep.txt",shell=True);  # This is the first time through- just get baseline plot to pick super-master.

    # Automatically decide on super-master and pop it to the front of data.in. 
    masterid = sentinel_utilities.choose_master_image(config_params.master, config_params.swath);
    print("master image is...")
    print(masterid);

    # MODE 2: Now you have picked a super-master. 
    # sentinel_utilities.write_super_master_batch_config(masterid);
    write_xml_prep(config_params.polarization, config_params.swath);      # writes the beginning, common part of README_prep.txt
    write_preproc_mode2(config_params.swath);                            # This writes the bottom of README_prep
    call("./README_prep.txt",shell=True); # This calls preproc_batch_tops.csh the second time.  Aligning will happen!  
    # NOTE: We automatically put the super-master into batch_tops.config (format is like [S1A20160829_ALL_F2])   
    return;

def write_xml_prep(polarization, swath):
    list_of_xml=sentinel_utilities.get_all_xml_names('F'+str(swath)+'/raw_orig',polarization,swath);
    print("Writing xmls in README_prep.txt");
    outfile=open("README_prep.txt",'w');
    outfile.write("# First, prepare the files.\n");
    outfile.write("cd F"+str(swath)+"\n");
    outfile.write("mkdir -p raw\n");
    outfile.write("cd raw\n");
    outfile.write("# in order to correct for Elevation Antenna Pattern Change, cat the manifest and aux files to the xmls\n");
    outfile.write("# delete the first line of the manifest file as it's not a typical xml file.\n\n");
    for item in list_of_xml:
        mydate=sentinel_utilities.get_date_from_xml(item);
        item=item.split("/")[-1];   # getting rid of the directory names
        outfile.write("awk 'NR>1 {print $0}' < ../raw_orig/"+mydate+"_manifest.safe > tmp_file\n");
        outfile.write("cat ../raw_orig/"+item+" tmp_file ../raw_orig/s1a-aux-cal.xml > ./"+item+"\n");
    outfile.write("rm tmp_file\n");
    outfile.write("ln -s ../raw_orig/*EOF .\nln -s ../raw_orig/*tiff .\nln -s ../topo/dem.grd .\n");
    outfile.write("cd ../../\n\n\n")
    outfile.close();
    return;

def write_preproc_mode1(swath):
    outfile=open("README_prep.txt",'a')
    outfile.write("cd F"+str(swath)+"\n");    
    outfile.write("cd raw\n");
    outfile.write("mv ../data.in .\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in ../topo/dem.grd 1'\n");
    outfile.write("preproc_batch_tops.csh data.in ../topo/dem.grd 1\n\n");
    outfile.write("cd ../../\n");
    outfile.close();
    print("Ready to call README_prep.txt in Mode 1.")
    call("chmod +x README_prep.txt",shell=True);
    return;

def write_preproc_mode2(swath):
    outfile=open("README_prep.txt",'a')
    outfile.write("cd F"+str(swath)+"\n");    
    outfile.write("cd raw\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in ../topo/dem.grd 2'\n");
    outfile.write("preproc_batch_tops.csh data.in ../topo/dem.grd 2\n\n");
    outfile.write("cd ../../\n");
    outfile.close();
    print("Ready to call README_prep.txt in Mode 2.")
    call("chmod +x README_prep.txt",shell=True);
    return;




# --------------- STEP 3: DEM topo2ra --------------- # 
def topo2ra(config_params):
    if config_params.startstage>3:  # if we're starting at intf, we don't do this. 
        return;
    if config_params.endstage<3:   # if we're ending at preproc, we don't do this. 
        return;
    call("sentinel_dem2topo_ra.csh "+config_params.config_file,shell=True);
    return;




# --------------- STEP 4: Make Interferograms ------------ # 
def make_interferograms(config_params):
    """
    1. form interferogram pairs from baseline_table
    2. make network plot
    3. write README_proc.txt
    """
    if config_params.startstage>4:  # if we're starting at sbas, we don't do this. 
        return;
    if config_params.endstage<4:   # if we're ending at topo, we don't do this. 
        return;

    [stems, times, baselines, missiondays] = sentinel_utilities.read_baseline_table('F'+str(config_params.swath)+'/raw/baseline_table.dat')
    intf_pairs=[];
    if config_params.ts_type=="SBAS" or config_params.ts_type=="NSBAS":
        intf_pairs_sbas = sentinel_utilities.get_small_baseline_subsets(stems, times, baselines, config_params.tbaseline, config_params.xbaseline, '', '');
        intf_pairs_manual = sentinel_utilities.get_manual_chain(stems, times, config_params.tbaseline); # ['20151118'] for mendocino
        intf_pairs = intf_pairs_sbas+intf_pairs_manual;
        print("README_proc.txt will be printed with tbaseline_max = "+str(config_params.tbaseline)+" days and xbaseline_max = "+str(config_params.xbaseline)+"m. ")
    elif config_params.ts_type=="CHAIN":
        intf_pairs = sentinel_utilities.get_chain_subsets(stems, times, baselines, config_params.bypass);
    elif config_params.ts_type=="MANUAL":
        intf_pairs = sentinel_utilities.get_manual_chain(stems, times, config_params.tbaseline, ['20151118']);
    else:
        print("config_params.ts_type is not a valid ts_type");
        sys.exit(1);




    # Here I will add the three pairs needed for the closure test manually. 
    intf_pairs=[];
    # intf_pairs = ["S1_20190411_ALL_F"+config_params.swath+":S1_20190423_ALL_F"+config_params.swath,
    #          "S1_20190330_ALL_F"+config_params.swath+":S1_20190411_ALL_F"+config_params.swath,
    #          "S1_20190330_ALL_F"+config_params.swath+":S1_20190423_ALL_F"+config_params.swath];

    # intf_pairs is the list of interferogram pairs made from SBAS or NSBAS or manual or chain 
    # Now we want to add the longer interferograms. 
    crit_days=30; crit_baseline=20;  # days, meters
    long_intfs_1 = rose_baseline_plot.compute_new_pairs(stems, times, baselines, crit_days, crit_baseline, 1);  # 1 year
    long_intfs_2 = rose_baseline_plot.compute_new_pairs(stems, times, baselines, crit_days, crit_baseline, 2);  # 2 years
    long_intfs_3 = rose_baseline_plot.compute_new_pairs(stems, times, baselines, crit_days, crit_baseline, 3);  # 2 years
    intf_all=intf_pairs + long_intfs_1 + long_intfs_2 + long_intfs_3;
    print(intf_all);

    # Make the stick plot of baselines 
    sentinel_utilities.make_network_plot(intf_all,stems,times, baselines, "F"+str(config_params.swath)+"/Total_Network_Geometry.eps");
    
    # Write the intf.in files
    # Writing to process interferograms. 
    outfile=open("README_proc.txt",'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch process Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("# First, create the files needed for intf_tops.csh\n\n");
    outfile.write("cd F"+str(config_params.swath)+"\n");
    outfile.write("ln -s ../batch.config .\n");
    outfile.write("rm intf*.in\n");
    for i,item in enumerate(intf_all):
        outfile.write('echo "' + item +'" >> intf_record.in\n');
    for i,item in enumerate(intf_all):
        outfile.write('echo "' + item +'" >> intf'+str(np.mod(i,config_params.numproc))+'.in\n');        
    outfile.write("\n# Process the interferograms.\n\n")
    outfile.write("ls intf?.in | parallel --eta 'intf_batch_tops_mod.csh {} "+config_params.config_file+"'\n\n\n");  # If you have parallel on your box
    # outfile.write("intf_batch_tops_mod.csh intf_record.in "+config_params.config_file+"\n\n\n");  # if you don't have parallel 
    outfile.write("cd ../\n");
    outfile.close();
    print("Ready to call README_proc.txt.")
    call("chmod +x README_proc.txt",shell=True);
    
    call("./README_proc.txt",shell=True);

    # print("Summarizing correlation for all interferograms.")
    # analyze_coherence.analyze_coherence_function();

    return;



# --------------- STEP 5: Unwrapping ------------ # 

def unwrapping(config_params):
    if config_params.startstage>5:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<5:   # if we're ending at intf, we don't do this. 
        return;   

    # Marie-Pierre's atmosphere correction 
    if config_params.detrend_atm_topo==1:
        flattentopo_driver.main_function();

    # Make plots of phasefilt.grd files. 
    # phasefilt_plot.top_level_driver('manual_remove.txt');
    # If you want to do this before and after flattentopo, you have to do it separately. 

    unwrap_sh_file="README_unwrap.txt";
    # sentinel_utilities.write_unordered_unwrapping(config_params.numproc, config_params.swath, unwrap_sh_file, config_params.config_file);
    sentinel_utilities.write_long_unwrapping(config_params.numproc, config_params.swath, unwrap_sh_file, config_params.config_file);
    # sentinel_utilities.write_select_unwrapping(config_params.numproc, config_params.swath, unwrap_sh_file, config_params.config_file);

    print("Ready to call "+unwrap_sh_file)
    call(['chmod','+x',unwrap_sh_file],shell=False);
    call("./"+unwrap_sh_file,shell=True);

    # call from the processing directory to place all unwrap.grd into single directory. 
    # print("Putting all unwrap.grd into a separate directory...")
    # file_of_interest='unwrap.grd';
    # call(['coalesce_intf_all_files.sh',file_of_interest],shell=False)  
    # sentinel_utilities.check_intf_all_sanity();
    return;




# --------------- STEP 6: Make SBAS ------------ # 
def do_timeseries(config_params):
    if config_params.startstage>6:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<6:   # if we're ending at intf, we don't do this. 
        return;

    # For later plotting, we want to project available GPS into LOS. 
    rowref=1524-621;
    colref=32;
    gps_into_LOS.top_level_driver(config_params, rowref, colref);

    # NOTE: 
    # Should copy batch.config into the nsbas directory
    # Should implement reference pixel at end of GACOS and unwrapping_errors

    return;





