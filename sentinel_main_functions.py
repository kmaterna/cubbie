import collections
import os,sys,argparse,time,configparser
import numpy as np
from subprocess import call
import glob
import sentinel_utilities

Params=collections.namedtuple('Params',['config_file','SAT','startstage','endstage','master','align_file','intf_file','orbit_dir','tbaseline','xbaseline','restart','mode','swath','polarization']);

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
    
    # Setup MPI (optional) or python multiprocessing pool
    if args.mpi:
        from mpi4py import MPI
        import mpi4py_map
        comm = MPI.COMM_WORLD
        ver  = MPI.Get_version()
        numproc = comm.Get_size()
        rank = comm.Get_rank()
        fstSec  = MPI.Wtime()
    else:
        import multiprocessing
        numproc=config.getint('py-config','num_processors')
        rank = 0    

    # get options from config file
    config_file_orig=args.config;
    SAT=config.get('py-config','satellite')
    startstage=config.getint('py-config','startstage')
    endstage=config.getint('py-config','endstage')    
    master=config.get('csh-config','master_image')
    align_file=config.get('py-config','align_file')
    intf_file=config.get('py-config','intf_file')
    orbit_dir=config.get('py-config','orbit_dir')
    tbaseline=config.getint('py-config','max_timespan')
    xbaseline=config.getint('py-config','max_baseline') 
    restart=config.getboolean('py-config','restart')
    mode=config.get('py-config','mode') 
    swath=config.get('py-config','swath') 
    polarization=config.get('py-config','polarization') 
    
    # print config options
    if args.debug:
        print('Running gmtsar_app.py:')
        if args.mpi:
            print('  Using MPI')
        else:
            print('  Not using MPI')
        print('  config file:',args.config,'contains the following options')    
        print(config.write(sys.stdout))
    
    # check config options
    
    #default names for files
    if align_file == '':
        align_file = 'align_batch.in'
        
    if intf_file == '':
        intf_file = 'intf_batch.in'
    
    # only valid option for mode is 'scan'    
    if mode != '' and mode != 'scan':
        print('warning: invalid option mode = %s in config file, using blank'%mode)
        mode = ''

    # if master specified in the config file disagrees with existing data.in, we must re-do the pre-processing.
    if master and startstage > 1 and os.path.isfile('raw/data.in'):
        #check master in data.in
        dataDotIn=np.genfromtxt('raw/data.in',dtype='str')
        oldmaster=dataDotIn[0]
        if '-'+master[3:11]+'t' not in oldmaster:  # For sentinel, oldmaster is formatted like s1a-iw1-slc-vv-20171201t142317-...; master is formatted like S1A20171213_ALL_F1
            print('Warning: The master specified in the config file disagrees with the old master in data.in.')
            print('We will re-run starting from pre-processing with the master from the config file.')
            startstage = 1
            restart = True
    
    # if data.in is not found, we must do pre-processing.
    if startstage > 1 and not os.path.isfile('raw/data.in'):
        print('Warning: Pre-processing has not been run, changing startstage to 1')
        startstage = 1
            
    # enforce startstage <= endstage
    if endstage < startstage:
        print('Warning: endstage is less than startstage. Setting endstage = startstage.')
        endstage = startstage
    
    # write the new config file to use for processing

    if rank == 0:
        #logtime is the timestamp added to all logfiles created during this run
        logtime=time.strftime("%Y_%m_%d-%H_%M_%S")
        config_file='batch.run.'+ logtime +'.cfg'
        with open(config_file, 'w') as configfilehandle:
            config.write(configfilehandle)
    else:
        logtime = None
        config_file = None
    if args.mpi:
        logtime = comm.bcast(logtime,root=0)
        config_file = comm.bcast(config_file,root=0)

    config_params=Params(config_file=config_file_orig, SAT=SAT,startstage=startstage,endstage=endstage,master=master,align_file=align_file,intf_file=intf_file,orbit_dir=orbit_dir,tbaseline=tbaseline, xbaseline=xbaseline,restart=restart,mode=mode,swath=swath,polarization=polarization);

    return config_params; 


def manifest2raw_orig_eof(config_params):
	# This will set up the raw_orig directory from the DATA/.SAFE directories
	# Will also go into orbit directory and make copies of the right orbit files into the raw_orig directory. 

        if config_params.startstage>1:  # don't need to set up if we're starting mid-stream. 
            return;

        # Unpack the .SAFE directories into raw_orig
        call(["mkdir","-p","raw_orig"],shell=False);
	file_list = glob.glob("DATA/*.SAFE");
        print file_list;
        print "Copying xml files into raw_orig..."
        print "Copying manifest.safe files into raw_orig..."
        print "Copying tiff files into raw_orig..."
        # Copying these files is a lot of space, but it breaks if you only put the links to the files in the space. 
        for onefile in file_list:
            xml_files = sentinel_utilities.get_all_xml_names(onefile+'/annotation',config_params.polarization, config_params.swath);
            call(['cp',xml_files[0],'raw_orig'],shell=False);
            manifest_safe_file = sentinel_utilities.get_manifest_safe_names(onefile);
            yyyymmdd=sentinel_utilities.get_date_from_xml(xml_files[0]);
            call(['cp',manifest_safe_file[0],'raw_orig/'+yyyymmdd+'_manifest.safe'],shell=False);
            tiff_files = sentinel_utilities.get_all_tiff_names(onefile+'/measurement',config_params.polarization, config_params.swath);
            one_tiff_file=tiff_files[0].split("/")[-1];
            if not os.path.isfile('raw_orig/'+one_tiff_file):
                call(['cp',tiff_files[0],'raw_orig'],shell=False);

        # STEP 2: get orbit files into the raw_orig directory
        for onefile in file_list:
            xml_name = glob.glob(onefile+'/annotation/*vv*.xml')[0];
            eof_name = sentinel_utilities.get_eof_from_xml(xml_name, config_params.orbit_dir);
            print "Copying %s to raw_orig..." % eof_name;
            call(['cp',eof_name,'raw_orig'],shell=False);
	print "copying s1a-aux-cal.xml to raw_orig..."
        call(['cp',config_params.orbit_dir+'/s1a-aux-cal.xml','raw_orig'],shell=False);
	return;


# --------------- STEP 1: Pre-processing (also aligning for Sentinel) ------------ # 
def preprocess(config_params):

    if config_params.startstage>1:  # don't need to pre-process if we're starting at topo or intf. 
        return;

    # MODE 1: Before you know your super-master
    write_xml_prep(config_params.polarization, config_params.swath);   # writes the beginning, common part of README_prep.txt
    sentinel_utilities.make_data_in(config_params.polarization, config_params.swath, config_params.master);  # makes data.in the first time, with no super_master
    write_preproc_mode1();              # writes the bottom of README_prep
    call("./README_prep.txt",shell=True);  # This is the first time through- just get baseline plot to pick super-master.

    # Automatically decide on super-master and pop it to the front of data.in. 
    masterid = sentinel_utilities.choose_master_image();
    print "master image is..."
    print masterid;

    # MODE 2: Now you have picked a super-master. 
    sentinel_utilities.write_super_master_batch_config(masterid);
    write_xml_prep(config_params.polarization, config_params.swath);                 # writes the beginning, common part of README_prep.txt
    write_preproc_mode2();                            # This writes the bottom of README_prep
    call("./README_prep.txt",shell=True); # This calls preproc_batch_tops.csh the second time.  Aligning will happen!  
    # NOTE: We automatically put the super-master into batch_tops.config (format is like [S1A20160829_ALL_F2])   
    return;



def write_xml_prep(polarization, swath):
    list_of_xml=sentinel_utilities.get_all_xml_names('raw_orig',polarization,swath);
    print "Writing xmls in README_prep.txt";	
    outfile=open("README_prep.txt",'w');
    outfile.write("# First, prepare the files.\n");
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
    outfile.write("cd ../\n\n\n")
    outfile.close();
    return;

def write_preproc_mode1():
    outfile=open("README_prep.txt",'a')
    outfile.write("cd raw\n");
    outfile.write("mv ../data.in .\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in dem.grd 1'\n");
    outfile.write("preproc_batch_tops.csh data.in dem.grd 1\n\n");
    outfile.write("cd ../\n");
    #outfile.write("cp raw/baseline_table.dat .\n"); # I find these two lines sort of annoying. 
    #outfile.write("cp raw/baseline.ps .\n\n")
    outfile.close();
    print "Ready to call README_prep.txt in Mode 1."
    call("chmod +x README_prep.txt",shell=True);
    return;

def write_preproc_mode2():
    outfile=open("README_prep.txt",'a')
    outfile.write("cd raw\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in dem.grd 2'\n");
    outfile.write("preproc_batch_tops.csh data.in dem.grd 2\n\n");
    outfile.write("cd ../\n");
    outfile.close();
    print "Ready to call README_prep.txt in Mode 2."
    call("chmod +x README_prep.txt",shell=True);
    return;




# --------------- STEP 3: DEM topo2ra --------------- # 
def topo2ra(config_params):
    if config_params.startstage>3:  # if we're starting at intf, we don't do this. 
        return;
    if config_params.endstage<3:   # if we're ending at preproc, we don't do this. 
        return;
    return;



# --------------- STEP 4: Make Interferograms ------------ # 
def make_interferograms(config_params):
    """
    1. form interferogram pairs from baseline_table
    2. make network plot
    3. write README_proc.txt
    """
    if config_params.endstage<4:  # if we're ending at topo, we don't need to do this. 
        return;

    baselineFile = np.genfromtxt('raw/baseline_table.dat',dtype=str)
    stems = baselineFile[:,0].astype(str)
    times = baselineFile[:,1].astype(float)
    baselines = baselineFile[:,4].astype(float)
    intf_pairs = sentinel_utilities.get_small_baseline_subsets(stems, times, baselines, config_params.tbaseline, config_params.xbaseline)

    # Make the stick plot of baselines 
    sentinel_utilities.make_network_plot(intf_pairs,stems,times, baselines);

    # Writing to process interferograms. 
    outfile=open("README_proc.txt",'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch process Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("# First, create the files needed for intf_tops.csh\n\n");
    outfile.write("rm -f intf.in\nrm -r intf intf_all\n\n");
    for item in intf_pairs:
        outfile.write('echo "' + item +'" >> intf.in\n');
    outfile.write("\n# Process the interferograms, remember to set your super master in the batch.config file.\n\n")
    outfile.write("intf_tops.csh intf.in "+config_params.config_file+"\n\n\n");
    outfile.close();
    print "README_proc.txt printed with tbaseline_max = "+str(config_params.tbaseline)+" days and xbaseline_max = "+str(config_params.xbaseline)+"m. "
    print "Ready to call README_proc.txt."
    call("chmod +x README_proc.txt",shell=True);
    call("./README_proc.txt",shell=True);

    return;



