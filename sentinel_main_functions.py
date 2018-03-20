import collections
import os,sys,argparse,time,configparser
import numpy as np
from subprocess import call
import glob
import sentinel_utilities

Params=collections.namedtuple('Params',['config_file','SAT','startstage','endstage','master','align_file','intf_file','orbit_dir','tbaseline','xbaseline','restart','mode','swath','polarization','frame1','frame2','numproc','ts_type','bypass']);

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
    frame_nearrange1=config.get('py-config','frame_nearrange1')
    frame_nearrange2=config.get('py-config','frame_nearrange2')
    ts_type=config.get('timeseries-config','ts_type')
    bypass=config.get('timeseries-config','bypass')
    
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
            #print('We will re-run starting from pre-processing with the master from the config file.')
            #startstage = 1
            #restart = True
    
    # if data.in is not found, we must do pre-processing.
    if startstage > 1 and not os.path.isfile('raw/data.in'):
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

    config_params=Params(config_file=config_file_orig, SAT=SAT,startstage=startstage,endstage=endstage,master=master,align_file=align_file,intf_file=intf_file,orbit_dir=orbit_dir,tbaseline=tbaseline, xbaseline=xbaseline,restart=restart,mode=mode,swath=swath,polarization=polarization,frame1=frame_nearrange1, frame2=frame_nearrange2, numproc=numproc, ts_type=ts_type, bypass=bypass);

    return config_params; 


def manifest2raw_orig_eof(config_params):
	# This will set up the raw_orig directory from the DATA/.SAFE directories
	# Will also go into orbit directory and make copies of the right orbit files into the raw_orig directory. 

    if config_params.startstage>1:  # don't need to set up if we're starting mid-stream. 
        return;

    file_list = get_frames_for_raw_orig(config_params);  # will assemble frames if necessary. Otherwise will just return the DATA/*.SAFE files
    print file_list;

    # Unpack the .SAFE directories into raw_orig
    call(["mkdir","-p","raw_orig"],shell=False);
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
            call(['cp',tiff_files[0],'raw_orig'],shell=False);  # only copy the tiff files if they don't already exist. 

        # STEP 2: get orbit files into the raw_orig directory
    for onefile in file_list:
        xml_name = glob.glob(onefile+'/annotation/*vv*.xml')[0];
        mydate = sentinel_utilities.get_date_from_xml(xml_name);
        sat    = sentinel_utilities.get_sat_from_xml(xml_name);
        eof_name = sentinel_utilities.get_eof_from_date_sat(mydate, sat, config_params.orbit_dir);
        print "Copying %s to raw_orig..." % eof_name;
        call(['cp',eof_name,'raw_orig'],shell=False);
    print "copying s1a-aux-cal.xml to raw_orig..."
    call(['cp',config_params.orbit_dir+'/s1a-aux-cal.xml','raw_orig'],shell=False);
    return;


def get_frames_for_raw_orig(config_params):
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
        outfile.write("cd FRAMES\n");
        outfile.write("if [ -z \"$(ls -A $1 )\" ]; then\n"); # if the directory is empty, then we make more frames. 
        outfile.write("  readlink -f ../DATA/*.SAFE > data.list\n");
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
        call(['rm','make_frame_commands.sh'],shell=False)

        # Copy the scenes where only one scene is exactly covering the pre-defined frame (otherwise will be skipped because there's no combining to do)
        # It turns out that sometimes, the second scene that covers the frame doesn't exist, so there's only one scene for that given date.  
        # Other times, one scene covers the whole frame. 
        # I haven't figured out a way to automate this quite yet. 
        # Thankfully, make_s1a_frame.csh already copies the orbit files into the FRAMES directory, even if the .SAFE isn't copied. 

        # Might as well make a list of dates in FRAMES/*.safe and compare with dates in the data directories. 
        # We already have the GMT script that plots this... 

        call('compare_frames_acquisitions.sh',shell=True);
        print "Please check the frames and acquisitions and see if all your data has been included. "


        file_list = glob.glob("FRAMES/FRAME_1/*.SAFE");  # if we're assembling frames, we use the FRAMES directory. 
    else: 
        file_list = glob.glob("DATA/*.SAFE");   # if we're not assembling frames, we use the DATA directory.     
    return file_list;


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

    [stems, times, baselines, missiondays] = sentinel_utilities.read_baseline_table('raw/baseline_table.dat')
    if config_params.ts_type=="SBAS":
        startdate="2016243";
        enddate="2016293";
        #startdate="2016291"
        #enddate="";
        intf_pairs = sentinel_utilities.get_small_baseline_subsets(stems, times, baselines, config_params.tbaseline, config_params.xbaseline, startdate, enddate);
        print "README_proc.txt will be printed with tbaseline_max = "+str(config_params.tbaseline)+" days and xbaseline_max = "+str(config_params.xbaseline)+"m. "
    
    elif config_params.ts_type=="CHAIN":
        intf_pairs = sentinel_utilities.get_chain_subsets(stems, times, baselines, config_params.bypass);
    else:
        print "config_params.ts_type is not a valid ts_type";
        sys.exit(1);

    # Make the stick plot of baselines 
    sentinel_utilities.make_network_plot(intf_pairs,stems,times, baselines);

    # Writing to process interferograms. 
    outfile=open("README_proc.txt",'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch process Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("# First, create the files needed for intf_tops.csh\n\n");
    outfile.write("rm intf*.in\n");
    for i,item in enumerate(intf_pairs):
        outfile.write('echo "' + item +'" >> intf_record.in\n');
    for i,item in enumerate(intf_pairs):
        outfile.write('echo "' + item +'" >> intf'+str(np.mod(i,config_params.numproc))+'.in\n');        
    outfile.write("\n# Process the interferograms.\n\n")
    outfile.write("ls intf?.in | parallel --eta 'intf_batch_tops_km.csh {} "+config_params.config_file+"'\n\n\n");
    outfile.close();
    print "Ready to call README_proc.txt."
    call("chmod +x README_proc.txt",shell=True);
    call("./README_proc.txt",shell=True);

    print "Summarizing correlation for all interferograms."
    #call("get_corr_all_intfs.sh",shell=True);

    return;



# --------------- STEP 5: Unwrapping ------------ # 

def unwrapping(config_params):
    if config_params.startstage>5:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<5:   # if we're ending at intf, we don't do this. 
        return;   

    call("rm intf?.in",shell=True);
    unwrap_sh_file="README_unwrap.txt";
    sentinel_utilities.write_ordered_unwrapping(config_params.numproc, unwrap_sh_file, config_params.config_file);

    print "Ready to call "+unwrap_sh_file
    call(['chmod','+x',unwrap_sh_file],shell=False);
    call("./"+unwrap_sh_file,shell=True);

    return;




# --------------- STEP 6: Make SBAS ------------ # 
def do_sbas(config_params):

    if config_params.startstage>6:  # if we're starting after, we don't do this. 
        return;
    if config_params.endstage<6:   # if we're ending at intf, we don't do this. 
        return;

    [stems,tbaseline,xbaseline,mission_days]=sentinel_utilities.read_baseline_table('raw/baseline_table.dat');
    t_int=[];
    for t in tbaseline:
        t_int.append(round(float(t)));  # a list of integers like 2016214 for 2016-day-214.

    mission_days_sorted=[x for (y, x) in sorted(zip(t_int,mission_days))];
    t_int.sort();
    tbaseline.sort();

    intf_computed=sentinel_utilities.glob_intf_computed();  # looks like a list of labels like 2016217_2016205
    n_intf=len(intf_computed);
    outfile=open("README_sbas.txt",'w');
    outfile.write("# First, prepare the input files needed for sbas\n#\n");
    outfile.write("rm -f SBAS\nmkdir SBAS\ncd SBAS\nrm intf.tab scene.tab\n\n\n");
    outfile.write("# based on baseline_table.dat create the intf.tab and scene.tab for sbas\n");

    # writing intf.tab
    outfile.write("# phase  corherence  ref_id  rep_id  baseline\n")
    for img_pair in intf_computed:
        first_image=img_pair[0:7]
        second_image=img_pair[8:]
        for a, b in zip(xbaseline, t_int):
            if abs(int(np.floor(b)) - int(first_image))<=1:
                # print "first image found";
                # print int(np.floor(b))
                # print int(first_image)
                master_xbaseline=a;
            if abs(int(np.floor(b)) - int(second_image))<=1:
                slave_xbaseline=a;
                # print "second image found";
                # print int(np.floor(b))
                # print int(second_image)                
        total_baseline=slave_xbaseline - master_xbaseline;
        outfile.write('echo "../intf_all/'+img_pair+'/unwrap.grd ')
        outfile.write('../intf_all/'+img_pair+'/corr.grd ')
        outfile.write(first_image+' '+second_image+' ')
        outfile.write(str(total_baseline))
        outfile.write('" >> intf.tab\n');
    outfile.write("#\n\n");

    # writing scene.tab (only the scenes that are actually used in SBAS processing)
    # Right now this isn't producing anything. 
    outfile.write("# scene_id  day\n");
    scenes_used='';
    n_scenes=0;
    scenes_used=[];
    for intf in intf_computed:
        scenes_used.append(intf[0:7]);  # catch which scenes are actually used in SBAS processing. 
        scenes_used.append(intf[8:15]);
    for x in range(len(tbaseline)):
        temp = tbaseline[x];
        tempint = int(np.round(temp))
        if str(tempint) in scenes_used or str(tempint+1) in scenes_used or str(tempint-1) in scenes_used:
            outfile.write('echo "'+str(tempint)+' '+mission_days_sorted[x]+'" >> scene.tab\n');
            n_scenes+=1;
    outfile.write("\n\n");

    intf_ex=intf_computed[0];  # an example interferogram where we get the geographic coordinates for grdinfo
    outfile.write("xdim=`gmt grdinfo -C ../intf_all/"+intf_ex+"/unwrap.grd | awk '{print $10}'`\n");
    outfile.write("ydim=`gmt grdinfo -C ../intf_all/"+intf_ex+"/unwrap.grd | awk '{print $11}'`\n\n\n");

    outfile.write("# run sbas\n");
    outfile.write("sbas intf.tab scene.tab "+str(n_intf)+" "+str(n_scenes)+" $xdim $ydim -smooth 1.0 -wavelength 0.0554658 -incidence 30 -range 800184.946186 -rms -dem\n\n\n")

    outfile.write("# project the velocity to Geocooridnates\n");
    outfile.write('echo "writing to georeferenced coordinates..."\n');
    outfile.write("ln -s ../topo/trans.dat .\n");
    outfile.write("proj_ra2ll.csh trans.dat vel.grd vel_ll.grd\n");
    outfile.write("gmt grd2cpt vel_ll.grd -T= -Z -Cjet > vel_ll.cpt\n");
    outfile.write("grd2kml.csh vel_ll vel_ll.cpt\n");
    outfile.write("cd ..\n\n");
    outfile.write('echo "SBAS operation performed!"\n\n')

    outfile.close();
    print "README_sbas.txt written. Ready to call README_sbas.txt."
    call("chmod u+x README_sbas.txt",shell=True);
    call("./README_sbas.txt",shell=True); # Make sbas!
    return;





