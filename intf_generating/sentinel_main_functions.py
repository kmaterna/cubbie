import collections
import os, sys, argparse, configparser, glob
import numpy as np
from subprocess import call
from intf_generating import sentinel_utilities, rose_baseline_plot
from intf_atm_tools import flattentopo_driver
from stack_metrics import analyze_coherence

Params = collections.namedtuple('Params',
                                ['config_file', 'SAT', 'wavelength', 'startstage', 'endstage', 'master',
                                 'orbit_dir', 'DATA_dir', 'FRAMES_dir', 'intf_type', 'starttime', 'endtime',
                                 'tbaseline', 'xbaseline', 'annual_crit_days', 'annual_crit_baseline',
                                 'swath', 'polarization', 'atm_topo_detrend', 'desired_swaths',
                                 'frame1', 'frame2', 'numproc', 'threshold_snaphu']);

def read_config_argument_parsing():
    ################################################
    # Stage 0: Read and check config parameters
    #
    # read command line arguments and parse config file.
    parser = argparse.ArgumentParser(
        description='Run GMTSAR batch processing. Default automatically determines master, does alignment, '
                    'and runs desired interferograms.')
    parser.add_argument('config', type=str, help='supply name of config file to setup processing options. Required.')
    parser.add_argument('--debug', action='store_true', help='Print extra debugging messages (default: false)')
    args = parser.parse_args();
    config_file = args.config;
    config, config_params = read_config(config_file)
    return config_params;


def read_config(config_file):

    # read config file
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(config_file)

    # Setup gnu parallel multiprocessing tool
    numproc = config.getint('py-config', 'num_processors')

    # get options from config file
    config_file_orig = config_file;
    SAT = config.get('py-config', 'satellite')
    wavelength = config.getfloat('py-config', 'wavelength')
    startstage = config.getint('py-config', 'startstage')
    endstage = config.getint('py-config', 'endstage')
    master = config.get('csh-config', 'master_image').split()[0] \
        if len(config.get('csh-config', 'master_image')) > 0 else '';
    orbit_dir = config.get('py-config', 'orbit_dir')
    DATA_dir = config.get('py-config', 'DATA_dir')
    FRAMES_dir = config.get('py-config', 'FRAMES_dir')
    tbaseline = config.getint('py-config', 'max_timespan')
    xbaseline = config.getint('py-config', 'max_baseline')
    starttime = config.get('py-config', 'starttime');
    endtime = config.get('py-config', 'endtime');
    intf_type = config.get('py-config', 'intf_type');
    annual_crit_days = config.getint('py-config', 'annual_crit_days')
    annual_crit_baseline = config.getint('py-config', 'annual_crit_baseline')
    swath = config.get('py-config', 'swath')
    polarization = config.get('py-config', 'polarization')
    desired_swaths_temp = config.get('py-config', 'desired_subs')
    atm_topo_detrend = config.getint('py-config', 'atm_topo_detrend')
    frame_nearrange1 = config.get('py-config', 'frame_nearrange1')
    frame_nearrange2 = config.get('py-config', 'frame_nearrange2')
    threshold_snaphu = config.getfloat('csh-config', 'threshold_snaphu')
    threshold_geocode = config.getfloat('csh-config', 'threshold_geocode')

    print("Running sentinel batch processing, starting with stage %d" % startstage);

    # if master specified in the config file disagrees with existing data.in, we must re-do the pre-processing.
    if master and startstage > 1 and os.path.isfile('F' + str(swath) + '/raw/data.in'):
        master_datestr = sentinel_utilities.format_image_name_as_datestr(master);
        dataDotIn = np.genfromtxt('F' + str(swath) + '/raw/data.in', dtype='str')
        oldmaster = dataDotIn[0].split(':')[0];
        if master_datestr not in oldmaster:
            # sometimes, oldmaster is formatted like s1a-iw1-slc-vv-20171201t142317-...;
            # master_datestr is formatted like 20171213; master is S1_20171213_F1_ALL
            print('Error: The master specified in the config file disagrees with the old master in data.in. Exiting.')
            print('You should re-run starting from pre-processing with the master from the config file.')
            sys.exit(0);

    # if data.in is not found, we must do pre-processing.
    if startstage > 1 and not os.path.isfile('F' + str(swath) + '/raw/data.in'):
        print('Warning: Pre-processing has not been run, changing startstage to 1')
        startstage = 1

    # enforce startstage <= endstage
    if endstage < startstage:
        print('Warning: endstage is less than startstage. Setting endstage = startstage.')
        endstage = startstage

    # Turn '1,2,3' into ['1', '2', '3']
    desired_swaths = desired_swaths_temp.split(',');

    assert(threshold_geocode == 0), ValueError("Threshold_geocode should be 0 to skip geocoding. ")

    config_params = Params(config_file=config_file_orig, SAT=SAT, wavelength=wavelength, startstage=startstage,
                           endstage=endstage, master=master,
                           orbit_dir=orbit_dir, DATA_dir=DATA_dir, FRAMES_dir=FRAMES_dir,
                           tbaseline=tbaseline, xbaseline=xbaseline, starttime=starttime, endtime=endtime,
                           intf_type=intf_type, atm_topo_detrend=atm_topo_detrend,
                           desired_swaths=desired_swaths,
                           annual_crit_days=annual_crit_days, annual_crit_baseline=annual_crit_baseline,
                           swath=swath, polarization=polarization, frame1=frame_nearrange1,
                           frame2=frame_nearrange2, numproc=numproc, threshold_snaphu=threshold_snaphu);

    return config, config_params;


# --------------- STEP -1: Making frames from bursts ------------ #
def compile_frame_from_bursts(config_params):
    # Will assemble frames if pins are provided.
    # Would be good to check the total status of the download before doing this
    # (report_on_s1_data_holdings.py)
    # Works for all 3 swaths at once.
    if config_params.startstage > -1:  # don't need to set up if we're starting mid-stream.
        return;
    if config_params.endstage < -1:  # don't need to do this step
        return;

    if config_params.frame1:
        print("Assembling frames based on pins. ")
        pins_filename = config_params.FRAMES_dir+'/pins_ll.txt';
        auto_script = config_params.FRAMES_dir+'/auto_frame_commands.sh';

        # Here we want a frame to be made. We write the near-range pins.
        call(['mkdir', '-p', config_params.FRAMES_dir], shell=False);
        if len(glob.glob(config_params.FRAMES_dir+'/*.SAFE')) > 0:
            print("Looks like we already have frames assembled. End stage -1.");
            sentinel_utilities.compare_frames_with_safes(config_params);
            return;

        pins_file = open(pins_filename, 'w');
        pins_file.write(config_params.frame1.replace('/', ' ')+'\n');
        pins_file.write(config_params.frame2.replace('/', ' ')+'\n');
        pins_file.close();
        if config_params.polarization == 'vh':
            mode = 2;
        else:
            mode = 1;

        # Write:
        # --- data in a file, in chronological order, for each date
        # --- eof file
        # --- polarization
        outfile = open(auto_script, 'w');
        outfile.write("#!/bin/bash\n\n");
        dirlist, datelist = sentinel_utilities.get_all_safes_in_dir(config_params.DATA_dir);
        unique_datelist = sorted(set(datelist));
        for onedate in unique_datelist:
            ordered_safes = sentinel_utilities.get_safes_of_date(config_params.DATA_dir, onedate);
            satellite = ordered_safes[0].split('/')[-1][0:3]
            orbit_file = sentinel_utilities.get_eof_from_date_sat(onedate, satellite, config_params.orbit_dir);
            outfile.write("rm safes.txt\n")
            for item in ordered_safes:
                outfile.write("echo " + item + " >> safes.txt\n");
            outfile.write("create_frame_tops.csh safes.txt " + "../" + orbit_file + " pins_ll.txt "+str(mode)+"\n\n");
            # Putting ../ in front of orbit dir because we're down from the processing directory.
        outfile.close();

        call(['chmod', '+x', auto_script], shell=False);
        os.chdir(config_params.FRAMES_dir);
        call(['./auto_frame_commands.sh'], shell=False)
        os.chdir('../')

        # Manually copy the scenes where only one scene is exactly covering the pre-defined frame
        # (otherwise will be skipped because there's no combining to do)
        # It turns out that sometimes, the second scene that covers the frame doesn't exist, so there's only
        # one scene for that given date.  At other times, one scene covers the whole frame.

        # Output verification.
        sentinel_utilities.compare_frames_with_safes(config_params);
    else:
        print("No frames requested; moving on.")
    print("End stage -1.");
    return;


# --------------- STEP 0: Setting up raw_orig with safe, eof, xml, tiff ------------ #
def manifest2raw_orig_eof(config_params):
    # This will set up the raw_orig directory from the DATA/.SAFE directories
    # Will also go into orbit directory and make copies of the right orbit files into the raw_orig directory.
    # Happens for each swath.

    if config_params.startstage > 0:  # don't need to set up if we're starting mid-stream.
        return;
    if config_params.endstage < 0:  # don't need to do this step
        return;

    # get SAFEs from FRAMES/ or DATA/ and their associated datetimes
    safe_file_list, dt_list = sentinel_utilities.get_SAFE_list_for_raw_orig(config_params);

    # When starting from preprocess, the system will often find a new super-master and re-align
    swath = config_params.swath

    # Unpack the .SAFE directories into raw_orig
    call(["mkdir", "-p", "F" + swath + "/raw_orig"], shell=False);
    print("Copying xml files into raw_orig...")
    print("Copying manifest.safe files into raw_orig...")
    print("Copying tiff files into raw_orig...")
    # Copying these files is a lot of space, but it breaks if you only put the links to the files in the space.
    for onefile, onedt in zip(safe_file_list, dt_list):

        # Step 1: Get the names for tiff, xml, and eof files
        xmls, yyyymmdd = sentinel_utilities.get_all_xml_tiff_names(onefile + '/annotation', config_params.polarization,
                                                                   swath, filetype='xml');
        manifest_safe_file = onefile+'/manifest.safe';
        sat = sentinel_utilities.get_sat_from_xml(xmls[0]);
        eof_name = sentinel_utilities.get_eof_from_date_sat(onedt, sat, config_params.orbit_dir);

        # Copy various files
        call(['cp', xmls[0], 'F' + swath + '/raw_orig'], shell=False);
        call(['cp', manifest_safe_file, 'F' + swath + '/raw_orig/' + yyyymmdd[0] + '_manifest.safe'], shell=False);
        tiff_files, _ = sentinel_utilities.get_all_xml_tiff_names(onefile + '/measurement', config_params.polarization,
                                                                  swath, filetype='tiff');
        # only copy the tiff files if they don't already exist.
        one_tiff_file = tiff_files[0].split("/")[-1];
        if not os.path.isfile('F' + swath + '/raw_orig/' + one_tiff_file):
            call(['cp', tiff_files[0], 'F' + swath + '/raw_orig'], shell=False);

        # copy orbit files into the raw_orig directory
        print("Copying %s and associated tiff/xml/manifest.safe to raw_orig..." % eof_name);
        call(['cp', eof_name, 'F' + swath + '/raw_orig'], shell=False);
    print("copying s1a-aux-cal.xml to raw_orig...");
    call(['cp', config_params.orbit_dir + '/s1a-aux-cal.xml', 'F' + swath + '/raw_orig'], shell=False);
    sentinel_utilities.check_raw_orig_sanity(swath);
    return;


# --------------- STEP 1: Pre-processing (also aligning for Sentinel) ------------ #
def preprocess(config_params):
    if config_params.startstage > 1:  # don't need to pre-process if we're starting at topo or intf.
        return;
    if config_params.endstage < 1:  # don't need to pre-process if we're doing stage 0
        return;

    # Check: Stop if the number of SLCs is equal to the number of images
    outdir = 'F'+config_params.swath + '/raw'  # the output of this process
    if os.path.isdir(outdir):
        slc_list = glob.glob(outdir+'/*.SLC');
        datadotin = list(sentinel_utilities.np.genfromtxt(outdir+'/data.in', dtype='str'));
        if len(datadotin) > 0:
            assert(len(slc_list) < len(datadotin)), AssertionError("You have "+str(len(slc_list)) +
                                                                   " SLCs already; you want to keep them.");

    # Make data.in. Proper master not required.
    sentinel_utilities.write_data_in(config_params.polarization, config_params.swath, config_params.master,
                                     "F" + config_params.swath);
    write_xml_prep(config_params.polarization, config_params.swath);  # writes common part of README_prep.txt
    write_preproc_mode1(config_params.swath);  # writes the bottom of README_prep
    call(["./README_prep.txt"], shell=False);  # First time through- just get baseline plot and baseline table.

    # MODE 1: Before you know your super-master
    if config_params.master == "":
        # Automatically decide on super-master and pop it to the front of data.in.
        baseline_table_file = 'F'+config_params.swath+'/raw/baseline_table.dat';
        masterid = sentinel_utilities.choose_master_image(baseline_table_file);
        sentinel_utilities.write_data_in(config_params.polarization, config_params.swath, masterid,
                                         target_dir="F"+config_params.swath);
        sentinel_utilities.write_data_in(config_params.polarization, config_params.swath, masterid,
                                         target_dir="F"+config_params.swath+'/raw');
        print("master image is... "+masterid);
        sentinel_utilities.write_super_master_batch_config(masterid);  # automatically put super-master in batch.config
        print("Please check to confirm that you are happy with this. ")
        sys.exit(0);

    # MODE 2: Now you have picked a super-master.
    write_xml_prep(config_params.polarization, config_params.swath);  # writes common part of README_prep.txt
    write_preproc_mode2(config_params.swath);  # This writes the bottom of README_prep
    call(["./README_prep.txt"], shell=False);  # preproc_batch_tops.csh the second time.  Aligning will happen!
    return;


def write_xml_prep(polarization, swath):
    list_xml, list_datestrs = sentinel_utilities.get_all_xml_tiff_names('F'+swath+'/raw_orig', polarization, swath);
    print("Writing xmls in README_prep.txt");
    outfile = open("README_prep.txt", 'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# First, prepare the files.\n");
    outfile.write("cd F" + swath + "\n");
    outfile.write("mkdir -p raw\n");
    outfile.write("cd raw\n");
    outfile.write(
        "# in order to correct for Elevation Antenna Pattern Change, cat the manifest and aux files to the xmls\n");
    outfile.write("# delete the first line of the manifest file as it's not a typical xml file.\n\n");
    for item, mydate in zip(list_xml, list_datestrs):
        item = item.split("/")[-1];  # getting rid of the directory names
        outfile.write("awk 'NR>1 {print $0}' < ../raw_orig/" + mydate + "_manifest.safe > tmp_file\n");
        outfile.write("cat ../raw_orig/" + item + " tmp_file ../raw_orig/s1a-aux-cal.xml > ./" + item + "\n");
    outfile.write("rm tmp_file\n");
    outfile.write("ln -s ../raw_orig/*EOF .\nln -s ../raw_orig/*tiff .\nln -s ../topo/dem.grd .\n");
    outfile.write("cd ../../\n\n\n")
    outfile.close();
    return;


def write_preproc_mode1(swath):
    outfile = open("README_prep.txt", 'a');
    outfile.write("cd F"+swath+"\n");
    outfile.write("cd raw\n");
    outfile.write("mv ../data.in .\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in ../topo/dem.grd 1'\n");
    outfile.write("preproc_batch_tops.csh data.in ../topo/dem.grd 1\n\n");
    outfile.write("cd ../../\n");
    outfile.close();
    print("Ready to call README_prep.txt in Mode 1.")
    call(["chmod", "+x", "README_prep.txt"], shell=False);
    return;


def write_preproc_mode2(swath):
    outfile = open("README_prep.txt", 'a')
    outfile.write("cd F"+swath+"\n");
    outfile.write("cd raw\n");
    outfile.write("echo 'Calling preproc_batch_tops.csh data.in ../topo/dem.grd 2'\n");
    outfile.write("preproc_batch_tops.csh data.in ../topo/dem.grd 2\n\n");
    outfile.write("cd ../../\n");
    outfile.close();
    print("Ready to call README_prep.txt in Mode 2.")
    call(["chmod", "+x", "README_prep.txt"], shell=False);
    return;


# --------------- STEP 3: DEM topo2ra --------------- #
def topo2ra(config_params):
    if config_params.startstage > 3:  # if we're starting at intf, we don't do this.
        return;
    if config_params.endstage < 3:  # if we're ending at preproc, we don't do this.
        return;
    call(["sentinel_dem2topo_ra.csh", config_params.config_file], shell=False);
    return;


# --------------- STEP 4: Make Interferograms ------------ #

def get_total_intf_all(config_params):
    """Make a selection of interferograms to form.
    Hard coding which swath for consistency between swaths."""
    baseline_tuple_list = sentinel_utilities.read_baseline_table('F1/raw/baseline_table.dat');

    # Retrieving interferogram pairs based on settings in config_files. 
    intf_pairs = [];
    if "SBAS" in config_params.intf_type:
        intf_pairs = intf_pairs + sentinel_utilities.get_small_baseline_subsets(baseline_tuple_list,
                                                                                config_params.tbaseline,
                                                                                config_params.xbaseline);
    if "CHAIN" in config_params.intf_type:
        intf_pairs = intf_pairs + sentinel_utilities.get_chain_subsets(baseline_tuple_list);
    if "1YR" in config_params.intf_type:
        intf_pairs = intf_pairs + rose_baseline_plot.compute_new_pairs(baseline_tuple_list,
                                                                       config_params.annual_crit_days,
                                                                       config_params.annual_crit_baseline, 1);  # 1 yr
    if "2YR" in config_params.intf_type:
        intf_pairs = intf_pairs + rose_baseline_plot.compute_new_pairs(baseline_tuple_list,
                                                                       config_params.annual_crit_days,
                                                                       config_params.annual_crit_baseline, 2);  # 2 yrs
    if "3YR" in config_params.intf_type:
        intf_pairs = intf_pairs + rose_baseline_plot.compute_new_pairs(baseline_tuple_list,
                                                                       config_params.annual_crit_days,
                                                                       config_params.annual_crit_baseline, 3);  # 3 yrs
    if not intf_pairs:
        print("No intf_pairs found. Cannot make any interferograms.");
        print("Make sure config_params.intf_type is [combos of SBAS, CHAIN, 1YR, SBAS+CHAIN, etc.]");
        sys.exit(1);
    intf_pairs = sentinel_utilities.filter_intf_start_end(intf_pairs, config_params.starttime, config_params.endtime);
    intf_all = list(set(intf_pairs));  # removing duplicates. 
    print("Finding %d unique interferograms. " % len(intf_all));

    # Make the stick plot of baselines 
    sentinel_utilities.make_network_plot(intf_all, baseline_tuple_list,
                                         "F"+config_params.swath+"/Network_Geometry.eps");
    return intf_all;


def make_interferograms(config_params):
    """
    1. form interferogram pairs from baseline_table
    2. write README_proc.txt
    """
    if config_params.startstage > 4:  # if we're starting at sbas, we don't do this.
        return;
    if config_params.endstage < 4:  # if we're ending at topo, we don't do this.
        return;

    intf_all = get_total_intf_all(config_params);  # Make selection of interferograms to form.

    # Write the intf.in files
    outdir = "F" + str(config_params.swath + "/intf_all/");
    outfile = open("README_proc.txt", 'w');
    outfile.write("#!/bin/bash\n");
    outfile.write("# Script to batch process Sentinel-1 TOPS mode data sets.\n\n");
    outfile.write("# First, create the files needed for intf_tops.csh\n\n");
    outfile.write("cd F" + str(config_params.swath) + "\n");
    outfile.write("ln -s ../batch.config .\n");
    outfile.write("rm intf*.in\n");
    count = 0;
    for i, item in enumerate(intf_all):
        # Will only create interferograms that don't already exist.
        date1 = item[3:11];
        date2 = item[22:30];
        expected_folder = sentinel_utilities.ymd2yj(date1) + "_" + sentinel_utilities.ymd2yj(date2);
        if os.path.isfile(outdir + expected_folder + "/phasefilt.grd"):
            continue;
        else:
            # in case we're using one swath to generate intf_all for other swaths
            new_item = item.replace("_F1", "_F" + config_params.swath);
            outfile.write('echo "' + new_item + '" >> intf_record.in\n');
            outfile.write('echo "' + new_item + '" >> intf' + str(np.mod(count, config_params.numproc)) + '.in\n');
            count = count+1;
        # Write out in all cases.
    outfile.write("\n# Process the interferograms.\n\n")
    if int(config_params.numproc) > 1:   # parallel processing if you have GNU parallel on your box.
        outfile.write("ls intf?.in | parallel --eta 'intf_batch_tops_mod.csh {} "+config_params.config_file+"'\n\n\n");
    else:   # you don't have GNU parallel on your box
        outfile.write("intf_batch_tops_mod.csh intf_record.in "+config_params.config_file+"\n\n\n");
    outfile.write("cd ../\n");
    outfile.close();
    print("Ready to call README_proc.txt.")
    call(["chmod", "+x", "README_proc.txt"], shell=False);

    # The money line
    call(["./README_proc.txt"], shell=False);

    # print("Summarizing correlation for all interferograms.")
    # analyze_coherence.analyze_coherence_function();

    return;


# --------------- STEP 5: Unwrapping ------------ #

def unwrapping(config_params):
    """
    Unfortunately the parameters here often interact and cause code breaks, so you should proceed with caution
    """
    if config_params.startstage > 5:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 5:  # if we're ending at intf, we don't do this.
        return;

    unwrap_sh_file = "README_unwrap.txt";

    # Marie-Pierre's atmosphere correction, done before unwrapping
    if config_params.atm_topo_detrend == 1:
        if len(config_params.desired_swaths) == 1:  # for single swath
            merge_directory = "F" + str(config_params.desired_swaths[0]) + "/intf_all/";
            flattentopo_directory = "F" + str(config_params.desired_swaths[0]) + "/flattentopo/";
        else:   # Make merged wrapped files.  Doesn't need any metadata files except dem.grd.
            merge_directory = 'merged/'
            flattentopo_directory = 'merged_flattentopo/'
            sentinel_utilities.merge_wrapped(config_params.desired_swaths, config_params.master, merge_directory);

        flattentopo_example_rsc = flattentopo_directory + "/example_sd.int.rsc";    # set this up manually
        flattentopo_topora_grd = flattentopo_directory + "/topo_ra.grd";           # set this up manually
        flattentopo_driver.main_function(merge_directory, flattentopo_directory,
                                         flattentopo_topora_grd, flattentopo_example_rsc);
        # Then unwrap.  Some quick and dirty code here to get things running
        sentinel_utilities.write_unordered_unwrapping(config_params.numproc, config_params.swath, unwrap_sh_file,
                                                      config_params.config_file,
                                                      multiple_swaths=True if len(config_params.desired_swaths) > 1
                                                      else False);

    else:   # not using Marie-Pierre's code
        if len(config_params.desired_swaths) == 1:  # for single swath (might not be working)
            sentinel_utilities.write_unordered_unwrapping(config_params.numproc, config_params.swath, unwrap_sh_file,
                                                          config_params.config_file);

        else:  # For merging multiple swaths, write intfs and PRM files into inputfile, prepare for merging.
            common_intfs = sentinel_utilities.set_up_merge_unwrap(config_params.desired_swaths, 'merged');  # check path
            sentinel_utilities.write_merge_batch_input(common_intfs, config_params.master,
                                                       config_params.desired_swaths);
            sentinel_utilities.write_merge_unwrap(unwrap_sh_file);

    print("Ready to call " + unwrap_sh_file)
    call(['chmod', '+x', unwrap_sh_file], shell=False);
    call(["./"+unwrap_sh_file], shell=False);
    return;

# --------------- STEP 6: View Metrics ------------ #

def metrics(config_params):
    if config_params.startstage > 6:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 6:  # if we're ending at intf, we don't do this.
        return;

    print("Performing metrics on the stack.")
    igram_dir = 'merged_flattentopo/';
    raw_dir = 'F1/'
    analyze_coherence.analyze_coherence_function(igram_dir, raw_dir);
    return;
