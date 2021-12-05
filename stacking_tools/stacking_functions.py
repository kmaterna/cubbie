import re, glob
from subprocess import call
from . import stacking_utilities, nsbas_accessing, coseismic_stack, stack_corr, workflow_isce_with_uavsar
from . import Super_Simple_Stack as sss


# --------------- STEP 0: Setting up ------------ # 
def set_up_output_directories(config_params):
    if config_params.startstage > 0:
        return;
    if config_params.endstage < 0:
        return;
    print("\nStart Stage 0 - Setting up output directories");
    call(['mkdir', '-p', config_params.ts_output_dir], shell=False);
    print('calling: mkdir -p %s' % config_params.ts_output_dir);
    call(['cp', config_params.config_file, config_params.ts_output_dir], shell=False);
    print('calling: cp %s %s' % (config_params.config_file, config_params.ts_output_dir));
    if config_params.skip_file:
        call(['cp', config_params.skip_file, config_params.ts_output_dir], shell=False);
    print("End Stage 0 - Setting up output directories (%s) \n" % config_params.ts_output_dir);
    return;


# --------------- STEP 1: Make corrections ------------ # 
def make_corrections(config_params):
    if config_params.startstage > 1:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 1:  # if we're ending at intf, we don't do this.
        return;
    print("Start Stage 1 - optional atm and unwrapping corrections");
    if config_params.custom_unwrapping:
        workflow_isce_with_uavsar.custom_isce_unwrapping(config_params);
    # This is where we would implement GACOS, APS, topo-detrending, or unwrapping errors if we had them.
    print("End Stage 1 - optional atm corrections\n");
    return;


# --------------- STEP 2: Get Reference Pixel ------------ #
def get_ref(config_params):
    if config_params.startstage > 2:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 2:  # if we're ending at intf, we don't do this.
        return;
    print("Start Stage 2 - Finding Files and Reference Pixel");

    # Very general, returns the (d1, d2, filenames) tuples of all interferograms; doesn't discriminate
    intfs = stacking_utilities.get_list_of_intf_all(config_params, returnval='intf_files');
    intfs = stacking_utilities.exclude_intfs_manually(intfs, config_params.skip_file);

    # Here we get ref_idx if we don't have it already
    stacking_utilities.get_ref_index(config_params.ref_loc, config_params.ref_idx, config_params.geocoded_intfs, intfs,
                                     config_params.ts_output_dir+config_params.signal_spread_filename);

    print("End Stage 2 - Finding Files and Reference Pixel\n");
    return;


# --------------- STEP 3: Velocities and Time Series! ------------ # 
def vels_and_ts(config_params):
    if config_params.startstage > 3:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 3:  # if we're ending at intf, we don't do this.
        return;

    print("Start Stage 3 - Velocities and Time Series");
    call(['cp', config_params.config_file, config_params.ts_output_dir], shell=False);

    # This is where hand-picking takes place: manual excludes, long intfs only, ramp-removed, atm-removed, etc.
    intf_files, corr_files = stacking_utilities.make_selection_of_intfs(config_params);

    # Make signal spread after excludes have taken place.
    # Beginning of refactor is here.
    if config_params.make_signal_spread:
        stack_corr.drive_signal_spread_calculation(corr_files, config_params.signal_coh_cutoff,
                                                   config_params.ts_output_dir, config_params.signal_spread_filename);

    if config_params.ts_type == "STACK":
        print("\nRunning velocities by simple stack.")
        sss.drive_velocity_simple_stack(config_params, intf_files);
    if config_params.ts_type == "COSEISMIC":
        print("\nMaking a simple coseismic stack");
        coseismic_stack.drive_coseismic_stack(config_params, intf_files);
    if config_params.ts_type == "NSBAS" or config_params.ts_type == "WNSBAS":
        print("\nRunning velocities and time series by NSBAS or WNSBAS");
        nsbas_accessing.nsbas_ts_format_selector(config_params, intf_files, corr_files);

    print("End Stage 3 - Velocities and Time Series\n");
    return;


# --------------- STEP 4: Geocoding Velocities ------------ # 
def geocode_vels(config_params):
    if config_params.startstage > 4:  # if we're starting after, we don't do this.
        return;
    if config_params.endstage < 4:  # if we're ending at intf, we don't do this.
        return;
    print("Start Stage 4 - Geocoding");
    if config_params.SAT == "UAVSAR":
        workflow_isce_with_uavsar.geocode_isce_uavsar(config_params);

    # Then, quickly geocode all the time series files. 
    # Call from the processing directory
    grd_source_directory = config_params.ts_output_dir + "/"
    filelist = glob.glob(grd_source_directory + "????????.grd");
    for i in range(len(filelist)):
        datestr = re.findall(r"\d\d\d\d\d\d\d\d", filelist[i])[0];
        print(datestr);
        call(["quick_geocode.csh", grd_source_directory, "merged", datestr + ".grd", datestr + "_ll"], shell=False);

    print("End Stage 4 - Geocoding");
    return;
