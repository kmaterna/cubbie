# Stacking config parser
import argparse
import configparser
import sys
import datetime as dt
import collections

Params = collections.namedtuple('Params',
                                ['config_file', 'SAT', 'wavelength', 'startstage', 'endstage', 'ref_loc', 'ref_idx',
                                 'ts_type', 'file_format',
                                 'custom_unwrapping', 'detrend_atm_topo', 'gacos', 'aps', 'dem_error',
                                 'sbas_smoothing', 'ts_format', 'make_signal_spread', 'signal_coh_cutoff', 
                                 'nsbas_min_intfs', 'intf_filename', 'corr_filename', 'geocoded_intfs', 'baseline_file',
                                 'start_time', 'end_time', 'coseismic', 'intf_timespan', 'gps_file', 'flight_angle',
                                 'look_angle', 'skip_file', 'signal_spread_filename',
                                 'intf_dir', 'ts_points_file', 'ts_output_dir'])


Params_custom = collections.namedtuple('Params_custom', ['config_file', 'rlks', 'alks', 'filt', 'cor_cutoff_mask',
                                                         'xbounds', 'ybounds', 'llh_file', 'lkv_file',
                                                         'ts_output_dir'])


# ----------------------------- #

def parse_cmd_and_config():
    ################################################
    # Stage 0: Read and check config parameters
    # read command line arguments and parse config file.
    parser = argparse.ArgumentParser(description='Run stack processing. Either config or print_config are required.')
    parser.add_argument('--config', type=str, help='supply name of config file to setup processing options.')
    parser.add_argument('--print_config', help='file system location where to print example config file.')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    if args.config is not None:
        config_file = args.config
        config, config_params = read_config_general(config_file)
        return config, config_params
    elif args.print_config is not None:
        directory = args.print_config
        print_default_config(directory+"/stacking_default_config.txt")
        sys.exit(0)


def read_config_isce(config_file):
    """
    ISCE-specific parameters
    read config file
    """
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(config_file)

    rlks = config.getint('py-config', 'rlks') if config.has_option('py-config', 'rlks') else 0
    alks = config.getint('py-config', 'alks') if config.has_option('py-config', 'alks') else 0
    filt = config.getfloat('py-config', 'filt') if config.has_option('py-config', 'filt') else 0
    cor_cutoff_mask = config.getfloat('py-config', 'cor_cutoff_mask') if \
        config.has_option('py-config', 'cor_cutoff_mask') else 1
    xbounds = config.get('py-config', 'xbounds') if config.has_option('py-config', 'xbounds') else '0/100'
    ybounds = config.get('py-config', 'ybounds') if config.has_option('py-config', 'ybounds') else '0/100'
    llh_file = config.get('py-config', 'llh_file') if config.has_option('py-config', 'llh_file') else ''
    lkv_file = config.get('py-config', 'lkv_file') if config.has_option('py-config', 'lkv_file') else ''
    ts_output_dir = config.get('py-config', 'ts_output_dir') if config.has_option('py-config', 'ts_output_dir') else ''
    ts_output_dir = ts_output_dir + "_" + str(rlks) + "_" + str(alks) + "_" + str(filt)  # custom output directory
    my_Params = Params_custom(config_file=config_file, rlks=rlks, alks=alks, filt=filt,
                              cor_cutoff_mask=cor_cutoff_mask, xbounds=xbounds, ybounds=ybounds, llh_file=llh_file,
                              lkv_file=lkv_file, ts_output_dir=ts_output_dir)
    return my_Params


def read_config_general(config_file):
    """
    Common to both GMTSAR and ISCE.
    read config file
    """
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(config_file)

    # get options from config file
    SAT = config.get('py-config', 'satellite')
    wavelength = config.getfloat('py-config', 'wavelength')
    startstage = config.getint('py-config', 'startstage')
    endstage = config.getint('py-config', 'endstage')
    ref_loc = config.get('py-config', 'ref_loc') if (
        config.has_option('py-config', 'ref_loc')) else ''
    ref_idx = config.get('py-config', 'ref_idx') if (
        config.has_option('py-config', 'ref_idx')) else ''
    custom_unwrapping = config.getint('py-config', 'custom_unwrapping') if (
        config.has_option('py-config', 'custom_unwrapping')) else 0
    gacos = config.getint('py-config', 'gacos') if (
        config.has_option('py-config', 'gacos')) else 0
    aps = config.getint('py-config', 'aps') if (
        config.has_option('py-config', 'aps')) else 0
    dem_error = config.getint('py-config', 'dem_error') if (
        config.has_option('py-config', 'dem_error')) else 0
    detrend_atm_topo = config.getint('py-config', 'detrend_atm_topo') if (
        config.has_option('py-config', 'detrend_atm_topo')) else 0
    nsbas_min_intfs = config.getfloat('py-config', 'nsbas_min_intfs') if (
        config.has_option('py-config', 'nsbas_min_intfs')) else 50
    sbas_smoothing = config.getfloat('py-config', 'sbas_smoothing') if (
        config.has_option('py-config', 'sbas_smoothing')) else 1
    ts_type = config.get('py-config', 'ts_type')
    ts_format = config.get('py-config', 'ts_format')
    file_format = config.get('py-config', 'file_format')
    intf_dir = config.get('py-config', 'intf_dir')
    intf_filename = config.get('py-config', 'intf_filename')
    corr_filename = config.get('py-config', 'corr_filename')
    signal_spread_filename = config.get('py-config', 'signal_spread_filename')
    make_signal_spread = config.getint('py-config', 'make_signal_spread') if (
        config.has_option('py-config', 'make_signal_spread')) else 1
    signal_coh_cutoff = config.getfloat('py-config', 'signal_coh_cutoff') if (
        config.has_option('py-config', 'signal_coh_cutoff')) else 0.1
    baseline_file = config.get('py-config', 'baseline_file') if config.has_option('py-config', 'baseline_file') else 0
    geocoded_intfs = config.getint('py-config', 'geocoded_intfs') if (
        config.has_option('py-config', 'geocoded_intfs')) else 0
    start_time = config.get('py-config', 'start_time') if config.has_option('py-config', 'start_time') else '19900101'
    end_time = config.get('py-config', 'end_time') if config.has_option('py-config', 'end_time') else '20500101'
    coseismic = config.get('py-config', 'coseismic') if config.has_option('py-config', 'coseismic') else ''
    intf_timespan = config.get('py-config', 'intf_timespan') if config.has_option('py-config', 'intf_timespan') else ''
    gps_file = config.get('py-config', 'gps_file') if config.has_option('py-config', 'gps_file') else ''
    flight_angle = config.getfloat('py-config', 'flight_angle') if config.has_option('py-config', 'flight_angle') else 0
    look_angle = config.getfloat('py-config', 'look_angle') if config.has_option('py-config', 'look_angle') else 0
    skip_file = config.get('py-config', 'skip_file') if config.has_option('py-config', 'skip_file') else ''
    ts_points_file = config.get('py-config', 'ts_points_file') if (
        config.has_option('py-config', 'ts_points_file')) else ''
    ts_output_dir = config.get('py-config', 'ts_output_dir')

    # Start time and end times in datetime format. 
    start_time = dt.datetime.strptime(start_time, "%Y%m%d")
    end_time = dt.datetime.strptime(end_time, "%Y%m%d")
    if coseismic != "":
        coseismic = dt.datetime.strptime(coseismic, "%Y%m%d")

    # enforce startstage <= endstage
    if endstage < startstage:
        print('Warning: endstage is less than startstage. Setting endstage = startstage.')
        endstage = startstage

    config_params = Params(config_file=config_file, SAT=SAT, wavelength=wavelength, startstage=startstage,
                           endstage=endstage,
                           ref_loc=ref_loc, ref_idx=ref_idx, ts_type=ts_type, custom_unwrapping=custom_unwrapping,
                           detrend_atm_topo=detrend_atm_topo, gacos=gacos, aps=aps, dem_error=dem_error,
                           sbas_smoothing=sbas_smoothing, ts_format=ts_format, file_format=file_format,
                           nsbas_min_intfs=nsbas_min_intfs, intf_filename=intf_filename, corr_filename=corr_filename,
                           baseline_file=baseline_file, geocoded_intfs=geocoded_intfs,
                           start_time=start_time, end_time=end_time, coseismic=coseismic, intf_timespan=intf_timespan,
                           gps_file=gps_file, flight_angle=flight_angle, look_angle=look_angle,
                           skip_file=skip_file, signal_spread_filename=signal_spread_filename, 
                           make_signal_spread=make_signal_spread, signal_coh_cutoff=signal_coh_cutoff, 
                           ts_points_file=ts_points_file, intf_dir=intf_dir,
                           ts_output_dir=ts_output_dir)

    return config, config_params


def print_default_config(filename):
    print("Writing %s " % filename)
    ifile = open(filename, 'w')
    ifile.write("# Config file for Stacking methods\n")
    ifile.write("#############################################\n")
    ifile.write("# Python configparser module section header \n")
    ifile.write("[py-config]\n\n")
    ifile.write("# satellite options: S1\n")
    ifile.write("satellite = S1\n")
    ifile.write("wavelength = 56\n\n")
    ifile.write("# stage0: setup\n")
    ifile.write("# stage1: corrections for atmosphere and errors etc.\n")
    ifile.write("# stage2: get reference pixel\n")
    ifile.write("# stage3: velocity and time series formation\n")
    ifile.write("# stage4: geocode\n")
    ifile.write("startstage = 0\n")
    ifile.write("endstage  =  0\n\n")
    ifile.write("# Reference pixel (lon/lat, and/or row/col)\n")
    ifile.write("ref_loc = \n")
    ifile.write("ref_idx = \n\n")
    ifile.write("# timeseries type: STACK or NSBAS or WNSBAS or COSEISMIC\n")
    ifile.write("# timeseries format: velocity, points, timeseries, velocities_from_timeseries\n")
    ifile.write("ts_type = \n")
    ifile.write("ts_format = \n\n")
    ifile.write("# File I/O Options\n")
    ifile.write("# intf_dir should have all the folders with intfs (format YYYYJJJ_YYYYJJJ)\n")
    ifile.write("# geocoded_intfs: are the interferograms already geocoded? \n")
    ifile.write("# ts_output_dir is where all outputs will live\n")
    ifile.write("# file_format is either 'isce' or 'gmtsar'\n")
    ifile.write("intf_dir = \n")
    ifile.write("intf_filename = unwrap.grd\n")
    ifile.write("corr_filename = corr.grd\n")
    ifile.write("file_format = gmtsar\n")
    ifile.write("geocoded_intfs = 1\n")
    ifile.write("ts_output_dir = Output/\n")
    ifile.write("make_signal_spread = 1\n")
    ifile.write("signal_coh_cutoff = 0\n")
    ifile.write("signal_spread_filename = signalspread.nc\n")
    ifile.write("baseline_file = \n\n")
    ifile.write("# sbas parameters\n")
    ifile.write("sbas_smoothing = 1\n\n")
    ifile.write("# nsbas parameters: minimum % of good igrams for nsbas, or -1 for full-rank pixels only\n")
    ifile.write("nsbas_min_intfs = 50\n\n")
    ifile.write("# Do you want to choose a subset of your images to generate a time series? \n")
    ifile.write("# intf_timespan is the duration of interferograms you want to use (300- means less than 300 days "
                "300+ means greater than 300 days)\n")
    ifile.write("start_time = 20141001\n")
    ifile.write("end_time = 20190704\n")
    ifile.write("coseismic = \n")
    ifile.write("intf_timespan = \n\n")
    ifile.write("# choose which interferograms to skip (bad intfs)\n")
    ifile.write("skip_file = \n\n")
    ifile.write("# Choose points to reverse-geocode and get their velocities and time series\n")
    ifile.write("ts_points_file = \n\n")
    ifile.write("#de-trending options, including GPS velocity file, optional to be used for de-trending\n")
    ifile.write("detrend = False\n")
    ifile.write("trendparams = \n")
    ifile.write("flight_angle = 194\n")
    ifile.write("look_angle = 30\n\n")
    ifile.write("# Options for fancy corrections: phase unwrapping errors, and atmospheric corrections\n")
    ifile.write("# DEM error: Fattahi and Amelung, IEEE, 2013. Assumes linear deformation. \n")
    ifile.write("detrend_atm_topo = 0\n")
    ifile.write("dem_error = 0\n")
    ifile.close()
    return
