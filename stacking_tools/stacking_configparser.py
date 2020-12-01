# Stacking config parser
import argparse, configparser
import datetime as dt
import collections

Params = collections.namedtuple('Params',
                                ['config_file', 'SAT', 'wavelength', 'startstage', 'endstage', 'ref_loc', 'ref_idx',
                                 'ts_type', 'solve_unwrap_errors', 'detrend_atm_topo', 'gacos', 'aps', 'dem_error',
                                 'sbas_smoothing', 'ts_format',
                                 'nsbas_min_intfs', 'intf_filename', 'corr_filename', 'geocoded_intfs', 'baseline_file',
                                 'start_time', 'end_time', 'coseismic', 'intf_timespan', 'gps_file', 'flight_angle',
                                 'look_angle', 'skip_file', 'signal_spread_filename',
                                 'intf_dir', 'ts_points_file', 'ts_output_dir']);


Params_custom = collections.namedtuple('Params_custom', ['config_file', 'rlks', 'alks', 'filt', 'cor_cutoff_mask',
                                                         'xbounds', 'ybounds', 'llh_file', 'lkv_file',
                                                         'ts_output_dir']);


# ----------------------------- #

def read_config_isce(config_file):
    # ISCE-specific parameters

    # read config file
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(config_file)

    rlks = config.getint('py-config', 'rlks')
    alks = config.getint('py-config', 'alks')
    filt = config.getfloat('py-config', 'filt')
    cor_cutoff_mask = config.getfloat('py-config', 'cor_cutoff_mask');
    xbounds = config.get('py-config', 'xbounds');
    ybounds = config.get('py-config', 'ybounds');
    llh_file = config.get('py-config', 'llh_file');
    lkv_file = config.get('py-config', 'lkv_file');
    ts_output_dir = config.get('py-config', 'ts_output_dir');
    ts_output_dir = ts_output_dir + "_" + str(rlks) + "_" + str(alks) + "_" + str(filt);  # custom output directory
    Params = Params_custom(config_file=config_file, rlks=rlks, alks=alks, filt=filt,
                           cor_cutoff_mask=cor_cutoff_mask, xbounds=xbounds, ybounds=ybounds, llh_file=llh_file,
                           lkv_file=lkv_file, ts_output_dir=ts_output_dir);
    return Params;


def read_config_general():
    ################################################
    # Stage 0: Read and check config parameters
    #
    # read command line arguments and parse config file.
    # Common to both GMTSAR and ISCE. 
    parser = argparse.ArgumentParser(description='Run stack processing. ')
    parser.add_argument('config', type=str, help='supply name of config file to setup processing options. Required.')
    parser.add_argument('--debug', action='store_true', help='Print extra debugging messages (default: false)')
    args = parser.parse_args()

    # read config file
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(args.config)

    # get options from config file
    config_file_orig = args.config;
    SAT = config.get('py-config', 'satellite')
    wavelength = config.getfloat('py-config', 'wavelength')
    startstage = config.getint('py-config', 'startstage')
    endstage = config.getint('py-config', 'endstage')
    ref_loc = config.get('py-config', 'ref_loc')
    ref_idx = config.get('py-config', 'ref_idx')
    ts_type = config.get('py-config', 'ts_type')
    solve_unwrap_errors = config.getint('py-config', 'solve_unwrap_errors');
    gacos = config.getint('py-config', 'gacos');
    aps = config.getint('py-config', 'aps');
    dem_error = config.getint('py-config', 'dem_error');
    detrend_atm_topo = config.getint('py-config', 'detrend_atm_topo');
    nsbas_min_intfs = config.getint('py-config', 'nsbas_min_intfs');
    sbas_smoothing = config.getfloat('py-config', 'sbas_smoothing');
    ts_format = config.get('py-config', 'ts_format');
    intf_dir = config.get('py-config', 'intf_dir');
    intf_filename = config.get('py-config', 'intf_filename');
    corr_filename = config.get('py-config', 'corr_filename');
    signal_spread_filename = config.get('py-config', 'signal_spread_filename');
    baseline_file = config.get('py-config', 'baseline_file');
    geocoded_intfs = config.getint('py-config', 'geocoded_intfs');
    start_time = config.get('py-config', 'start_time');
    end_time = config.get('py-config', 'end_time');
    coseismic = config.get('py-config', 'coseismic');
    intf_timespan = config.get('py-config', 'intf_timespan');
    gps_file = config.get('py-config', 'gps_file');
    flight_angle = config.getfloat('py-config', 'flight_angle');
    look_angle = config.getfloat('py-config', 'look_angle');
    skip_file = config.get('py-config', 'skip_file');
    ts_points_file = config.get('py-config', 'ts_points_file');
    ts_output_dir = config.get('py-config', 'ts_output_dir');

    # Start time and end times in datetime format. 
    start_time = dt.datetime.strptime(start_time, "%Y%m%d");
    end_time = dt.datetime.strptime(end_time, "%Y%m%d");
    if coseismic != "":
        coseismic = dt.datetime.strptime(coseismic, "%Y%m%d");

    print("Running velocity and time series processing, starting with stage %d" % startstage);

    # enforce startstage <= endstage
    if endstage < startstage:
        print('Warning: endstage is less than startstage. Setting endstage = startstage.')
        endstage = startstage

    config_params = Params(config_file=config_file_orig, SAT=SAT, wavelength=wavelength, startstage=startstage,
                           endstage=endstage,
                           ref_loc=ref_loc, ref_idx=ref_idx, ts_type=ts_type, solve_unwrap_errors=solve_unwrap_errors,
                           detrend_atm_topo=detrend_atm_topo, gacos=gacos, aps=aps, dem_error=dem_error,
                           sbas_smoothing=sbas_smoothing, ts_format=ts_format,
                           nsbas_min_intfs=nsbas_min_intfs, intf_filename=intf_filename, corr_filename=corr_filename,
                           baseline_file=baseline_file, geocoded_intfs=geocoded_intfs,
                           start_time=start_time, end_time=end_time, coseismic=coseismic, intf_timespan=intf_timespan,
                           gps_file=gps_file, flight_angle=flight_angle, look_angle=look_angle,
                           skip_file=skip_file, signal_spread_filename=signal_spread_filename, 
                           ts_points_file=ts_points_file, intf_dir=intf_dir,
                           ts_output_dir=ts_output_dir);

    return config, config_params;
