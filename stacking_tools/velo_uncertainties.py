import numpy as np
from . import nsbas


def compute_empirical_uncertainties(i, j, ts_tuple, x_axis_days):
    """
    One way of generating uncertainty: 1/2 the RMS of the residuals from a linear fit to the time series.
    This function is passed into an iterator.
    """
    ts_values = ts_tuple.zvalues[:, i, j]
    if sum(np.isnan(ts_values)) > 30:   # under certain degenerate conditions, nanflag=1
        return np.nan, 1, {}
    vel, _, _ = nsbas.compute_velocity_from_ts(i, j, ts_tuple, x_axis_days)   # in mm/yr
    vel_days = vel / 365.24  # vel in mm/day
    # Generate a residual time series.
    ts_model = ts_values[0] + np.multiply(vel_days, x_axis_days)
    ts_residual = np.subtract(ts_values, ts_model)
    unc_empirical = 0.5 * np.sqrt(np.nanmean(np.square(ts_residual)))
    if unc_empirical <= 2:
        unc_empirical = 2
    return unc_empirical, 0, {}


def empirical_uncertainty(ts_tuple):
    """Find empirical uncertainties from a 2D grid of time series"""
    retval_main = np.zeros([len(ts_tuple.yvalues), len(ts_tuple.xvalues)])
    retval_metrics = [[{} for _i in range(len(ts_tuple.xvalues))] for _j in range(len(ts_tuple.yvalues))]
    x_axis_days = [(i - ts_tuple.ts_dates[0]).days for i in ts_tuple.ts_dates]

    def packager_function(i, j, ts_tuple):
        # Giving access to all these variables.
        return compute_empirical_uncertainties(i, j, ts_tuple, x_axis_days)

    retval_main, retval_metrics = nsbas.iterator_func(ts_tuple, packager_function, retval_main, retval_metrics)

    return retval_main
