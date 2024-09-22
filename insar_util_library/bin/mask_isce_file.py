#!/usr/bin/env python

"""Apply a mask to an isce file"""

import argparse
from s1_batches.insar_util_library.read_write_insar_utilities import isce_read_write
from s1_batches.insar_util_library.math_tools import mask_and_interpolate as mask_tools
import numpy as np

help_message = "Produce a masked version of an ISCE binary data file (single band file only). \nUsage: " \
               "mask_isce_file.py --data_file phase.int --outfile phase_masked.int " \
                "--coherence filt_fine.cor --cutoff 0.37"


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', '--data_file', type=str,
                   help='''filename for raster information, REQUIRED.''', required=True)
    p.add_argument('-o', '--outfile', type=str,
                   help='''filename for output data file, REQUIRED.''', required=True)
    p.add_argument('-c', '--maskfile', type=str,
                   help='''filename for coherence data file, REQUIRED.''', required=True)
    p.add_argument('-t', '--threshold', type=float,
                   help='''Coherence threshold for masking (0-1).''', default=0.37)
    p.add_argument('-v', '--mask_value', type=float,
                   help='''Value used for masked data [default nan].''', default=np.nan)
    exp_dict = vars(p.parse_args())
    return exp_dict


def main_body(paramdict):
    x, y, scalars = isce_read_write.read_scalar_data(paramdict['data_file'], band=1)
    _, _, cor = isce_read_write.read_scalar_data(paramdict['maskfile'], band=1)
    mask = mask_tools.make_coherence_mask(cor, paramdict['threshold'])
    masked = mask_tools.apply_coherence_mask(scalars, mask, is_float32=True, mask_value=paramdict['mask_value'])
    ny, nx = np.shape(masked)
    isce_read_write.write_isce_data(masked, nx, ny, "FLOAT", paramdict['outfile'])
    return


if __name__ == "__main__":
    args_dict = parse_arguments()
    main_body(args_dict)
