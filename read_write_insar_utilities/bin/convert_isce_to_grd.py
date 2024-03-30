#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Driver program to convert an isce interferogram to gmt GRD file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import sys
import argparse
from s1_batches.read_write_insar_utilities import insar_format_convert

help_message = "Convert geocoded isce interferogram into gmt grdfile. \nUsage: " \
               "convert_isce_to_grd --isce_file geo_fine_filt.int --grd_file phase_ll.grd"


def cmd_parser(cmdargs):
    """Simple command line parser for converting insar filetypes. Returns two params."""
    p = argparse.ArgumentParser(
          description="\n"+help_message+"\n")
    if len(cmdargs) < 2:
        print("\n"+help_message+"\n")
        print("For full help arguments list, try --help")
        sys.exit(0)
    p.add_argument('--isce_file', type=str, required=True,
                   help='''filename of an isce file, required''')
    p.add_argument('--grd_file', type=str, required=True,
                   help='''filename of output grd file, required''')
    config_default = {}
    p.set_defaults(**config_default)
    config = vars(p.parse_args())
    isce_filename = config['isce_file']
    grd_filename = config['grd_file']
    return isce_filename, grd_filename


if __name__ == "__main__":
    isce_filename, grd_filename = cmd_parser(cmdargs=sys.argv)
    insar_format_convert.isce_to_grd(isce_filename, grd_filename)
