#!/usr/bin/env python

import sys
import argparse
from s1_batches.read_write_insar_utilities import insar_format_convert

help_message = "Convert isce two-band binary interferogram into single-band phase file. \nUsage: " \
               "isceintf_to_phase.py --isce_file filt_fine.int --outfile phase.int"


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
    p.add_argument('--outfile', type=str, required=True,
                   help='''filename of output binary file, required''')
    config_default = {}
    p.set_defaults(**config_default)
    config = vars(p.parse_args())
    iscefilename = config['isce_file']
    outfilename = config['outfile']
    return iscefilename, outfilename


if __name__ == "__main__":
    isce_filename, outfile = cmd_parser(cmdargs=sys.argv)
    insar_format_convert.convert_intf_phase(isce_filename, outfile)
