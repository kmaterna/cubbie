#!/usr/bin/env python

import sys
import argparse
from cubbie.read_write_insar_utilities import insar_format_convert

help_message = "Convert isce two-band binary unw datafile into single-band phase file. \nUsage: " \
               "unw_to_phase.py --unw_file filt_fine.unw --outfile filt_fine_single.unw"


def cmd_parser(cmdargs):
    """Simple command line parser for converting insar filetypes. Returns two params."""
    p = argparse.ArgumentParser(
          description="\n"+help_message+"\n")
    if len(cmdargs) < 2:
        print("\n"+help_message+"\n")
        print("For full help arguments list, try --help")
        sys.exit(0)
    p.add_argument('--unw_file', type=str, required=True,
                   help='''filename of an isce file, required''')
    p.add_argument('--outfile', type=str, required=True,
                   help='''filename of output binary file, required''')
    config_default = {}
    p.set_defaults(**config_default)
    config = vars(p.parse_args())
    iscefilename = config['unw_file']
    outfilename = config['outfile']
    return iscefilename, outfilename


if __name__ == "__main__":
    isce_filename, outfile = cmd_parser(cmdargs=sys.argv)
    insar_format_convert.extract_unw_phase(isce_filename, outfile)
