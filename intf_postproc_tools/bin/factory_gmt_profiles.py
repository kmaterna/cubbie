#!/usr/bin/env python
"""
Convenient wrapper for writing auto-generated GMT code that extracts multiple profiles on the same data file.
It encapsulates the data (lon, lat, width, edges) of the profiles.
Example runstring: factory_gmt_profiles.py --profile_file profile_input_data/profile_specs.txt
--data_file los_igram1.grd --outdir profiles_taken
"""

import numpy as np
import argparse

def arg_parser():
    p = argparse.ArgumentParser()
    p.add_argument('-p', '--profile_file', type=str,
                   help='''filename for specifying the desired profiles, REQUIRED''', required=True)
    p.add_argument('-i', '--data_file', type=str,
                   help='''filename for phase information, grd file, REQUIRED''', required=True)
    p.add_argument('-o', '--outdir', type=str,
                   help='''Output directory for profiles, string, REQUIRED''', required=True)
    exp_dict = vars(p.parse_args())
    return exp_dict


def program_main(exp_dict):
    profile_data = read_specs(exp_dict['profile_file'])
    write_gmt_lines(profile_data, exp_dict['data_file'], exp_dict['outdir'])
    return


def read_specs(infile):
    lons, lats, lengths, widths, azimuths, names = np.loadtxt(infile, unpack=True, dtype=
    {'names': ('lon', 'lat', 'length', 'width', 'azimuth', 'name'),
     'formats': (float, float, float, float, float, 'U20')})
    profiles = []
    for i in range(len(lons)):
        new_profile = (lons[i], lats[i], lengths[i], widths[i], azimuths[i], names[i])
        profiles.append(new_profile)
    return profiles


def write_gmt_lines(profile_data, grdfile, outdir):
    print("# BEGIN AUTO GENERATED CODE")
    for profile in profile_data:
        grdfile_name = grdfile.split('.grd')[0]
        grdfile_specific_name = grdfile_name.split('/')[-1]
        lby2 = profile[2] / 2
        wby2 = profile[3] / 2
        # Example: echo "-115.80 33.003 30 8 0.4" | gmt psxy -R -J -SJ -Wthick,purple -O -K >> $output
        print("echo \"" + str(profile[0]) + " " + str(profile[1]) + " " + str(profile[4]) + " " +
              str(profile[2]) + " " + str(profile[3]) + "\" | gmt psxy -R -J -SJ -Wthick,purple -O -K >> $output")
        print("gmt grd2xyz " + grdfile + " | gmt project -C" + str(profile[0]) + "/" + str(profile[1]) + " -A" + str(
            profile[4]) +
              " -L" + str(-lby2) + "/" + str(lby2) + " -W" + str(-wby2) + "/" + str(wby2) + " -Q > " +
              outdir + "/" + grdfile_specific_name + "_"+profile[5]+".txt")
    print("# END AUTO GENERATED CODE")
    return


if __name__ == "__main__":
    exp_dict = arg_parser()
    program_main(exp_dict)
