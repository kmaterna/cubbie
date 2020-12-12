# This is a function to write and call the GMTSAR SBAS control files
# Because it only does the GMTSAR standard options, I think it's better
# to use nsbas.py in almost all cases. 
# The staging_directory functionality doesn't work right now because this authomatically reads intf_all/*/unwrap.grd

import numpy as np
from intf_generating import sentinel_utilities
from subprocess import call


def do_sbas(config_params, staging_directory):
    [stems, tbaseline, xbaseline, mission_days] = sentinel_utilities.read_baseline_table('raw/baseline_table.dat');
    t_int = [];
    for t in tbaseline:
        t_int.append(round(float(t)));  # a list of integers like 2016214 for 2016-day-214.

    mission_days_sorted = [x for (y, x) in sorted(zip(t_int, mission_days))];
    t_int.sort();
    tbaseline.sort();

    intf_computed = sentinel_utilities.glob_intf_computed();  # looks like a list of labels like 2016217_2016205
    n_intf = len(intf_computed);
    outfile = open("README_sbas.txt", 'w');
    outfile.write("# First, prepare the input files needed for sbas\n#\n");
    outfile.write("rm -f SBAS\nmkdir SBAS\ncd SBAS\n\n\n");
    outfile.write("# based on baseline_table.dat create the intf.tab and scene.tab for sbas\n");

    # writing intf.tab
    outfile.write("# phase  corherence  ref_id  rep_id  baseline\n")
    for img_pair in intf_computed:
        first_image = img_pair[0:7]
        second_image = img_pair[8:]
        for a, b in zip(xbaseline, t_int):
            if abs(int(np.floor(b)) - int(first_image)) <= 1:
                # print "first image found";
                # print int(np.floor(b))
                # print int(first_image)
                master_xbaseline = a;
            if abs(int(np.floor(b)) - int(second_image)) <= 1:
                slave_xbaseline = a;
                # print "second image found";
                # print int(np.floor(b))
                # print int(second_image)                
        total_baseline = slave_xbaseline - master_xbaseline;
        outfile.write('echo "../intf_all/' + img_pair + '/unwrap.grd ')
        outfile.write('../intf_all/' + img_pair + '/corr.grd ')
        outfile.write(first_image + ' ' + second_image + ' ')
        outfile.write(str(total_baseline))
        outfile.write('" >> intf.tab\n');
    outfile.write("#\n\n");

    # writing scene.tab (only the scenes that are actually used in SBAS processing)
    # Right now this isn't producing anything. 
    outfile.write("# scene_id  day\n");
    scenes_used = '';
    n_scenes = 0;
    scenes_used = [];
    for intf in intf_computed:
        scenes_used.append(intf[0:7]);  # catch which scenes are actually used in SBAS processing. 
        scenes_used.append(intf[8:15]);
    for x in range(len(tbaseline)):
        temp = tbaseline[x];
        tempint = int(np.round(temp))
        if str(tempint) in scenes_used or str(tempint + 1) in scenes_used or str(tempint - 1) in scenes_used:
            outfile.write('echo "' + str(tempint) + ' ' + mission_days_sorted[x] + '" >> scene.tab\n');
            n_scenes += 1;
    outfile.write("\n\n");

    intf_ex = intf_computed[0];  # an example interferogram where we get the geographic coordinates for grdinfo
    outfile.write("xdim=`gmt grdinfo -C ../intf_all/" + intf_ex + "/unwrap.grd | awk '{print $10}'`\n");
    outfile.write("ydim=`gmt grdinfo -C ../intf_all/" + intf_ex + "/unwrap.grd | awk '{print $11}'`\n\n\n");

    outfile.write("# run sbas\n");
    outfile.write("sbas intf.tab scene.tab " + str(n_intf) + " " + str(
        n_scenes) + " $xdim $ydim -smooth 1.0 -wavelength 0.0554658 -incidence 30 -range 800184.946186 -rms -dem\n\n\n")

    outfile.write("# project the velocity to Geocooridnates\n");
    outfile.write('echo "writing to georeferenced coordinates..."\n');
    outfile.write("ln -s ../topo/trans.dat .\n");
    outfile.write("proj_ra2ll.csh trans.dat vel.grd vel_ll.grd\n");
    outfile.write("gmt grd2cpt vel_ll.grd -T= -Z -Cjet > vel_ll.cpt\n");
    outfile.write("grd2kml.csh vel_ll vel_ll.cpt\n");
    outfile.write("cd ..\n\n");
    outfile.write('echo "SBAS operation performed!"\n\n')

    outfile.close();
    print("README_sbas.txt written. Ready to call README_sbas.txt.")
    call("chmod u+x README_sbas.txt", shell=True);
    call("./README_sbas.txt", shell=True);  # Make sbas!
    return;
