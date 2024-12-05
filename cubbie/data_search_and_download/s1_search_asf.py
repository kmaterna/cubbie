#!/usr/bin/env python

"""
This script queries the ASF for Sentinel-1 SLC TOPS-mode images.

As of October 2023, searching Sentinel-1 data has migrated away from SciHub
"""

import asf_search as asf
import datetime as dt
import pygmt
import matplotlib.pyplot as plt
import argparse
import sys
import os
import numpy as np
from argparse import RawTextHelpFormatter

# Help message
help_message = "\n"\
  "This script queries ASF for Sentinel-1 SLC TOPS-mode images.\n"\
  "Remember to surround numerical values in quotations and possibly with a space "\
               "if they begin with negative numbers (such as longitudes).\n\n\n"\
  "Usage: s1_search_asf.py -options\n"\
  "Example: s1_search_asf.py -s 2015-08-01 -e NOW -r \"-123.0/-123.3/40.0/40.2\" -d Descending\n\n" \
  "Example: s1_search_asf.py -s 2015-08-01 -e 2017-01-01 -c \" -115.0/32\" -d Descending\n\n" \
  "\n"


def cmd_parse():
    p = argparse.ArgumentParser(description=help_message, formatter_class=RawTextHelpFormatter)
    p.add_argument('-s', '--start_time', type=str, help='''Starting date [yyyy-mm-dd]''')
    p.add_argument('-e', '--end_time', type=str, help='''Ending date [yyyy-mm-dd]''', default='NOW')
    p.add_argument('-b', '--bbox', type=str, help='''Bounding box [lonW/lonE/latS/latN]''')
    p.add_argument('-r', '--region', type=str, help='''Four-sided polygon [lon1/lat1/lon2/lat2/lon3/lat3/lon4/lat4]''')
    p.add_argument('-c', '--coordinate', type=str, help='''A single point [lon/lat]''')
    p.add_argument('-o', '--orbit_number', type=int, help='''Track or relative orbit [0-175]''')
    p.add_argument('-d', '--orbit_direction', type=str, help='''[Ascending/Descending]''', default='both')
    p.add_argument('-m', '--sar_mode', type=str, help='''[IW/other/ALL]''', default='IW')
    p.add_argument('-z', '--output_file', type=str, help='''Specified output file''', default="search_results.txt")
    p.add_argument('-v', '--write_footprints', help='''Write the footprints in GMT-friendly format''', action='store_true')
    exp_dict = vars(p.parse_args(args=None if sys.argv[1:] else ['--help']))
    return exp_dict


def asf_query(args):
    # This is a basic search

    # Defensive programming for basic inputs
    if args['end_time'] == "NOW":
        args['end_time'] = dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d")
    if args['start_time'] is None:
        args['start_time'] = '2014-01-01'
    if args['bbox'] is not None:
        bbox_list = args['bbox'].split('/')
        args['bbox'] = [float(bbox_list[0]), float(bbox_list[1]), float(bbox_list[2]), float(bbox_list[3])]
    if args['coordinate'] is not None:
        coords_list = args['coordinate'].split('/')
        args['coordinate'] = [float(coords_list[0]), float(coords_list[1])]

    opts = {
        'platform': asf.PLATFORM.SENTINEL1,
        'start': args['start_time'] + 'T00:00:00Z',
        'end': args['end_time'] + 'T23:59:59Z',
        'processingLevel': 'SLC'
    }

    # Parse the orbit direction
    if args['orbit_direction'] == 'Ascending':
        opts['flightDirection'] = 'ASCENDING'
    elif args['orbit_direction'] == 'Descending':
        opts['flightDirection'] = 'DESCENDING'
    elif args['orbit_direction'] == 'both':
        print("Searching both orbit directions.")
    else:
        print("Provided orbit direction is invalid: %s " % args['orbit_direction'])
        sys.exit(0)
    if args['orbit_number'] is not None:
        opts['relativeOrbit'] = args['orbit_number']
    if args['sar_mode'] != 'ALL':
        opts['beamMode'] = args['sar_mode']

    # Searching for bounding box
    if args['bbox'] is not None:
        aoi = ('POLYGON((' + str(args['bbox'][0]) + ' ' + str(args['bbox'][2]) + ','
                           + str(args['bbox'][0]) + ' ' + str(args['bbox'][3]) + ','
                           + str(args['bbox'][1]) + ' ' + str(args['bbox'][3]) + ','
                           + str(args['bbox'][1]) + ' ' + str(args['bbox'][2]) + ','
                           + str(args['bbox'][0]) + ' ' + str(args['bbox'][2]) + '))')
        results = asf.geo_search(intersectsWith=aoi, **opts)
    # Searching for point
    elif args['coordinate'] is not None:
        point = 'POINT('+str(args['coordinate'][0])+' '+str(args['coordinate'][1])+')'
        results = asf.geo_search(intersectsWith=point, **opts)
    # Searching for polygon
    elif args['region'] is not None:
        print("Searching for arbitrary polygon not yet supported. Exiting.")
        aoi = 'POLYGON((-152.81 58.49,-154.90 57.49,-155.08 56.30,-153.82 56.34,-151.99 57.30,-152.81 58.49))'
        results = asf.geo_search(intersectsWith=aoi, **opts)
    else:
        print("Either a point or a bounding box must be set. Exiting.")
        sys.exit(0)

    for item in results:
        print(item.properties['fileID'])

    total_dates = [item.properties['startTime'].split('T')[0] for item in results]
    print("Total number of results:", len(results))
    print("On total number of dates:", len(set(total_dates)))
    if len(results) == 0:
        print("Exiting.")
        sys.exit(0)

    test = results.jsonlite()
    with open('raw_search_results.txt', 'w') as f:
        f.writelines(list(test))

    with open(args['output_file'], 'w') as f:
        f.write("# Search Terms: "+str(args)+"\n")
        for item in results:
            full_name = item.properties["fileID"]
            acq_date = item.properties['startTime'].split('T')[0]
            acq_time = item.properties['startTime'].split('T')[-1]
            direction = item.properties['flightDirection']
            track = item.properties['pathNumber']
            slc_name = full_name.split('/')[-1]
            slc_name = slc_name.replace('-SLC', '')
            f.write(slc_name+' '+acq_date+' '+acq_time+' '+direction+' '+str(track) + '\n')
    print("Writing output file %s " % args['output_file'])

    # Pretty plots
    pygmt_plots(results, args)
    timing_plots(args['output_file'], "timing.png")
    return


def get_wesn_from_coordinates(returned_object):
    busy_array = returned_object.geometry['coordinates'][0]
    lons = [x[0] for x in busy_array]
    lats = [x[1] for x in busy_array]
    wesn = [np.min(lons), np.max(lons), np.min(lats), np.max(lats)]
    return wesn


def get_general_bbox(results):
    """Returns in WESN"""
    overall_bbox = get_wesn_from_coordinates(results[0])
    for item in results:
        new_bbox = get_wesn_from_coordinates(item)
        if new_bbox[0] < overall_bbox[0]:
            overall_bbox[0] = new_bbox[0]
        if new_bbox[1] > overall_bbox[1]:
            overall_bbox[1] = new_bbox[1]
        if new_bbox[2] < overall_bbox[2]:
            overall_bbox[2] = new_bbox[2]
        if new_bbox[3] > overall_bbox[3]:
            overall_bbox[3] = new_bbox[3]
    return overall_bbox[0], overall_bbox[1], overall_bbox[2], overall_bbox[3]


def get_bbox_drawing_points(bbox):
    lons = [bbox[0], bbox[0], bbox[1], bbox[1], bbox[0]]
    lats = [bbox[2], bbox[3], bbox[3], bbox[2], bbox[2]]
    return lons, lats


def pygmt_plots(results, args):
    region = get_general_bbox(results)
    proj = "M6i"
    fig = pygmt.Figure()
    title = "+t\"Search results: "+str(len(results))+" acquisitions\""  # must put escaped quotations around title.
    fig.basemap(region=region, projection=proj, frame=title)
    fig.coast(shorelines="1.0p,black", region=region, borders="1", projection=proj, frame="1.0")  # the boundary.
    fig.coast(region=region, projection=proj, borders='2', shorelines='0.5p,black', water='white')

    # Write the footprint polygons
    if args['write_footprints']:
        footprint_file = "footprints.txt"
    else:
        footprint_file = "tmp.txt"
    with open(footprint_file, 'w') as f:
        for im in results:
            len_coords = len(im.geometry['coordinates'][0])
            f.write(">\n")
            for i in range(len_coords):
                f.write("%f %f\n" % (im.geometry['coordinates'][0][i][0], im.geometry['coordinates'][0][i][1]))
    fig.plot(data=footprint_file, pen="0.4p,red")
    if not args['write_footprints']:
        os.remove(footprint_file)

    # Write the search coordinate or search region
    if args['coordinate'] is not None:
        fig.plot(x=args['coordinate'][0], y=args['coordinate'][1], pen="0.4p,black", style='c0.2c', fill='black')
    if args['bbox'] is not None:
        bbox_lons, bbox_lats = get_bbox_drawing_points(args['bbox'])
        fig.plot(x=bbox_lons, y=bbox_lats, pen="0.4p,black,dashed")

    print("Mapping results in %s " % "footprints.png")
    fig.savefig("footprints.png")
    return


def timing_plots(results_file, output_plot='timing.png'):
    """
    Build a plot of the acquisition timing for a set of ASF search results. Works from a text file directly now.

    :param results_file: string, name of a text file with human-readable results from an ASF query
    :param output_plot: string, name of the output png
    """
    borders = [dt.datetime.strptime("2014-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2017-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2020-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2023-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2026-07-01", "%Y-%m-%d")]

    # Divide the results into sub-plots
    a1, d1, a2, d2, a3, d3, a4, d4 = [], [], [], [], [], [], [], []

    # Read the data in from the file
    datestrs, flight_directions, tracks = np.loadtxt(results_file, usecols=(1, 3, 4), unpack=True,
                                                     dtype={'names': ('dts', 'direction', 'track'),
                                                            'formats': ('U10', 'U9', float)})

    for datestr, direction, track in zip(datestrs, flight_directions, tracks):
        acq_date = dt.datetime.strptime(datestr, "%Y-%m-%d")
        if borders[0] < acq_date < borders[1]:
            if direction == "ASCENDING":
                a1.append(acq_date)
            else:
                d1.append(acq_date)
        if borders[1] < acq_date < borders[2]:
            if direction == "ASCENDING":
                a2.append(acq_date)
            else:
                d2.append(acq_date)
        if borders[2] < acq_date < borders[3]:
            if direction == "ASCENDING":
                a3.append(acq_date)
            else:
                d3.append(acq_date)
        if borders[3] < acq_date < borders[4]:
            if direction == "ASCENDING":
                a4.append(acq_date)
            else:
                d4.append(acq_date)

    fig, axarr = plt.subplots(4, 1, figsize=(14, 10), dpi=300)
    ms = 4

    axarr[0].set_title("Search Results: %s acquisitions" % (len(datestrs)), fontsize=20)
    axarr[0].set_xlim([borders[0], borders[1]])
    axarr[0].set_ylim([0, 1])
    axarr[0].set_yticks([])
    axarr[0].plot(a1, [0.3 for _x in a1], color='red', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[0].plot(d1, [0.7 for _x in d1], color='blue', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[0].grid(True)

    axarr[1].set_xlim([borders[1], borders[2]])
    axarr[1].set_ylim([0, 1])
    axarr[1].set_yticks([])
    axarr[1].plot(a2, [0.3 for _x in a2], color='red', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[1].plot(d2, [0.7 for _x in d2], color='blue', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[1].grid(True)

    axarr[2].set_xlim([borders[2], borders[3]])
    axarr[2].set_ylim([0, 1])
    axarr[2].set_yticks([])
    axarr[2].plot(a3, [0.3 for _x in a3], color='red', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[2].plot(d3, [0.7 for _x in d3], color='blue', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[2].grid(True)

    axarr[3].set_xlim([borders[3], borders[4]])
    axarr[3].set_ylim([0, 1])
    axarr[3].set_yticks([])
    axarr[3].plot(a4, [0.3 for _x in a4], color='red', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[3].plot(d4, [0.7 for _x in d4], color='blue', marker='o', markersize=ms, linestyle=None, linewidth=0)
    axarr[3].grid(True)
    axarr[3].legend(['Ascending', 'Descending'], fontsize=15)
    print("Plotting acquisitions in "+output_plot)
    fig.savefig(output_plot)
    return


if __name__ == "__main__":
    my_args = cmd_parse()
    asf_query(my_args)
