#!/usr/bin/env python

"""
This script queries the Odata Sentinel-hub for Sentinel-1 SLC TOPS-mode images.

As of October 2023, searching Sentinel-1 data has migrated from SciHub to Odata
https://dataspace.copernicus.eu/news/2023-9-28-accessing-sentinel-mission-data-new-copernicus-data-space-ecosystem-apis

"""

# Utilities
import datetime as dt
import pygmt
import matplotlib.pyplot as plt
import json
import argparse
import sys
import os
from argparse import RawTextHelpFormatter
from sentinelhub import SHConfig, DataCollection, SentinelHubCatalog, BBox, CRS


# Help message
help_message = "\n"\
  "This script queries the copernicus scihub for Sentinel-1 SLC TOPS-mode images.\n"\
  "As of October 2023, Sentinel-1 data has migrated from SciHub to Odata.\n"\
  "Note: To use OData, you will generally need a login and a token (get from https://dataspace.copernicus.eu/ and \n" \
               "https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Overview/Authentication.html).\n"\
  "    Put your credentials into a file called ~/.wget_cred with format:\n"\
  "    Odata_client_id: your_client_id\n"\
  "    Odata_client_secret: your_client_secret\n\n\n"\
  "Remember to surround numerical values in quotations and possibly with a space "\
               "if they begin with negative numbers (such as longitudes).\n\n\n"\
  "Usage: s1_search_Odata.py -options\n"\
  "Example: s1_search_Odata.py -s 2015-08-01 -e NOW -r \"-123.0/-123.3/40.0/40.2\" -d Descending\n\n" \
  "Example: s1_search_Odata.py -s 2015-08-01 -e 2017-01-01 -c \" -115.0/32\" -d Descending\n\n" \
  "\n"


def cmd_parse():
    p = argparse.ArgumentParser(description=help_message, formatter_class=RawTextHelpFormatter)
    p.add_argument('-s', '--start_time', type=str, help='''Starting date [yyyy-mm-dd]''')
    p.add_argument('-e', '--end_time', type=str, help='''Ending date [yyyy-mm-dd]''', default='NOW')
    p.add_argument('-b', '--bbox', type=str, help='''Bounding box [lonW/lonE/latS/latN]''')
    p.add_argument('-r', '--region', type=str, help='''Four-sided polygon [lon1/lat1/lon2/lat2/lon3/lat3/lon4/lat4]''')
    p.add_argument('-c', '--coordinate', type=str, help='''A single point [lon/lat]''')
    p.add_argument('-o', '--orbit_number', type=int, help='''Track or relative orbit [0-175]''')
    p.add_argument('-d', '--orbit_direction', type=str, help='''[Ascending/Descending]''')
    p.add_argument('-m', '--sar_mode', type=str, help='''[IW/other/ALL]''', default='IW')
    p.add_argument('-z', '--output_file', type=str, help='''Specified output file''', default="search_results.txt")
    exp_dict = vars(p.parse_args(args=None if sys.argv[1:] else ['--help']))
    return exp_dict


def SentinelHub_query(args):
    # This is a basic search
    print(args)
    config = SHConfig()
    config.sh_client_id, config.sh_client_secret = get_auth_token_and_pw()
    config.sh_token_url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    config.sh_base_url = "https://sh.dataspace.copernicus.eu"
    config.save("cdse")  # Saved config can be later accessed with config = SHConfig("cdse")

    config = SHConfig("cdse")
    catalog = SentinelHubCatalog(config=config)
    search_iterator = []

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

    time_interval = args['start_time'], args['end_time']

    # Parse the orbit direction
    if args['orbit_direction'] == 'Ascending':
        filter_str = "sat:orbit_state='ascending'"
    elif args['orbit_direction'] == 'Descending':
        filter_str = "sat:orbit_state='descending'"
    elif args['orbit_direction'] is None:
        filter_str = ""
    else:
        print("Provided orbit direction is invalid: %s " % args['orbit_direction'])
        sys.exit(0)

    # Searching for bounding box
    if args['bbox'] is not None:
        bbox_coords = (args['bbox'][0], args['bbox'][2], args['bbox'][1], args['bbox'][3])
        aoi_bbox = BBox(bbox=bbox_coords, crs=CRS.WGS84)
        search_iterator = catalog.search(
            DataCollection.SENTINEL1,
            bbox=aoi_bbox,
            time=time_interval,
            filter=filter_str
        )

    # Searching for point
    elif args['coordinate'] is not None:
        point_dict = {
            "type": "Point",
            "coordinates": [
                args['coordinate'][0],
                args['coordinate'][1],
            ]}
        search_iterator = catalog.search(
            DataCollection.SENTINEL1,
            intersects=point_dict,
            time=time_interval,
            filter=filter_str
        )

    # Searching for polygon
    elif args['coordinate'] is not None:
        print("Searching for arbitrary polygon not yet supported. Exiting.")
    else:
        print("Either a point or a bounding box must be set. Exiting.")
        sys.exit(0)

    results = list(search_iterator)

    # LAST FILTERS: Filter on mode and relative orbit, which aren't provided directly by the API
    if args['orbit_number'] is not None:
        results = filter_by_relative_orbit(results, args['orbit_number'])
    if args['sar_mode'] != 'ALL':
        results = filter_by_sar_mode(results, args['sar_mode'])

    print("Total number of results:", len(results))
    if len(results) == 0:
        print("Exiting.")
        sys.exit(0)

    with open('raw_search_results.txt', 'w') as f:
        json.dump(results, f, indent=2)

    with open(args['output_file'], 'w') as f:
        f.write("# Search Terms: "+str(args)+"\n")
        for item in results:
            full_name = item['assets']['data']['href']
            acq_date = item['properties']['datetime'].split('T')[0]
            acq_time = item['properties']['datetime'].split('T')[-1]
            direction = item['properties']['sat:orbit_state']
            track = item['properties']['sat:relative_orbit']
            slc_name = full_name.split('/')[-1]
            slc_name = slc_name.replace('_COG.SAFE', '')
            f.write(slc_name+' '+acq_date+' '+acq_time+' '+direction+' '+str(track) + '\n')
    print("Writing output file %s " % args['output_file'])

    # Pretty plots
    pygmt_plots(results, args)
    timing_plots(results)
    return


def filter_by_relative_orbit(search_results, track_number):
    track_results = [item for item in search_results if int(item['properties']['sat:relative_orbit']) == track_number]
    return track_results


def filter_by_sar_mode(search_results, mode='IW'):
    track_results = [item for item in search_results if item['properties']['sar:instrument_mode'] == mode]
    return track_results


def get_auth_token_and_pw():
    token, pw = '', ''
    with open(os.path.expanduser('~/.wget_cred'), 'r') as f:
        for line in f:
            if 'Odata_client_id' in line:
                token = line.split()[1]
            if 'Odata_client_secret' in line:
                pw = line.split()[1]
    return token, pw


def get_general_bbox(results):
    """Returns in WESN"""
    overall_bbox = results[0]['bbox']  # WSEN
    for item in results:
        new_bbox = item['bbox']
        if new_bbox[0] < overall_bbox[0]:
            overall_bbox[0] = new_bbox[0]
        if new_bbox[1] < overall_bbox[1]:
            overall_bbox[1] = new_bbox[1]
        if new_bbox[2] > overall_bbox[2]:
            overall_bbox[2] = new_bbox[2]
        if new_bbox[3] > overall_bbox[3]:
            overall_bbox[3] = new_bbox[3]
    return overall_bbox[0], overall_bbox[2], overall_bbox[1], overall_bbox[3]


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
    with open("tmp.txt", 'w') as f:
        for im in results:
            len_coords = len(im['geometry']['coordinates'][0])
            f.write(">\n")
            for i in range(len_coords):
                f.write("%f %f\n" % (im['geometry']['coordinates'][0][i][0], im['geometry']['coordinates'][0][i][1]))
    fig.plot(data='tmp.txt', pen="0.4p,red")
    # os.remove('tmp.txt')

    # Write the search coordinate or search region
    if args['coordinate'] is not None:
        fig.plot(x=args['coordinate'][0], y=args['coordinate'][1], pen="0.4p,black", style='c0.2c', fill='black')
    if args['bbox'] is not None:
        bbox_lons, bbox_lats = get_bbox_drawing_points(args['bbox'])
        fig.plot(x=bbox_lons, y=bbox_lats, pen="0.4p,black,dashed")

    print("Mapping results in %s " % "footprints.png")
    fig.savefig("footprints.png")
    return


def timing_plots(results):
    borders = [dt.datetime.strptime("2014-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2017-01-01", "%Y-%m-%d"),
               dt.datetime.strptime("2019-07-01", "%Y-%m-%d"),
               dt.datetime.strptime("2022-01-01", "%Y-%m-%d"),
               dt.datetime.strptime("2024-07-01", "%Y-%m-%d")]

    # Divide the results into sub-plots
    a1, d1, a2, d2, a3, d3, a4, d4 = [], [], [], [], [], [], [], []
    for item in results:
        acq_date = dt.datetime.strptime(item['properties']['datetime'].split('T')[0], "%Y-%m-%d")
        direction = item['properties']['sat:orbit_state']
        _track = item['properties']['sat:relative_orbit']
        if borders[0] < acq_date < borders[1]:
            if direction == "ascending":
                a1.append(acq_date)
            else:
                d1.append(acq_date)
        if borders[1] < acq_date < borders[2]:
            if direction == "ascending":
                a2.append(acq_date)
            else:
                d2.append(acq_date)
        if borders[2] < acq_date < borders[3]:
            if direction == "ascending":
                a3.append(acq_date)
            else:
                d3.append(acq_date)
        if borders[3] < acq_date < borders[4]:
            if direction == "ascending":
                a4.append(acq_date)
            else:
                d4.append(acq_date)

    fig, axarr = plt.subplots(4, 1, figsize=(12, 10), dpi=300)
    ms = 4

    axarr[0].set_title("Search Results: %s acquisitions" % (len(results)), fontsize=20)
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
    print("Plotting acquisitions in timing.png")
    fig.savefig("timing.png")
    return


if __name__ == "__main__":
    my_args = cmd_parse()
    SentinelHub_query(my_args)
