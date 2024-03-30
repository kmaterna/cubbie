#! /usr/bin/env python
###############################################################################
#  auig2_download.py
#
#  Purpose:  Command line download from AUIG2
#  Author:   Scott Baker
#  Created:  Apr 2015
#
###############################################################################
#  Copyright (c) 2015, Scott Baker
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################

import os
import sys
import time
import datetime
import urllib
import urllib2
import cookielib
import argparse


BASE_URL = 'https://auig2.jaxa.jp/openam/UI/Login'
USERNAME = 'P1298002'  # YOUR USERNAME CAN ALOS BE HARDWIRED HERE
PASSWORD = 'Cba@8580247802645520813'  # YOUR PASSWORD CAN ALOS BE HARDWIRED HERE


def loginToAUIG2(opener, inps):
    """
    Handle login. This should populate our cookie jar.
    """
    login_data = urllib.urlencode({
        'IDToken1': inps.username,
        'IDToken2': inps.password,
    })
    response = opener.open(BASE_URL, login_data)
    return ''.join(response.readlines())


def parse():
    """Command line parser."""
    desc = """Command line client for downloading from AUIG2
For questions or comments, contact Scott Baker: baker@unavco.org
    """
    epi = """You can hardwire your AUIG2 USERNAME and PASSWORD in this file (it's near the top), or use command line args"""
    usage = """Example:
auig2_download.py -o ORDER_ID -u USERNAME -p PASSWORD

If you have your credentials hardwired in this file, just do:
auig2_download.py -o ORDER_ID
"""
    parser = argparse.ArgumentParser(description=desc, epilog=epi, usage=usage)
    parser.add_argument('-o', '--orderid', action="store", dest="order_id", metavar='<ORDERID>', required=True, help='This is your AUIG2 Order ID')
    parser.add_argument('-u', '--username', action="store", dest="username", metavar='<USERNAME>', default=USERNAME, help='AUIG2 Login')
    parser.add_argument('-p', '--password', action="store", dest="password", metavar='<PASSWORD>', default=PASSWORD, help='AUIG2 Login')
    inps = parser.parse_args()
    return inps


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    # ## READ IN PARAMETERS FROM THE COMMAND LINE ###
    inps = parse()

    # ## OPEN A CONNECTION TO AUIG2 AND LOG IN ###
    cookie_filename = "cookie_%s.txt" % datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")
    cj = cookielib.MozillaCookieJar(cookie_filename)
    if os.access(cookie_filename, os.F_OK):
        cj.load()
    opener = urllib2.build_opener(
        urllib2.HTTPRedirectHandler(),
        urllib2.HTTPHandler(debuglevel=0),
        urllib2.HTTPSHandler(debuglevel=0),
        urllib2.HTTPCookieProcessor(cj)
    )
    opener.addheaders = [ ('User-agent', ('Mozilla/4.0 (compatible; MSIE 6.0; ' 'Windows NT 5.2; .NET CLR 1.1.4322)')) ]
    # need this twice - once to set cookies, once to log in...
    loginToAUIG2(opener, inps)
    loginToAUIG2(opener, inps)

    # ## DOWNLOAD THE FILE WITH THE GIVEN ORDER ID ###
    url = "http://auig2.jaxa.jp/pp/service/download?downloadurl=/start/download/file&itemname=%s&itemtype=1" % inps.order_id
    f = opener.open(url)
    filename = f.headers['Content-Disposition'].split("=")[-1].strip()
    print("ALOS-2 AUIG2 Download:", filename)
    start = time.time()
    CHUNK = 256 * 1024
    with open(filename, 'wb') as fp:
        while True:
            chunk = f.read(CHUNK)
            if not chunk:
                break
            fp.write(chunk)
    f.close()
    total_time = time.time() - start
    mb_sec = (os.path.getsize(filename)/(1024*1024.0))/total_time
    print("%s download time: %.2f secs (%.2f MB/sec)" % (filename, total_time, mb_sec))
