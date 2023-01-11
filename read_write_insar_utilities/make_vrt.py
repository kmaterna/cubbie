#!/usr/bin/env python
# Take an isce-formated file and produce a VRT
# Do this instead of gdalbuildvrt, which is tempting but doesn't produce working ISCE dems for some reason

import sys
import isce
from applications.gdal2isce_xml import gdal2isce_xml 


if __name__ == "__main__":
    filename = sys.argv[1];
    xml_file = gdal2isce_xml(filename);
