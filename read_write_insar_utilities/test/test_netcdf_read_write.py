import os
import unittest

# Do the read/write functions give the same data afterwards?
# Do the example configs in the repo actually parse through the configparser?
# etc.

from ..netcdf_read_write import read_grd_lonlatz


class NetCDFTests(unittest.TestCase):

    def test_simple_netcdf3(self):
        # HACK because I can't understand why conda's 3.7 dosen't have importlib.resources
        filename = os.path.join(os.path.dirname(__file__), 'USGS_vel.nc3')
        lon, lat, z = read_grd_lonlatz(filename)
        self.assertTrue(lon.any(), 'Could not read longitude of test file')
        self.assertTrue(lat.any(), 'Could not read latitude of test file')
        self.assertTrue(z.any(), 'Could not read z of test file')


if __name__ == '__main__':
    unittest.main()
