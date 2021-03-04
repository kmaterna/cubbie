# Do the example configs in the repo actually parse through the configparser?

import unittest
from intf_generating import sentinel_main_functions


class ConfigTests(unittest.TestCase):

    def test_reading_config(self):
        config_file = "configs_and_setup/batch.config"
        print("Parsing example config file %s " % config_file)
        opt1, opt2 = sentinel_main_functions.read_config(config_file);
        self.assertIsNotNone(opt1, "ERROR reading example stacking config file");
        self.assertIsNotNone(opt2, "ERROR reading example stacking config file");


if __name__ == "__main__":
    unittest.main();
