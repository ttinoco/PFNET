#***************************************************#
#
#
#
#
#***************************************************#

import unittest
import numpy as np
import pickle
import tempfile

# local
import pfnet as pf
import utils

from . import test_cases

class TestSerialization(unittest.TestCase):

    def setUp(self):

        pass

    def test_network_pickle(self):

        for case in test_cases.CASES:

            # load network
            net1 = pf.Parser(case).parse(case)

            # Testing pickle string
            pkld_net_string = pickle.dumps(net1, protocol=-1)       #protocol=-1 means HIGHEST_PROTOCOL; 2 or higher needed
            net2 = pickle.loads(pkld_net_string)
            utils.compare_two_networks(self, net1, net2)

            # Testing pickle with file
            with tempfile.TemporaryFile() as pkld_net_file:
                pickle.dump(net1, pkld_net_file, protocol=-1)
                pkld_net_file.seek(0)
                net3 = pickle.load(pkld_net_file)
                pkld_net_file.close()
                utils.compare_two_networks(self, net1, net3)

    # tests the method from utils.py
    def test_compare_two_networks(self):
        case = test_cases.CASES[0]
        net1 = pf.Parser(case).parse(case)
        net2 = pf.Parser(case).parse(case)
        net2.buses[0].v_mag = 100
        self.assertRaises(AssertionError, utils.compare_two_networks, self, net1, net2)

    def tearDown(self):

        pass
