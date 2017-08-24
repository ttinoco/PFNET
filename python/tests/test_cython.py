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

class TestCython(unittest.TestCase):

    def setUp(self):

        pass

    def test_network_pickle(self):

        for case in test_cases.CASES:

            # load network from .mat or .raw  files
            if case.endswith('.mat'):
                net1 = pf.ParserMAT().parse(case)
            elif case.endswith('.raw'):
                net1 = pf.ParserRAW().parse(case)

            # Testing pickle string
            pkld_net_string = pickle.dumps(net1, protocol=-1)       #protocol=-1 means HIGHEST_PROTOCOL; 2 or higher needed
            net2 = pickle.loads(pkld_net_string)
            try:
                utils.compare_two_networks(net1, net2)
            except AssertionError:
                self.fail("Networks not identical")

            # Testing pickle with file
            with tempfile.TemporaryFile() as pkld_net_file:
                pickle.dump(net1, pkld_net_file, protocol=-1)
                pkld_net_file.seek(0)
                net3 = pickle.load(pkld_net_file)
                pkld_net_file.close()
                try:
                    utils.compare_two_networks(net1, net3)
                except AssertionError:
                    self.fail("Networks not identical")


    # tests the method from utils.py
    def test_compare_two_networks(self):
        # test with different networks from mat format
        case1 = test_cases.CASES[0]
        case2 = test_cases.CASES[1]
        net1 = pf.ParserMAT().parse(case1)
        net2 = pf.ParserMAT().parse(case2)
        self.assertRaises(AssertionError, utils.compare_two_networks, net1, net2)

    def tearDown(self):

        pass
