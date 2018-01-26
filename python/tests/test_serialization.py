#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import unittest
import numpy as np
import pickle
import tempfile

import pfnet as pf

from . import test_cases

class TestSerialization(unittest.TestCase):

    def setUp(self):

        pass

    def test_network_pickle(self):

        for case in test_cases.CASES:

            # Load network
            net1 = pf.Parser(case).parse(case)

            # Some modifications
            for gen in net1.generators:
                gen.Q_par = np.random.rand()

            # Testing pickle string
            pkld_net_string = pickle.dumps(net1, protocol=-1)
            net2 = pickle.loads(pkld_net_string)
            pf.tests.utils.compare_networks(self, net1, net2)

            # Testing pickle with file
            with tempfile.TemporaryFile() as pkld_net_file:
                pickle.dump(net1, pkld_net_file, protocol=-1)
                pkld_net_file.seek(0)
                net3 = pickle.load(pkld_net_file)
                pkld_net_file.close()
                pf.tests.utils.compare_networks(self, net1, net3)

    def test_compare_two_networks(self):

        for case in test_cases.CASES:
            net1 = pf.Parser(case).parse(case)
            net2 = pf.Parser(case).parse(case)
            net2.buses[0].v_mag = 100
            self.assertRaises(AssertionError, pf.tests.utils.compare_networks, self, net1, net2)

    def tearDown(self):

        pass
