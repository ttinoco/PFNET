#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import json
import unittest
import pfnet as pf
from . import test_cases

class TestJSON(unittest.TestCase):

    def setUp(self):

        # Networks
        self.T = 5

    def test_bus_json_string(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for bus in net.buses:
                text = bus.json_string
                try:
                    json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

    def test_branch_json_string(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for branch in net.branches:
                text = branch.json_string
                try:
                    json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)
                
    def tearDown(self):

        pass
