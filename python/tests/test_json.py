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
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],bus.index)
                # Add more

    def test_branch_json_string(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for branch in net.branches:
                text = branch.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],branch.index)
                # Add more

    def test_load_json_string(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for load in net.loads:
                text = load.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],load.index)
                # Add more

    def test_gen_json_string(self):

        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for gen in net.generators:
                text = gen.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],gen.index)
                # Add more
                
    def tearDown(self):

        pass
