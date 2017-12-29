#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
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
                self.assertEqual(json_model['index'], bus.index)
                self.assertEqual(json_model['v_base'], bus.v_base)
                self.assertEqual(json_model['name'], bus.name)
                self.assertTrue('sens_P_balance' in json_model)
                self.assertTrue('sens_Q_balance' in json_model)
                self.assertTrue('sens_v_mag_u_bound' in json_model)
                self.assertTrue('sens_v_mag_l_bound' in json_model)
                self.assertTrue('sens_v_ang_u_bound' in json_model)
                self.assertTrue('sens_v_ang_l_bound' in json_model)
                self.assertTrue('sens_v_reg_by_gen' in json_model)
                self.assertTrue('sens_v_reg_by_tran' in json_model)
                self.assertTrue('sens_v_reg_by_shunt' in json_model)
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
                self.assertEqual(json_model['name'], branch.name)
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
                self.assertEqual(json_model['name'], load.name)
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
                self.assertEqual(json_model['name'], gen.name)
                # Add more

    def test_shunt_json_string(self):
                        
        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for shunt in net.shunts:
                text = shunt.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],shunt.index)
                self.assertEqual(json_model['name'], shunt.name)
                # Add more

    def test_vargen_json_string(self):
                        
        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            net.add_var_generators_from_parameters(net.get_load_buses(),100.,50.,30.,5,0.05)
            self.assertGreaterEqual(net.num_var_generators,1)

            for gen in net.var_generators:
                text = gen.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],gen.index)
                self.assertEqual(json_model['name'], gen.name)
                # Add more

    def test_bat_json_string(self):
                        
        # Multiperiod
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            net.add_batteries_from_parameters(net.get_generator_buses(),20.,50.)
            self.assertGreaterEqual(net.num_batteries,1)

            for bat in net.batteries:
                text = bat.json_string
                try:
                    json_model = json.loads(text)
                    valid_json = True
                except ValueError:
                    valid_json = False
                self.assertTrue(valid_json)

                # Detailed checks
                self.assertEqual(json_model['index'],bat.index)
                self.assertEqual(json_model['name'], bat.name)
                # Add more

    def test_net_json_string(self):

        import time
        
        # Multiperiod
        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)
            
            net.add_batteries_from_parameters(net.get_generator_buses(),20.,50.)
            net.add_var_generators_from_parameters(net.get_load_buses(),100.,50.,30.,5,0.05)
            self.assertGreaterEqual(net.num_batteries,1)
            self.assertGreaterEqual(net.num_var_generators,1)

            t0 = time.time()
            text = net.json_string
            t1 = time.time()
            try:
                json_model = json.loads(text)
                valid_json = True
            except ValueError:
                valid_json = False
            self.assertTrue(valid_json)

            # Detailed checks
            self.assertEqual(json_model['num_periods'],self.T)
            self.assertEqual(json_model['base_power'],net.base_power)
            self.assertEqual(json_model['version'],pf.info['version'])
            self.assertEqual(len(json_model['buses']),net.num_buses)
            self.assertEqual(len(json_model['branches']),net.num_branches)
            self.assertEqual(len(json_model['generators']),net.num_generators)
            self.assertEqual(len(json_model['shunts']),net.num_shunts)
            self.assertEqual(len(json_model['var_generators']),net.num_var_generators)
            self.assertEqual(len(json_model['batteries']),net.num_batteries)
            
    def tearDown(self):

        pass
