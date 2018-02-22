#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import json
import unittest
import numpy as np
import pfnet as pf
from . import test_cases

class TestOutages(unittest.TestCase):

    def test_generators(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)

            net.clear_outages()
            
            for gen in net.generators:

                self.assertFalse(gen.is_on_outage())
                self.assertFalse(gen.outage)
                
                reg = gen.is_regulator()
                slack = gen.is_slack()

                # set
                gen.outage = True

                self.assertTrue(gen.is_on_outage())
                self.assertTrue(gen.outage)

                # bus
                self.assertTrue(gen.bus is not None)
                self.assertTrue(gen.index in [g.index for g in gen.bus.generators])

                # regulation
                if reg:
                    self.assertTrue(gen.is_regulator())
                    self.assertTrue(gen.reg_bus is not None)
                    self.assertTrue(gen.index in [g.index for g in gen.reg_bus.reg_generators])

                    if all([g.is_on_outage() for g in gen.reg_bus.reg_generators]):
                        self.assertFalse(gen.reg_bus.is_regulated_by_gen())
                    else:
                        self.assertTrue(gen.reg_bus.is_regulated_by_gen())

                # slack
                if slack:
                    self.assertTrue(gen.is_slack())
                    self.assertTrue(gen.bus.is_slack())

                # clear
                gen.outage = False
                self.assertFalse(gen.is_on_outage())

                # set
                gen.outage = True
            
                # json
                json_string = gen.json_string
                d = json.loads(json_string)
                self.assertTrue('outage' in d)
                self.assertTrue(d['outage'])

            # copy
            new_net = net.get_copy()
            for new_gen in new_net.generators:
                self.assertTrue(new_gen.is_on_outage())
                
    def test_branches(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)

            net.clear_outages()

            for branch in net.branches:

                self.assertFalse(branch.is_on_outage())
                self.assertFalse(branch.outage)

                reg = branch.is_tap_changer_v()
                
                # set
                branch.outage = True
                self.assertTrue(branch.is_on_outage())
                self.assertTrue(branch.outage)

                # buses
                self.assertTrue(branch.bus_k is not None)
                self.assertTrue(branch.bus_m is not None)
                self.assertTrue(branch.index in [br.index for br in branch.bus_k.branches_k])
                self.assertTrue(branch.index in [br.index for br in branch.bus_m.branches_m])

                # regulation
                if reg:
                    self.assertTrue(branch.is_tap_changer_v())
                    self.assertTrue(branch.reg_bus is not None)
                    self.assertTrue(branch.index in [br.index for br in branch.reg_bus.reg_trans])

                    if all([br.is_on_outage() for br in branch.reg_bus.reg_trans]):
                        self.assertFalse(branch.reg_bus.is_regulated_by_tran())
                    else:
                        self.assertTrue(branch.reg_bus.is_regulated_by_tran())

                # clear
                branch.outage = False
                self.assertFalse(branch.is_on_outage())
                self.assertFalse(branch.outage)
                
                # set
                branch.outage = True
                
                # json
                json_string = branch.json_string
                d = json.loads(json_string)
                self.assertTrue('outage' in d)
                self.assertTrue(d['outage'])
                
            # copy  
            new_net = net.get_copy()
            for new_branch in new_net.branches:
                self.assertTrue(new_branch.is_on_outage())
                
    def test_buses(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)

            for bus in net.buses:

                net.clear_outages()

                # total gen P
                total = bus.get_total_gen_P()
                for gen in bus.generators:
                    gen.outage = True
                    new_total = bus.get_total_gen_P()
                    self.assertLess(np.abs(new_total-(total-gen.P)),1e-8)
                    total = new_total
                    
                net.clear_outages()
                
                # total gen Q
                total = bus.get_total_gen_Q()
                for gen in bus.generators:
                    gen.outage = True
                    new_total = bus.get_total_gen_Q()
                    self.assertLess(np.abs(new_total-(total-gen.Q)),1e-8)
                    total = new_total

                net.clear_outages()
                
                # total gen Qmin
                total = bus.get_total_gen_Q_min()
                for gen in bus.generators:
                    gen.outage = True
                    new_total = bus.get_total_gen_Q_min()
                    self.assertLess(np.abs(new_total-(total-gen.Q_min)),1e-8)
                    total = new_total

                net.clear_outages()
                
                # total gen Qmax
                total = bus.get_total_gen_Q_max()
                for gen in bus.generators:
                    gen.outage = True
                    new_total = bus.get_total_gen_Q_max()
                    self.assertLess(np.abs(new_total-(total-gen.Q_max)),1e-8)
                    total = new_total

                net.clear_outages()

                # toatl reg gen Q
                total = bus.get_total_reg_gen_Q()
                self.assertLess(np.abs(total-sum([g.Q for g in bus.reg_generators])), 1e-8)
                for gen in bus.reg_generators:
                    gen.outage = True
                    new_total = bus.get_total_reg_gen_Q()
                    self.assertLess(np.abs(new_total-(total-gen.Q)),1e-8)
                    total = new_total

                net.clear_outages()
                
                # total reg gen Qmin
                total = bus.get_total_reg_gen_Q_min()
                self.assertLess(np.abs(total-sum([g.Q_min for g in bus.reg_generators])), 1e-8)
                for gen in bus.reg_generators:
                    gen.outage = True
                    new_total = bus.get_total_reg_gen_Q_min()
                    self.assertLess(np.abs(new_total-(total-gen.Q_min)),1e-8)
                    total = new_total

                net.clear_outages()
                
                # total reg gen Qmax
                total = bus.get_total_reg_gen_Q_max()
                self.assertLess(np.abs(total-sum([g.Q_max for g in bus.reg_generators])), 1e-8)
                for gen in bus.reg_generators:
                    gen.outage = True
                    new_total = bus.get_total_reg_gen_Q_max()
                    self.assertLess(np.abs(new_total-(total-gen.Q_max)),1e-8)
                    total = new_total

                net.clear_outages()
                
                # reg by gen
                if bus.is_regulated_by_gen():
                    for i in range(len(bus.reg_generators)):
                        gen = bus.reg_generators[i]
                        gen.outage = True
                        if i < len(bus.reg_generators)-1:
                            self.assertTrue(bus.is_regulated_by_gen())
                        else:
                            self.assertFalse(bus.is_regulated_by_gen())

                net.clear_outages()
                            
                # reg by tran
                if bus.is_regulated_by_tran():
                    for i in range(len(bus.reg_trans)):
                        br = bus.reg_trans[i]
                        br.outage = True
                        if i < len(bus.reg_trans)-1:
                            self.assertTrue(bus.is_regulated_by_tran())
                        else:
                            self.assertFalse(bus.is_regulated_by_tran())

                net.clear_outages()
                            
                # slack
                if bus.is_slack():
                    for gen in bus.generators:
                        gen.outage = True
                        self.assertTrue(gen.is_on_outage())
                        self.assertTrue(gen.is_slack())
                        self.assertTrue(bus.is_slack())
    
    def test_network(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)

            # clear outages
            net.clear_outages()

            self.assertEqual(net.get_num_branches_not_on_outage(), net.num_branches)
            self.assertEqual(net.get_num_branches_on_outage(), 0)
            self.assertEqual(net.get_num_generators_not_on_outage(), net.num_generators)
            self.assertEqual(net.get_num_generators_on_outage(), 0)
            
            # num branches on outage
            for branch in net.branches:
                branch.outage = True
            self.assertEqual(net.get_num_branches_not_on_outage(), 0)
            self.assertEqual(net.get_num_branches_on_outage(), net.num_branches)
            
            # num gens on outage
            for gen in net.generators:
                gen.outage = True
            self.assertEqual(net.get_num_generators_not_on_outage(), 0)
            self.assertEqual(net.get_num_generators_on_outage(), net.num_generators)

    def test_network_properties(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            
            # tap ratio vio
            net.clear_outages()
            for branch in net.branches:
                if branch.is_tap_changer():
                    branch.ratio = branch.ratio_max + 10.
                    net.update_properties()
                    self.assertEqual(net.tran_r_vio, 10.)
                    branch.outage = True
                    net.update_properties()
                    self.assertNotEqual(net.tran_r_vio, 10.)
                    break
            
            # phase shift vio
            net.clear_outages() 
            for branch in net.branches:
                if branch.is_phase_shifter():
                    branch.phase = branch.phase_max + 20.
                    net.update_properties()
                    self.assertEqual(net.tran_p_vio, 20.)
                    branch.outage = True
                    net.update_properties()
                    self.assertNotEqual(net.tran_p_vio, 20.)
                    break
            
            # mismatches (branch and gen outages)
            net.clear_outages()
            net.update_properties()
            bus = net.get_bus(0)
            p_mis = bus.P_mismatch
            q_mis = bus.Q_mismatch
            for branch in bus.branches:
                branch.outage = True
            for gen in bus.generators:
                gen.outage = True
            net.update_properties()
            self.assertLess(np.abs(p_mis-(bus.P_mismatch + 
                                          sum([g.P for g in bus.generators]) -
                                          sum([br.P_km for br in bus.branches_k]) -
                                          sum([br.P_mk for br in bus.branches_m]))), 1e-8)
            self.assertLess(np.abs(q_mis-(bus.Q_mismatch +
                                          sum([g.Q for g in bus.generators]) -
                                          sum([br.Q_km for br in bus.branches_k]) -
                                          sum([br.Q_mk for br in bus.branches_m]))), 1e-8)
            
            # v reg limit violations
            net.clear_outages()
            for bus in net.buses:
                if bus.is_regulated_by_tran():
                    net.update_properties()
                    bus.v_mag = bus.v_max_reg + 15.
                    net.update_properties()
                    self.assertLess(np.abs(net.tran_v_vio-15.), 1e-8)
                    for branch in bus.reg_trans:
                        branch.outage = True
                    self.assertFalse(bus.is_regulated_by_tran())
                    net.update_properties()
                    self.assertGreater(np.abs(net.tran_v_vio-15.), 1e-8)
                    break
                            
            # v set deviations
            net.clear_outages()
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    bus.v_mag = bus.v_set + 33.
                    net.update_properties()
                    self.assertLess(np.abs(net.gen_v_dev-33.), 1e-8)
                    for gen in bus.reg_generators:
                        gen.outage = True
                    self.assertFalse(bus.is_regulated_by_gen())
                    net.update_properties()
                    self.assertGreater(np.abs(net.gen_v_dev-33.), 1e-8)
                    break
            
            # gen active power cost
            net.clear_outages()
            net.update_properties()
            cost = net.gen_P_cost
            for gen in net.generators:
                gen.outage = True
                net.update_properties()
                cost -= gen.cost_coeff_Q0 + gen.cost_coeff_Q1*gen.P + gen.cost_coeff_Q2*(gen.P**2.)
                self.assertLess(np.abs(cost-net.gen_P_cost), 1e-8)
            
            # gen Q vio
            net.clear_outages()
            for gen in net.generators:
                if gen.is_regulator():
                    gen.Q = gen.Q_max + 340.
                    net.update_properties()
                    self.assertLess(np.abs(net.gen_Q_vio-340.*net.base_power), 1e-8)
                    gen.outage = True
                    net.update_properties()
                    self.assertGreater(np.abs(net.gen_Q_vio-340.*net.base_power), 1e-8)
                    break
            
            # gen P vio
            net.clear_outages()
            for gen in net.generators:
                gen.P = gen.P_min - 540.
                net.update_properties()
                self.assertLess(np.abs(net.gen_P_vio-540.*net.base_power), 1e-8)
                gen.outage = True
                net.update_properties()
                self.assertGreater(np.abs(net.gen_P_vio-540.*net.base_power), 1e-8)
                break

    def test_functions(self):

        # gen cost
        
        # load util

        # netcon cost

        # reg phase

        # reg pq

        # reg ratio

        # reg susc

        # reg vang

        # reg var

        # reg vmag

        # reg slim vmag

        # sp controls

        pass

    def test_constraints(self):

        # ac flow limit

        # ac branch flow limit

        # ac pf

        # bat dyn

        # dc flow lim

        # dc pf

        # fix

        # gen ramp

        # lbound

        # linpf

        # load pf

        # nbound

        # par gen P

        # PVPQ switching

        # reg gen

        # reg shunt

        # reg tran

        pass

    
