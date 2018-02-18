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
                self.assertTrue(d.has_key('outage'))
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
                        self.assertFalse(br.reg_bus.is_regulated_by_tran())
                    else:
                        self.assertTrue(br.reg_bus.is_regulated_by_tran())

                # clear
                branch.outage = False
                self.assertFalse(branch.is_on_outage())
                self.assertFalse(branch.outage)
                
                # set
                branch.outage = True
                
                # json
                json_string = branch.json_string
                d = json.loads(json_string)
                self.assertTrue(d.has_key('outage'))
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

        # tap ratio vio and actions

        # phase shift vio and actions

        # mismatches (branch and gen outages)

        # v reg limit violations

        # v set deviations and actions

        # gen active power cost

        # gen Q vio

        # gen P vio and actions

        pass

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

    
