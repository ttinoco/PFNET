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
            
            net =pf.Parser(case).parse(case)

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
            
            net =pf.Parser(case).parse(case)

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

        # total gen P

        # total gen Q

        # total gen Qmin
        
        # total gen Qmax

        # total reg gen Qmin

        # total reg gen Qmax

        # reg by gen

        # reg by tran
        
        pass
    
    def test_network(self):

        # clear outages

        # num branches not on outage

        # num branches on outage

        # num gens not on outage

        # num gens on outage

        pass

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

    
