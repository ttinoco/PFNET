#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
import test_cases
import numpy as np

class TestFlags(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()
        
    def test_variables(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            # add vargens
            net.add_vargens(net.get_gen_buses(),50.,30.,5,0.05)

            num_vars = 0

            self.assertEqual(net.num_vars,0)

            # Buses
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_SLACK,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            num_vars += 2*net.get_num_slack_buses()
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_GEN,
                          [pf.BUS_VAR_VMAG,pf.BUS_VAR_VANG])
            num_vars += 2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            self.assertEqual(net.num_vars,num_vars)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VANG)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (bus.is_regulated_by_tran() and 
                    not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VANG)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (not bus.is_regulated_by_tran() and 
                    not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_vars,2*net.num_buses)

            # Generators
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_Q)
            num_vars += net.get_num_slack_gens()
            self.assertEqual(net.num_vars,num_vars)

            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            num_vars += 2*net.get_num_reg_gens() - net.get_num_slack_gens()
            self.assertEqual(net.num_vars,num_vars)

            # Branches
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
            num_vars += net.get_num_tap_changers_v()
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            num_vars += net.get_num_phase_shifters()
            self.assertEqual(net.num_vars,num_vars)
            
            # Shunts
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            num_vars += net.get_num_switched_shunts()
            self.assertEqual(net.num_vars,num_vars)

            # Vargens
            self.assertGreater(net.num_vargens,0)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            num_vars += net.num_vargens
            self.assertEqual(net.num_vars,num_vars)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_Q)
            num_vars += net.num_vargens
            self.assertEqual(net.num_vars,num_vars)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            num_vars += 0
            self.assertEqual(net.num_vars,num_vars)

            # Clear
            num_vars = 0
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            # Bus
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_REG_BY_SHUNT,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VVIO)
            num_vars += 3*net.get_num_buses_reg_by_shunt()
            self.assertEqual(net.num_vars,num_vars)
            
    def test_init_point(self):
        
        net = self.net
        
        for case in test_cases.CASES:
            
            net.load(case)

            # add vargens
            net.add_vargens(net.get_gen_buses(),50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = np.random.rand()
                vargen.Q = np.random.rand()

            self.assertEqual(net.num_vars,0)

            point = net.get_var_values()

            self.assertTrue(type(point) is np.ndarray)
            self.assertTupleEqual(point.shape,(0,))
            self.assertEqual(point.size,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_REG_BY_GEN,
                          pf.BUS_VAR_VMAG)
            num_vars = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if not bus.is_regulated_by_gen():
                    num_vars += 1
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_Q)
            num_vars += net.num_gens
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            num_vars += net.get_num_switched_shunts()
            self.assertEqual(num_vars,net.num_vars)
           
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            num_vars += net.num_vargens*2
            self.assertEqual(num_vars,net.num_vars)
 
            point = net.get_var_values()
            self.assertTrue(type(point) is np.ndarray)
            self.assertTupleEqual(point.shape,(num_vars,))
            self.assertEqual(point.size,num_vars)
            self.assertTrue(not np.any(np.isinf(point)))
            self.assertTrue(not np.any(np.isnan(point)))

            # check switched shunt vals
            for i in range(net.num_shunts):
                s = net.get_shunt(i)
                if s.is_switched_v():
                    self.assertLess(np.abs(point[s.index_b]-s.b),1e-10)

            # check vargens
            self.assertGreater(len(net.var_generators),0)
            for vargen in net.var_generators:
                self.assertNotEqual(vargen.P,0.)
                self.assertNotEqual(vargen.Q,0.)
                self.assertLess(np.abs(point[vargen.index_P]-vargen.P),1e-10)
                self.assertLess(np.abs(point[vargen.index_Q]-vargen.Q),1e-10)

    def test_tap_changer_v(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            num_vars = 0

            self.assertEqual(net.num_vars,0)

            # Tap changer v
            manual_num_vars = 0
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_tran():
                    self.assertGreater(len(bus.reg_trans),0)
                    self.assertTrue(all([t.is_tap_changer_v() for t in bus.reg_trans]))
                    for tran in bus.reg_trans:
                        manual_num_vars += 1
                        self.assertTrue(tran.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
            self.assertEqual(manual_num_vars,net.num_vars)
            
    def test_bounded(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_bounded,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_REG_BY_TRAN,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_BOUNDED,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
                    
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_BOUNDED,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)
            
            self.assertEqual(net.num_bounded,
                             (net.get_num_tap_changers_v() +
                              net.get_num_buses_reg_by_tran() + 
                              net.get_num_slack_gens() +
                              net.num_shunts))

    def test_fixed(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            net.add_vargens(net.get_gen_buses(),50.,30.,5,0.05)

            self.assertEqual(net.num_fixed,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_REG_BY_GEN,
                          [pf.BUS_VAR_VANG,pf.BUS_VAR_VMAG])

            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_FIXED,
                          pf.GEN_PROP_NOT_REG,
                          pf.GEN_VAR_Q)

            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)

            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_FIXED,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)

            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_FIXED,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            
            self.assertEqual(net.num_fixed,
                             (net.get_num_buses_reg_by_gen()*2 + 
                              net.num_gens - net.get_num_reg_gens() + 
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()+
                              net.num_vargens*2))

    def test_multiple_flags(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)            

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,0)
            
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS|pf.FLAG_SPARSE,
                          pf.BUS_PROP_REG_BY_GEN,
                          pf.BUS_VAR_VMAG)

            self.assertEqual(net.num_vars,net.get_num_buses_reg_by_gen())
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,net.get_num_buses_reg_by_gen())
            self.assertEqual(net.num_bounded,0)

            net.clear_flags()

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_FIXED|pf.FLAG_BOUNDED,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,2*net.get_num_slack_gens())
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,2*net.get_num_slack_gens())

    def test_custom_flags(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)            

            # add vargens
            net.add_vargens(net.get_gen_buses(),50.,30.,5,0.05)
            self.assertGreater(net.num_vargens,0)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,0)

            num_vars = 0
            num_bounded = 0
            num_fixed = 0

            # buses
            for bus in net.buses:
                if bus.index % 3 == 0:
                    net.set_flags_of_component(bus,
                                               pf.FLAG_VARS|pf.FLAG_BOUNDED,
                                               pf.BUS_VAR_VANG|pf.BUS_VAR_VMAG)
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG))
                    self.assertTrue(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))
                    self.assertFalse(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VVIO))
                    self.assertTrue(bus.has_flags(pf.FLAG_BOUNDED,pf.BUS_VAR_VANG))
                    self.assertTrue(bus.has_flags(pf.FLAG_BOUNDED,pf.BUS_VAR_VMAG))
                    self.assertFalse(bus.has_flags(pf.FLAG_BOUNDED,pf.BUS_VAR_VDEV))
                    self.assertFalse(bus.has_flags(pf.FLAG_FIXED,pf.BUS_VAR_VANG))
                    self.assertFalse(bus.has_flags(pf.FLAG_FIXED,pf.BUS_VAR_VMAG))
                    num_vars += 2
                    num_bounded += 2
                else:
                    self.assertFalse(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG))
                    self.assertFalse(bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_bounded,num_bounded)

            # branches
            for branch in net.branches:
                if branch.index == 5:
                    net.set_flags_of_component(branch,
                                               pf.FLAG_FIXED|pf.FLAG_VARS,
                                               pf.BRANCH_VAR_RATIO)
                    self.assertTrue(branch.has_flags(pf.FLAG_FIXED,pf.BRANCH_VAR_RATIO))
                    self.assertTrue(branch.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
                    self.assertFalse(branch.has_flags(pf.FLAG_FIXED,pf.BRANCH_VAR_PHASE))
                    self.assertFalse(branch.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_PHASE))
                    num_fixed += 1
                    num_vars += 1
                else:
                    self.assertFalse(branch.has_flags(pf.FLAG_FIXED,pf.BRANCH_VAR_RATIO))
                    self.assertFalse(branch.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

            # gens
            for gen in net.generators:
                if not gen.is_regulator():
                    net.set_flags_of_component(gen,
                                               pf.FLAG_VARS,
                                               pf.GEN_VAR_P)
                    self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                    self.assertFalse(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_Q))
                    self.assertFalse(gen.has_flags(pf.FLAG_FIXED,pf.GEN_VAR_P))
                else:
                    self.assertFalse(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
            num_vars += net.num_gens-net.get_num_reg_gens()
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

            # loads

            # shunts
            for shunt in net.shunts:
                if shunt.index % 2 == 0:
                    net.set_flags_of_component(shunt,
                                               pf.FLAG_BOUNDED,
                                               pf.SHUNT_VAR_SUSC)
                    self.assertTrue(shunt.has_flags(pf.FLAG_BOUNDED,pf.SHUNT_VAR_SUSC))
                    self.assertFalse(shunt.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC))
                    self.assertFalse(shunt.has_flags(pf.FLAG_FIXED,pf.SHUNT_VAR_SUSC))
                    num_bounded += 1
                else:
                    self.assertFalse(shunt.has_flags(pf.FLAG_BOUNDED,pf.SHUNT_VAR_SUSC))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)
            self.assertEqual(net.num_bounded,num_bounded)

            # vargens
            for vargen in net.var_generators:
                if vargen.index % 3 == 0:
                    net.set_flags_of_component(vargen,
                                               pf.FLAG_VARS,
                                               pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
                    net.set_flags_of_component(vargen,
                                               pf.FLAG_VARS,
                                               pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
                    self.assertTrue(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P))
                    self.assertTrue(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_Q))
                    self.assertFalse(vargen.has_flags(pf.FLAG_FIXED,pf.VARGEN_VAR_P))
                    self.assertFalse(vargen.has_flags(pf.FLAG_FIXED,pf.VARGEN_VAR_Q))
                    num_vars += 2
                else:
                    self.assertFalse(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P))
                    self.assertFalse(vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_Q))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)
            self.assertEqual(net.num_bounded,num_bounded)

    def test_errors(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertRaises(pf.NetworkError,
                              net.set_flags,
                              -1,
                              pf.FLAG_VARS,
                              pf.BUS_PROP_ANY,
                              pf.BUS_VAR_VANG)
            self.assertRaises(pf.NetworkError,
                              net.set_flags,
                              pf.OBJ_BUS,
                              -1,
                              pf.BUS_PROP_ANY,
                              pf.BUS_VAR_VANG)
            self.assertRaises(pf.NetworkError,
                              net.set_flags,
                              pf.OBJ_GEN,
                              -1,
                              pf.GEN_PROP_ANY,
                              pf.GEN_VAR_P)
            self.assertRaises(pf.NetworkError,
                              net.set_flags,
                              pf.OBJ_BRANCH,
                              -1,
                              pf.BRANCH_PROP_TAP_CHANGER_V,
                              pf.BRANCH_VAR_RATIO)
            
    def tearDown(self):
        
        pass




