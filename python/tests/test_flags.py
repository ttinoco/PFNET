#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
from . import test_cases
import numpy as np

class TestFlags(unittest.TestCase):
    
    def setUp(self):
        
        pass

    def test_variables(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # add vargens
            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)

            # add batteries
            net.add_batteries(net.get_generator_buses(),20.,50.)            

            # loads
            lcount = 0
            for load in net.loads:
                if load.index % 2 == 0:
                    load.P_min = -10.
                    load.P_max = 20.
                    lcount += 1
            self.assertEqual(lcount,net.get_num_P_adjust_loads())
            self.assertGreater(lcount,0)

            num_vars = 0

            self.assertEqual(net.num_vars,0)

            # Buses
            net.set_flags('bus',
                          'variable',
                          'slack',
                          ['voltage magnitude','voltage angle'])
            num_vars += 2*net.get_num_slack_buses()
            self.assertEqual(net.num_vars,num_vars)
            for bus in net.buses:
                self.assertEqual(bus.flags_fixed,0)
                self.assertEqual(bus.flags_sparse,0)
                self.assertEqual(bus.flags_bounded,0)
                if bus.is_slack():
                    self.assertTrue(bus.has_flags('variable',['voltage magnitude','voltage angle']))
                    self.assertNotEqual(bus.flags_vars,0)
                else:
                    self.assertEqual(bus.flags_vars,0)
            
            net.set_flags('bus',
                          'variable',
                          'regulated by generator',
                          ['voltage magnitude','voltage angle'])
            num_vars += 2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses())
            self.assertEqual(net.num_vars,num_vars)

            net.set_flags('bus',
                          'variable',
                          'regulated by transformer',
                          'voltage angle')
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (bus.is_regulated_by_tran() and 
                    not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            
            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage angle')
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if (not bus.is_regulated_by_tran() and 
                    not bus.is_regulated_by_gen() and
                    not bus.is_slack()):
                    num_vars +=1
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_vars,2*net.num_buses)

            # Generators
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'reactive power')
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'reactive power')
            num_vars += net.get_num_slack_gens()
            self.assertEqual(net.num_vars,num_vars)
            for gen in net.generators:
                self.assertEqual(gen.flags_fixed,0)
                self.assertEqual(gen.flags_sparse,0)
                self.assertEqual(gen.flags_bounded,0)
                if gen.is_slack():
                    self.assertTrue(gen.has_flags('variable','reactive power'))
                    self.assertNotEqual(gen.flags_vars,0)
                else:
                    self.assertEqual(gen.flags_vars,0)
                    
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          ['active power','reactive power'])
            num_vars += 2*net.get_num_reg_gens() - net.get_num_slack_gens()
            self.assertEqual(net.num_vars,num_vars)

            # Loads
            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          'active power')
            num_vars += net.get_num_P_adjust_loads()
            self.assertEqual(net.num_vars,num_vars)
            for load in net.loads:
                self.assertEqual(load.flags_fixed,0)
                self.assertEqual(load.flags_sparse,0)
                self.assertEqual(load.flags_bounded,0)
                if load.is_P_adjustable():
                    self.assertTrue(load.has_flags('variable','active power'))
                    self.assertNotEqual(load.flags_vars,0)
                else:
                    self.assertEqual(load.flags_vars,0)

            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          'active power')
            self.assertEqual(net.num_vars,num_vars)

            # Branches
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            num_vars += net.get_num_tap_changers_v()
            self.assertEqual(net.num_vars,num_vars)
            for branch in net.branches:
                self.assertEqual(branch.flags_fixed,0)
                self.assertEqual(branch.flags_sparse,0)
                self.assertEqual(branch.flags_bounded,0)
                if branch.is_tap_changer_v():
                    self.assertTrue(branch.has_flags('variable','tap ratio'))
                    self.assertNotEqual(branch.flags_vars,0)
                else:
                    self.assertEqual(branch.flags_vars,0)
            
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            num_vars += net.get_num_phase_shifters()
            self.assertEqual(net.num_vars,num_vars)
            
            # Shunts
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            num_vars += net.get_num_switched_shunts()
            self.assertEqual(net.num_vars,num_vars)
            for shunt in net.shunts:
                self.assertEqual(shunt.flags_fixed,0)
                self.assertEqual(shunt.flags_sparse,0)
                self.assertEqual(shunt.flags_bounded,0)
                if shunt.is_switched_v():
                    self.assertTrue(shunt.has_flags('variable','susceptance'))
                    self.assertNotEqual(shunt.flags_vars,0)
                else:
                    self.assertEqual(shunt.flags_vars,0)

            # Vargens
            self.assertGreater(net.num_var_generators,0)
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            num_vars += net.num_var_generators
            self.assertEqual(net.num_vars,num_vars)
            for vargen in net.var_generators:
                self.assertEqual(vargen.flags_fixed,0)
                self.assertEqual(vargen.flags_sparse,0)
                self.assertEqual(vargen.flags_bounded,0)
                self.assertNotEqual(vargen.flags_vars,0)
                self.assertTrue(vargen.has_flags('variable','active power'))
            
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'reactive power')
            num_vars += net.num_var_generators
            self.assertEqual(net.num_vars,num_vars)
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            num_vars += 0
            self.assertEqual(net.num_vars,num_vars)

            # Batteries
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            num_vars += 3*net.num_batteries
            self.assertEqual(net.num_vars,num_vars)
            for bat in net.batteries:
                self.assertEqual(bat.flags_fixed,0)
                self.assertEqual(bat.flags_sparse,0)
                self.assertEqual(bat.flags_bounded,0)
                self.assertNotEqual(bat.flags_vars,0)
                self.assertTrue(bat.has_flags('variable',['charging power','energy level']))

            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            num_vars += 0
            self.assertEqual(net.num_vars,num_vars)

            # Clear
            num_vars = 0
            net.clear_flags()
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            # Bus
            net.set_flags('bus',
                          'variable',
                          'regulated by shunt',
                          ['voltage magnitude'])
            num_vars += net.get_num_buses_reg_by_shunt()
            self.assertEqual(net.num_vars,num_vars)
            
    def test_init_point(self):
        
        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # add vargens
            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                vargen.P = np.random.rand()
                vargen.Q = np.random.rand()

            self.assertEqual(net.num_vars,0)

            point = net.get_var_values()

            self.assertTrue(type(point) is np.ndarray)
            self.assertTupleEqual(point.shape,(0,))
            self.assertEqual(point.size,0)

            net.set_flags('bus',
                          'variable',
                          'not regulated by generator',
                          'voltage magnitude')
            num_vars = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if not bus.is_regulated_by_gen():
                    num_vars += 1
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags('generator',
                          'variable',
                          'any',
                          'reactive power')
            num_vars += net.num_generators
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags('load',
                          'variable',
                          'any',
                          'active power')
            num_vars += net.num_loads
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            num_vars += net.get_num_switched_shunts()
            self.assertEqual(num_vars,net.num_vars)
           
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            num_vars += net.num_var_generators*2
            self.assertEqual(num_vars,net.num_vars)

            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            num_vars += 3*net.num_batteries
            self.assertEqual(net.num_vars,num_vars)
 
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

            # check loads
            for i in range(net.num_loads):
                l = net.get_load(i)
                self.assertTrue(l.has_flags('variable','active power'))
                self.assertEqual(point[l.index_P],l.P)

            # check bats
            for i in range(net.num_batteries):
                b = net.get_battery(i)
                self.assertTrue(b.has_flags('variable','charging power'))
                self.assertTrue(b.has_flags('variable','energy level'))
                if (b.P >= 0.):
                    self.assertEqual(point[b.index_Pc],b.P)
                    self.assertEqual(point[b.index_Pd],0.)
                else:
                    self.assertEqual(point[b.index_Pd],-b.P)
                    self.assertEqual(point[b.index_Pc],0.)
                self.assertEqual(point[b.index_E],b.E)

    def test_tap_changer_v(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            num_vars = 0

            self.assertEqual(net.num_vars,0)

            # Tap changer v
            manual_num_vars = 0
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.is_regulated_by_tran():
                    self.assertGreater(len(bus.reg_trans),0)
                    self.assertTrue(all([t.is_tap_changer_v() for t in bus.reg_trans]))
                    for tran in bus.reg_trans:
                        manual_num_vars += 1
                        self.assertTrue(tran.has_flags('variable','tap ratio'))
            self.assertEqual(manual_num_vars,net.num_vars)
            
    def test_bounded(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_bounded,0)

            net.set_flags('bus',
                          'bounded',
                          'regulated by transformer',
                          'voltage magnitude')
            net.set_flags('generator',
                          'bounded',
                          'slack',
                          'active power')
            net.set_flags('load',
                          'bounded',
                          'any',
                          'active power')
            net.set_flags('branch',
                          'bounded',
                          'tap changer - v',
                          'tap ratio')                    
            net.set_flags('shunt',
                          'bounded',
                          'any',
                          'susceptance')
            net.set_flags('battery',
                          'bounded',
                          'any',
                          'energy level')
            
            self.assertEqual(net.num_bounded,
                             (net.get_num_tap_changers_v() +
                              net.get_num_buses_reg_by_tran() + 
                              net.get_num_slack_gens() +
                              net.num_shunts +
                              net.num_loads+
                              net.num_batteries))
                             
            # loads
            for load in net.loads:
                self.assertTrue(load.has_flags('bounded','active power'))

            # bats
            for bat in net.batteries:
                self.assertTrue(bat.has_flags('bounded','energy level'))
                self.assertFalse(bat.has_flags('bounded','charging power'))

    def test_fixed(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)

            self.assertEqual(net.num_fixed,0)

            net.set_flags('bus',
                          'fixed',
                          'regulated by generator',
                          ['voltage angle','voltage magnitude'])

            net.set_flags('generator',
                          'fixed',
                          'not regulator',
                          'reactive power')

            net.set_flags('load',
                          'fixed',
                          'any',
                          'active power')

            net.set_flags('branch',
                          'fixed',
                          'phase shifter',
                          'phase shift')

            net.set_flags('shunt',
                          'fixed',
                          'switching - v',
                          'susceptance')

            net.set_flags('variable generator',
                          'fixed',
                          'any',
                          ['active power','reactive power'])

            net.set_flags('battery',
                          'fixed',
                          'any',
                          ['charging power','energy level'])
            
            self.assertEqual(net.num_fixed,
                             (net.get_num_buses_reg_by_gen()*2 + 
                              net.num_generators - net.get_num_reg_gens() + 
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts()+
                              net.num_var_generators*2 +
                              net.num_loads+
                              3*net.num_batteries))

            # loads
            for load in net.loads:
                self.assertFalse(load.has_flags('bounded','active power'))
                self.assertFalse(load.has_flags('variable','active power'))
                self.assertTrue(load.has_flags('fixed','active power'))

            # batteries
            for bat in net.batteries:
                self.assertFalse(bat.has_flags('bounded','charging power'))
                self.assertFalse(bat.has_flags('variable','charging power'))
                self.assertTrue(bat.has_flags('fixed','charging power'))
                self.assertFalse(bat.has_flags('bounded','energy level'))
                self.assertFalse(bat.has_flags('variable','energy level'))
                self.assertTrue(bat.has_flags('fixed','energy level'))

    def test_multiple_flags(self):

        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,0)
            
            net.set_flags('bus',
                          ['variable','sparse'],
                          'regulated by generator',
                          'voltage magnitude')

            self.assertEqual(net.num_vars,net.get_num_buses_reg_by_gen())
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,net.get_num_buses_reg_by_gen())
            self.assertEqual(net.num_bounded,0)

            net.clear_flags()

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags('generator',
                          ['fixed','bounded'],
                          'slack',
                          ['active power','reactive power'])

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,2*net.get_num_slack_gens())
            self.assertEqual(net.num_sparse,0)
            self.assertEqual(net.num_bounded,2*net.get_num_slack_gens())

    def test_custom_flags(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # add vargens
            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)

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
                                               ['variable','bounded'],
                                               ['voltage angle','voltage magnitude'])
                    self.assertTrue(bus.has_flags('variable','voltage angle'))
                    self.assertTrue(bus.has_flags('variable','voltage magnitude'))
                    self.assertTrue(bus.has_flags('bounded','voltage angle'))
                    self.assertTrue(bus.has_flags('bounded','voltage magnitude'))
                    self.assertFalse(bus.has_flags('fixed','voltage angle'))
                    self.assertFalse(bus.has_flags('fixed','voltage magnitude'))
                    num_vars += 2
                    num_bounded += 2
                else:
                    self.assertFalse(bus.has_flags('variable','voltage angle'))
                    self.assertFalse(bus.has_flags('variable','voltage magnitude'))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_bounded,num_bounded)

            # branches
            for branch in net.branches:
                if branch.index == 5:
                    net.set_flags_of_component(branch,
                                               ['fixed','variable'],
                                               'tap ratio')
                    self.assertTrue(branch.has_flags('fixed','tap ratio'))
                    self.assertTrue(branch.has_flags('variable','tap ratio'))
                    self.assertFalse(branch.has_flags('fixed','phase shift'))
                    self.assertFalse(branch.has_flags('variable','phase shift'))
                    num_fixed += 1
                    num_vars += 1
                else:
                    self.assertFalse(branch.has_flags('fixed','tap ratio'))
                    self.assertFalse(branch.has_flags('variable','tap ratio'))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

            # gens
            for gen in net.generators:
                if not gen.is_regulator():
                    net.set_flags_of_component(gen,
                                               'variable',
                                               'active power')
                    self.assertTrue(gen.has_flags('variable','active power'))
                    self.assertFalse(gen.has_flags('variable','reactive power'))
                    self.assertFalse(gen.has_flags('fixed','active power'))
                else:
                    self.assertFalse(gen.has_flags('variable','active power'))
            num_vars += net.num_generators-net.get_num_reg_gens()
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

            # loads
            lcount = 0
            for load in net.loads:
                if load.index % 5 == 0:
                    net.set_flags_of_component(load,
                                               'variable',
                                               'active power')
                    lcount += 1
            for load in net.loads:
                if load.index % 5 == 0:
                    self.assertTrue(load.has_flags('variable','active power'))
                    self.assertFalse(load.has_flags('fixed','active power'))
                else:
                    self.assertFalse(load.has_flags('variable','active power'))
            num_vars += lcount
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

            # shunts
            for shunt in net.shunts:
                if shunt.index % 2 == 0:
                    net.set_flags_of_component(shunt,
                                               'bounded',
                                               'susceptance')
                    self.assertTrue(shunt.has_flags('bounded','susceptance'))
                    self.assertFalse(shunt.has_flags('variable','susceptance'))
                    self.assertFalse(shunt.has_flags('fixed','susceptance'))
                    num_bounded += 1
                else:
                    self.assertFalse(shunt.has_flags('bounded','susceptance'))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)
            self.assertEqual(net.num_bounded,num_bounded)

            # vargens
            for vargen in net.var_generators:
                if vargen.index % 3 == 0:
                    net.set_flags_of_component(vargen,
                                               'variable',
                                               ['active power','reactive power'])
                    net.set_flags_of_component(vargen,
                                               'variable',
                                               ['active power','reactive power'])
                    self.assertTrue(vargen.has_flags('variable','active power'))
                    self.assertTrue(vargen.has_flags('variable','reactive power'))
                    self.assertFalse(vargen.has_flags('fixed','active power'))
                    self.assertFalse(vargen.has_flags('fixed','reactive power'))
                    num_vars += 2
                else:
                    self.assertFalse(vargen.has_flags('variable','active power'))
                    self.assertFalse(vargen.has_flags('variable','reactive power'))
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)
            self.assertEqual(net.num_bounded,num_bounded)

            # batteries
            bcount = 0
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    net.set_flags_of_component(bat,
                                               'variable',
                                               'charging power')
                    bcount += 2
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    self.assertTrue(bat.has_flags('variable','charging power'))
                    self.assertFalse(bat.has_flags('fixed','charging power'))
                else:
                    self.assertFalse(bat.has_flags('variable','charging power'))
            num_vars += bcount
            self.assertEqual(net.num_vars,num_vars)
            self.assertEqual(net.num_fixed,num_fixed)

    def test_errors(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertRaises(KeyError,
                              net.set_flags,
                              -1,
                              'variable',
                              'any',
                              'voltage angle')
            self.assertRaises(KeyError,
                              net.set_flags,
                              'bus',
                              -1,
                              'any',
                              'voltage angle')
            self.assertRaises(KeyError,
                              net.set_flags,
                              'generator',
                              -1,
                              'any',
                              'active power')
            self.assertRaises(KeyError,
                              net.set_flags,
                              'branch',
                              -1,
                              'tap changer - v',
                              'tap ratio')

    def test_var_info_strings(self):

        T = 5

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case,num_periods=T)
            self.assertEqual(net.num_periods,T)

            # Add vargens and batteries
            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)
            net.add_batteries(net.get_generator_buses(),20.,50.)
            self.assertGreaterEqual(net.num_var_generators,1)
            self.assertGreaterEqual(net.num_batteries,1)

            # Errors
            for bus in net.buses:
                self.assertRaises(pf.BusError,bus.get_var_info_string,10)
            for branch in net.branches:
                self.assertRaises(pf.BranchError,branch.get_var_info_string,10)
            for gen in net.generators:
                self.assertRaises(pf.GeneratorError,gen.get_var_info_string,10)
            for load in net.loads:
                self.assertRaises(pf.LoadError,load.get_var_info_string,10)
            for shunt in net.shunts:
                self.assertRaises(pf.ShuntError,shunt.get_var_info_string,10)
            for vargen in net.var_generators:
                self.assertRaises(pf.VarGeneratorError,vargen.get_var_info_string,10)
            for bat in net.batteries:
                self.assertRaises(pf.BatteryError,bat.get_var_info_string,10)
            self.assertRaises(pf.NetworkError,net.get_var_info_string,10)
                
            # Set flags
            net.set_flags('bus','variable','any',['voltage magnitude','voltage angle'])
            net.set_flags('generator','variable','any',['active power','reactive power'])
            net.set_flags('load','variable','any',['active power','reactive power'])
            net.set_flags('variable generator','variable','any',['active power','reactive power'])
            net.set_flags('branch','variable','tap changer','tap ratio')
            net.set_flags('branch','variable','phase shifter','phase shift')
            net.set_flags('shunt','variable','switching - v','susceptance')
            net.set_flags('battery','variable','any',['charging power','energy level'])
            self.assertEqual(net.num_vars,
                             (net.num_buses*2+net.num_generators*2+net.num_loads*2+
                              net.num_var_generators*2+net.get_num_tap_changers()+
                              net.get_num_phase_shifters()+
                              net.get_num_switched_shunts()+
                              net.num_batteries*3)*net.num_periods)

            # Info strings
            for t in range(T):
                for bus in net.buses:
                    target = "bus:%d:voltage magnitude:%d" %(bus.index,t)
                    self.assertEqual(bus.get_var_info_string(bus.index_v_mag[t]),target)
                    self.assertEqual(net.get_var_info_string(bus.index_v_mag[t]),target)
                    target = "bus:%d:voltage angle:%d" %(bus.index,t)
                    self.assertEqual(bus.get_var_info_string(bus.index_v_ang[t]),target)
                    self.assertEqual(net.get_var_info_string(bus.index_v_ang[t]),target)
                for gen in net.generators:
                    target = "generator:%d:active power:%d" %(gen.index,t)
                    self.assertEqual(gen.get_var_info_string(gen.index_P[t]),target)
                    self.assertEqual(net.get_var_info_string(gen.index_P[t]),target)
                    target = "generator:%d:reactive power:%d" %(gen.index,t)
                    self.assertEqual(gen.get_var_info_string(gen.index_Q[t]),target)
                    self.assertEqual(net.get_var_info_string(gen.index_Q[t]),target)
                for branch in net.branches:
                    if branch.is_tap_changer():
                        target = "branch:%d:tap ratio:%d" %(branch.index,t)
                        self.assertEqual(branch.get_var_info_string(branch.index_ratio[t]),target)
                        self.assertEqual(net.get_var_info_string(branch.index_ratio[t]),target)
                    if branch.is_phase_shifter():
                        target = "branch:%d:phase shift:%d" %(branch.index,t)
                        self.assertEqual(branch.get_var_info_string(branch.index_phase[t]),target)
                        self.assertEqual(net.get_var_info_string(branch.index_phase[t]),target)
                for load in net.loads:
                    target = "load:%d:active power:%d" %(load.index,t)
                    self.assertEqual(load.get_var_info_string(load.index_P[t]),target)
                    self.assertEqual(net.get_var_info_string(load.index_P[t]),target)
                    target = "load:%d:reactive power:%d" %(load.index,t)
                    self.assertEqual(load.get_var_info_string(load.index_Q[t]),target)
                    self.assertEqual(net.get_var_info_string(load.index_Q[t]),target)
                for vargen in net.var_generators:
                    target = "variable generator:%d:active power:%d" %(vargen.index,t)
                    self.assertEqual(vargen.get_var_info_string(vargen.index_P[t]),target)
                    self.assertEqual(net.get_var_info_string(vargen.index_P[t]),target)
                    target = "variable generator:%d:reactive power:%d" %(vargen.index,t)
                    self.assertEqual(vargen.get_var_info_string(vargen.index_Q[t]),target)
                    self.assertEqual(net.get_var_info_string(vargen.index_Q[t]),target)
                for shunt in net.shunts:
                    if shunt.is_switched_v():
                        target = "shunt:%d:susceptance:%d" %(shunt.index,t)
                        self.assertEqual(shunt.get_var_info_string(shunt.index_b[t]),target)
                        self.assertEqual(net.get_var_info_string(shunt.index_b[t]),target)
                for bat in net.batteries:
                    target = "battery:%d:charging power:%d" %(bat.index,t)
                    self.assertEqual(bat.get_var_info_string(bat.index_Pc[t]),target)
                    self.assertEqual(net.get_var_info_string(bat.index_Pc[t]),target)
                    target = "battery:%d:discharging power:%d" %(bat.index,t)
                    self.assertEqual(bat.get_var_info_string(bat.index_Pd[t]),target)
                    self.assertEqual(net.get_var_info_string(bat.index_Pd[t]),target)
                    target = "battery:%d:energy level:%d" %(bat.index,t)
                    self.assertEqual(bat.get_var_info_string(bat.index_E[t]),target)
                    self.assertEqual(net.get_var_info_string(bat.index_E[t]),target)
            
    def tearDown(self):
        
        pass
