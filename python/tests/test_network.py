#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import math
import unittest
import numpy as np
import pfnet as pf
from . import test_cases
from scipy.sparse import coo_matrix, bmat, triu

class TestNetwork(unittest.TestCase):

    def setUp(self):

        # Networks
        self.T = 5

    def test_network(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Network()

            self.assertAlmostEqual(net.total_load_P/net.base_power,sum([l.P for l in net.loads]))

            net.clear_properties()

            self.assertEqual(net.bus_v_max,0.)
            self.assertEqual(net.bus_v_min,0.)
            self.assertEqual(net.bus_P_mis,0.)
            self.assertEqual(net.bus_Q_mis,0.)

            self.assertEqual(net.gen_P_cost,0.)
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.gen_P_vio,0.)

            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.tran_p_vio,0.)

            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)

            self.assertEqual(net.load_P_util,0.)
            self.assertEqual(net.load_P_vio,0.)

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertGreater(net.num_buses,0)
            self.assertGreater(net.num_generators,0)
            self.assertGreater(net.num_branches,0)
            self.assertGreater(net.num_loads,0)
            self.assertGreaterEqual(net.num_shunts,0)
            self.assertGreaterEqual(net.num_var_generators,0)
            self.assertGreaterEqual(net.num_batteries,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            self.assertRaises(pf.NetworkError,net.get_bus,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_bus,net.num_buses)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_generator,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_generator,net.num_generators)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_branch,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_branch,net.num_branches)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_shunt,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_shunt,net.num_shunts)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_load,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_load,net.num_loads)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_var_generator,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_var_generator,net.num_var_generators)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_battery,net.num_batteries)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_battery,-1)
            net.clear_error()

            # Counters
            self.assertEqual(net.get_num_P_adjust_gens(),
                             len([g for g in net.generators if g.P_min < g.P_max]))
            self.assertEqual(net.get_num_P_adjust_loads(),
                             len([l for l in net.loads if l.P_min < l.P_max]))

        # Multi period
        for case in test_cases.CASES:

            net = pf.Network(self.T)

            for t in range(net.num_periods):
                self.assertAlmostEqual(net.total_load_P[t]/net.base_power,sum([l.P[t] for l in net.loads]))

            net.clear_properties()

            self.assertEqual(net.num_periods,self.T)

            # prop type
            self.assertTrue(isinstance(net.bus_v_max,np.ndarray))
            self.assertTrue(isinstance(net.bus_v_min,np.ndarray))
            self.assertTrue(isinstance(net.bus_v_vio,np.ndarray))
            self.assertTrue(isinstance(net.bus_P_mis,np.ndarray))
            self.assertTrue(isinstance(net.bus_Q_mis,np.ndarray))
            self.assertTrue(isinstance(net.gen_P_cost,np.ndarray))
            self.assertTrue(isinstance(net.gen_v_dev,np.ndarray))
            self.assertTrue(isinstance(net.gen_Q_vio,np.ndarray))
            self.assertTrue(isinstance(net.gen_P_vio,np.ndarray))
            self.assertTrue(isinstance(net.tran_v_vio,np.ndarray))
            self.assertTrue(isinstance(net.tran_r_vio,np.ndarray))
            self.assertTrue(isinstance(net.tran_p_vio,np.ndarray))
            self.assertTrue(isinstance(net.shunt_v_vio,np.ndarray))
            self.assertTrue(isinstance(net.shunt_b_vio,np.ndarray))
            self.assertTrue(isinstance(net.load_P_util,np.ndarray))
            self.assertTrue(isinstance(net.load_P_vio,np.ndarray))
            self.assertTrue(isinstance(net.num_actions,np.ndarray))

            # prop shape
            self.assertTupleEqual(net.bus_v_max.shape,(self.T,))
            self.assertTupleEqual(net.bus_v_min.shape,(self.T,))
            self.assertTupleEqual(net.bus_v_vio.shape,(self.T,))
            self.assertTupleEqual(net.bus_P_mis.shape,(self.T,))
            self.assertTupleEqual(net.bus_Q_mis.shape,(self.T,))
            self.assertTupleEqual(net.gen_P_cost.shape,(self.T,))
            self.assertTupleEqual(net.gen_v_dev.shape,(self.T,))
            self.assertTupleEqual(net.gen_Q_vio.shape,(self.T,))
            self.assertTupleEqual(net.gen_P_vio.shape,(self.T,))
            self.assertTupleEqual(net.tran_v_vio.shape,(self.T,))
            self.assertTupleEqual(net.tran_r_vio.shape,(self.T,))
            self.assertTupleEqual(net.tran_p_vio.shape,(self.T,))
            self.assertTupleEqual(net.shunt_v_vio.shape,(self.T,))
            self.assertTupleEqual(net.shunt_b_vio.shape,(self.T,))
            self.assertTupleEqual(net.load_P_util.shape,(self.T,))
            self.assertTupleEqual(net.load_P_vio.shape,(self.T,))
            self.assertTupleEqual(net.num_actions.shape,(self.T,))

    def test_variables(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            num_so_far = 0

            # P adjust gens
            num_adj = 0
            for gen in net.generators:
                if gen.P_min < gen.P_max:
                    num_adj += 1

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.get_num_P_adjust_gens(),num_adj)

            net.set_flags('generator',
                          'variable',
                          'adjustable active power',
                          'active power')
            num_so_far += net.get_num_P_adjust_gens()

            self.assertEqual(net.num_vars,
                             num_so_far)

            # P adjust loads
            for load in net.loads:
                if load.index % 2 == 0:
                    load.P_min = 0
                    load.P_max = 10
                    self.assertEqual(load.P_min,0.)
                    self.assertEqual(load.P_max,10.)
            num_adj = 0
            for load in net.loads:
                if load.P_min < load.P_max:
                    num_adj += 1

            self.assertEqual(net.get_num_P_adjust_loads(),num_adj)

            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          'active power')
            num_so_far += net.get_num_P_adjust_loads()

            self.assertEqual(net.num_vars,num_so_far)

            for load in net.loads:
                if load.index % 2 == 0:
                    self.assertTrue(load.has_flags('variable','active power'))
                else:
                    self.assertFalse(load.has_flags('variable','active power'))

            # loads reactive power
            net.set_flags('load',
                          'variable',
                          'any',
                          'reactive power')

            num_so_far += net.num_loads

            self.assertEqual(net.num_vars,num_so_far)

            # batter charging
            net.set_flags('battery',
                          'variable',
                          'any',
                          'charging power')
            num_so_far += 2*net.num_batteries

            self.assertEqual(net.num_vars,
                             num_so_far)

            for bat in net.batteries:
                self.assertTrue(bat.has_flags('variable','charging power'))
                self.assertFalse(bat.has_flags('variable','energy level'))

            # batter energy
            net.set_flags('battery',
                          'variable',
                          'any',
                          'energy level')
            num_so_far += net.num_batteries

            self.assertEqual(net.num_vars,
                             num_so_far)

            for bat in net.batteries:
                self.assertTrue(bat.has_flags('variable','charging power'))
                self.assertTrue(bat.has_flags('variable','energy level'))

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # add vargens
            net.add_var_generators(net.get_load_buses(),100.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)

            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude',
                           'voltage angle',
                           'voltage magnitude deviation',
                           'voltage magnitude violation'])

            net.set_flags('branch',
                          'variable',
                          'any',
                          ['tap ratio',
                           'tap ratio deviation',
                           'phase shift'])

            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power',
                           'reactive power'])

            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power'])

            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power',
                           'reactive power'])

            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power',
                           'energy level'])

            net.set_flags('shunt',
                          'variable',
                          'any',
                          ['susceptance',
                           'susceptance deviation'])

            self.assertEqual(net.num_vars,self.T*(net.num_buses*6+
                                                  net.num_branches*4+
                                                  net.num_generators*2+
                                                  net.num_loads*1+
                                                  net.num_batteries*3+
                                                  net.num_var_generators*2+
                                                  net.num_shunts*3))

    def test_buses(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertGreater(net.num_buses,0)

            self.assertEqual(len(net.buses),net.num_buses)

            # Comparisons
            bus0 = net.get_bus(0)
            br0 = net.get_branch(0)
            self.assertEqual(bus0.index,br0.index)
            self.assertFalse(bus0 is br0)
            self.assertFalse(bus0 == br0)
            self.assertTrue(bus0 != br0)

            self.assertTupleEqual(tuple(net.buses),tuple([net.get_bus(i) for i in range(net.num_buses)]))

            for i in range(net.num_buses):

                bus = net.get_bus(i)
                same_bus = net.get_bus(i)

                self.assertEqual(bus.index,i)

                # obj type
                self.assertEqual(bus.obj_type,'bus')
                self.assertNotEqual(bus.obj_type,'unknown')

                # vmag vang set get
                bus.v_mag = 1.234567
                self.assertEqual(bus.v_mag,1.234567)
                bus.v_ang = 0.123456
                self.assertEqual(bus.v_ang,0.123456)

                # Comparisons
                self.assertFalse(bus is same_bus)
                self.assertTrue(bus == same_bus)
                self.assertFalse(bus != same_bus)
                if i > 0:
                    j = 0
                else:
                    j = net.num_buses-1
                if i != j:
                    other_bus = net.get_bus(j)
                    self.assertFalse(bus is other_bus)
                    self.assertFalse(bus == other_bus)
                    self.assertTrue(bus != other_bus)

                # hash table (numbers)
                self.assertEqual(bus.number,net.get_bus_by_number(bus.number).number)

                # hash table (names)
                self.assertEqual(bus.name,net.get_bus_by_name(bus.name).name)

                # values
                self.assertGreater(bus.number,0)
                self.assertTrue(isinstance(bus.name,str) or isinstance(bus.name,unicode))
                self.assertGreater(bus.v_mag,0)
                self.assertGreater(bus.v_set,0)
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
                self.assertEqual(bus.sens_v_mag_u_bound,0.)
                self.assertEqual(bus.sens_v_mag_l_bound,0.)
                self.assertEqual(bus.sens_v_ang_u_bound,0.)
                self.assertEqual(bus.sens_v_ang_l_bound,0.)
                self.assertEqual(bus.sens_v_reg_by_gen,0.)
                self.assertEqual(bus.sens_v_reg_by_tran,0.)
                self.assertEqual(bus.sens_v_reg_by_shunt,0.)
                self.assertTrue(isinstance(bus.reg_generators,list))
                self.assertTrue(isinstance(bus.reg_trans,list))
                self.assertTrue(isinstance(bus.reg_shunts,list))
                self.assertTrue(isinstance(bus.generators,list))

                # name
                self.assertTrue(isinstance(bus.name,str) or isinstance(bus.name,unicode))
                bus.name = "some bus"
                self.assertTrue(isinstance(bus.name,str) or isinstance(bus.name,unicode))
                self.assertEqual(bus.name,"some bus")
                self.assertEqual(bus.name,'some bus')

                # generators
                for g in bus.generators:
                    self.assertTrue(isinstance(g,pf.Generator))
                    self.assertEqual(bus.index,g.bus.index)
                    self.assertTrue(bus == g.bus)

                # Load
                for l in bus.loads:
                    self.assertTrue(isinstance(l,pf.Load))
                    self.assertEqual(bus.index,l.bus.index)
                    self.assertTrue(bus == l.bus)

                # branches
                for b in bus.branches:
                    self.assertTrue(isinstance(b,pf.Branch))
                    self.assertTrue(bus == b.bus_k or bus == b.bus_m)

                # total injections
                self.assertEqual(bus.get_total_gen_P(),sum([g.P for g in bus.generators],0))
                self.assertEqual(bus.get_total_gen_Q(),sum([g.Q for g in bus.generators],0))
                self.assertEqual(bus.get_total_load_P(),sum([l.P for l in bus.loads],0))
                self.assertEqual(bus.get_total_load_Q(),sum([l.Q for l in bus.loads],0))

                # regulated by gen
                if bus.is_regulated_by_gen():
                    self.assertGreater(len(bus.reg_generators),0)
                    for gen in bus.reg_generators:
                        self.assertTrue(gen.is_regulator())
                        self.assertEqual(gen.reg_bus.number,bus.number)

                else:
                    self.assertTrue(len(bus.reg_generators) == 0)

                # slack
                if bus.is_slack():
                    self.assertTrue(bus.is_regulated_by_gen())
                    self.assertGreater(len(bus.generators),0)
                    self.assertGreater(len(bus.reg_generators),0)
                    self.assertEqual(len(bus.generators),len(bus.reg_generators))
                    self.assertEqual(set([g.index for g in bus.generators]),set([g.index for g in bus.reg_generators]))
                    for gen in bus.generators:
                        self.assertTrue(gen.is_slack())
                        self.assertEqual(gen.bus.number,bus.number)
                        self.assertTrue(gen.is_regulator())

                # regulated by tran
                if bus.is_regulated_by_tran():
                    self.assertGreater(len(bus.reg_trans),0)
                    self.assertTrue(any([t.is_tap_changer_v() for t in bus.reg_trans]))
                    self.assertGreaterEqual(bus.v_max,bus.v_min)
                    for tran in bus.reg_trans:
                        self.assertTrue(tran.is_tap_changer_v())
                        self.assertEqual(tran.reg_bus.number,bus.number)
                        if bus.is_regulated_by_gen():
                            self.assertGreaterEqual(bus.v_set,bus.v_min) # gen control set point
                            self.assertLessEqual(bus.v_set,bus.v_max)    # is inside tran control range
                    for tran in bus.reg_trans:
                        self.assertEqual(bus.number,tran.reg_bus.number)
                        if bus.number == tran.bus_k.number: # reg bus in "k"/ from side -> neg sensitivity
                            self.assertFalse(tran.has_pos_ratio_v_sens())
                        elif bus.number == tran.bus_m.number: # reg bus in "m"/to side -> pos sensitivity
                            self.assertTrue(tran.has_pos_ratio_v_sens())

                # regulated by shunt
                if bus.is_regulated_by_shunt():
                    self.assertGreater(len(bus.reg_shunts),0)
                    self.assertGreaterEqual(bus.v_max,bus.v_min)
                    for shunt in bus.reg_shunts:
                        self.assertTrue(shunt.is_switched_v())
                        self.assertEqual(shunt.reg_bus.number,bus.number)
                        if bus.is_regulated_by_gen():
                            self.assertGreaterEqual(bus.v_set,bus.v_min) # gen control set point
                            self.assertLessEqual(bus.v_set,bus.v_max)    # is inside tran control range

                # branches
                self.assertTrue(isinstance(bus.branches_k,list))
                for br in bus.branches_k:
                    self.assertEqual(bus.number,br.bus_k.number)
                self.assertTrue(isinstance(bus.branches_m,list))
                for br in bus.branches_m:
                    self.assertEqual(bus.number,br.bus_m.number)
                self.assertGreater(len(bus.branches),0)
                self.assertEqual(len(bus.branches),len(bus.branches_k)+len(bus.branches_m))
                self.assertEqual(len(bus.branches),bus.degree)

                # vargens
                for vg in bus.var_generators:
                    self.assertEqual(bus,vg.bus)

                # batters
                for b in bus.batteries:
                    self.assertTrue(isinstance(b,pf.Battery))
                    self.assertEqual(bus,b.bus)

                # prices
                self.assertEqual(bus.price,0.)
                bus.price = bus.index*4.33
                self.assertEqual(bus.price,bus.index*4.33)

            # sum of degrees
            sum_deg = 0
            for i in range(net.num_buses):
                sum_deg += net.get_bus(i).degree
            self.assertEqual(sum_deg,2*net.num_branches)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for bus in net.buses:

                self.assertEqual(bus.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(bus.v_mag[t],bus.v_mag[0])
                    self.assertEqual(bus.v_ang[t],bus.v_ang[0])
                    self.assertEqual(bus.v_set[t],bus.v_set[0])
                    self.assertEqual(bus.price[t],bus.price[0])

                # Set
                x = np.random.randn(self.T)
                bus.v_mag = x
                for t in range(self.T):
                    self.assertEqual(bus.v_mag[t],x[t])
                x = np.random.randn(self.T)
                bus.v_ang = x
                for t in range(self.T):
                    self.assertEqual(bus.v_ang[t],x[t])
                x = np.random.randn(self.T)
                bus.price = x
                for t in range(self.T):
                    self.assertEqual(bus.price[t],x[t])

                # Set (attribute array)
                for t in range(self.T):
                    mag = np.random.randn()
                    bus.v_mag[t] = mag
                    self.assertEqual(mag,bus.v_mag[t])
                    p = np.random.randn()
                    bus.price[t] = p
                    self.assertEqual(p,bus.price[t])
                    ang = np.random.randn()
                    bus.v_ang[t] = ang
                    self.assertEqual(ang,bus.v_ang[t])

                    mag = np.random.randn()
                    ar = bus.v_mag
                    self.assertTupleEqual(ar.shape,(self.T,))
                    ar[t] = mag
                    self.assertEqual(ar[t],mag)
                    self.assertEqual(bus.v_mag[t],mag)

            # Indexing
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude',
                           'voltage angle',
                           'voltage magnitude deviation',
                           'voltage magnitude violation'])

            index = 0
            for bus in net.buses:
                self.assertTrue(np.all(bus.index_v_mag == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(bus.index_v_ang == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(bus.index_y == range(index,index+2*self.T,2)))
                self.assertTrue(np.all(bus.index_z == range(index+1,index+2*self.T,2)))
                index += 2*self.T
                self.assertTrue(np.all(bus.index_vl == range(index,index+2*self.T,2)))
                self.assertTrue(np.all(bus.index_vh == range(index+1,index+2*self.T,2)))
                index += 2*self.T

    def test_generators(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertTrue(net.num_generators > 0)

            self.assertEqual(len(net.generators),net.num_generators)

            for i in range(net.num_generators):

                gen = net.get_generator(i)
                same_gen = net.get_generator(i)

                self.assertEqual(gen.index,i)

                self.assertEqual(gen.index,net.generators[i].index)

                self.assertTrue(gen.bus)

                # Comparisons
                self.assertFalse(gen == net.get_bus(0))
                self.assertFalse(gen is same_gen)
                self.assertTrue(gen == same_gen)
                self.assertFalse(gen != same_gen)
                if i > 0:
                    j = 0
                else:
                    j = net.num_generators-1
                if i != j:
                    other_gen = net.get_generator(j)
                    self.assertFalse(gen is other_gen)
                    self.assertFalse(gen == other_gen)
                    self.assertTrue(gen != other_gen)

                # obj type
                self.assertEqual(gen.obj_type,'generator')
                self.assertNotEqual(gen.obj_type,'unknown')

                # slack
                if gen.is_slack():
                    self.assertTrue(gen.is_regulator())
                    self.assertTrue(gen.bus.is_slack())
                else:
                    self.assertFalse(gen.bus.is_slack())

                # regulator
                if gen.is_regulator():
                    self.assertGreater(gen.Q_max,gen.Q_min)
                    self.assertTrue(gen.reg_bus)
                    self.assertTrue(gen.reg_bus.is_regulated_by_gen())
                else:
                    self.assertRaises(pf.BusError,lambda : gen.reg_bus)

                # p adjustable
                if gen.P_min < gen.P_max:
                    self.assertTrue(gen.is_P_adjustable())
                else:
                    self.assertFalse(gen.is_P_adjustable())

                # set/get P_min and P_max
                self.assertNotEqual(gen.P_max,2*np.pi)
                self.assertNotEqual(gen.P_min,np.pi)
                gen.P_min = np.pi
                gen.P_max = 2*np.pi
                self.assertEqual(gen.P_min,np.pi)
                self.assertEqual(gen.P_max,2*np.pi)

                # set/get P,Q
                self.assertNotEqual(gen.P,0.333)
                self.assertNotEqual(gen.Q,0.222)
                gen.P = 0.333
                gen.Q = 0.222
                self.assertEqual(gen.P,0.333)
                self.assertEqual(gen.Q,0.222)

                # set/get cost coeffs
                if case.split('.')[-1] == 'raw':
                    self.assertEqual(gen.cost_coeff_Q0,0.)
                    self.assertEqual(gen.cost_coeff_Q1,2000.)
                    self.assertEqual(gen.cost_coeff_Q2,100.)
                gen.cost_coeff_Q0 = 1.4
                gen.cost_coeff_Q1 = 42.
                gen.cost_coeff_Q2 = 2.5
                self.assertEqual(gen.cost_coeff_Q0,1.4)
                self.assertEqual(gen.cost_coeff_Q1,42.)
                self.assertEqual(gen.cost_coeff_Q2,2.5)

                # P cost
                self.assertEqual(gen.P_cost,
                                 (gen.cost_coeff_Q0+
                                  gen.cost_coeff_Q1*gen.P+
                                  gen.cost_coeff_Q2*(gen.P**2.)))

                # sens
                self.assertEqual(gen.sens_P_u_bound,0.)
                self.assertEqual(gen.sens_P_l_bound,0.)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for gen in net.generators:

                self.assertEqual(gen.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(gen.P[t],gen.P[0])
                    self.assertEqual(gen.Q[t],gen.Q[0])

                # Set
                x = np.random.randn(self.T)
                gen.P = x
                for t in range(self.T):
                    self.assertEqual(gen.P[t],x[t])
                x = np.random.randn(self.T)
                gen.Q = x
                for t in range(self.T):
                    self.assertEqual(gen.Q[t],x[t])
                x = np.random.randn(self.T)

                # Set (attribute array)
                for t in range(self.T):
                    p = np.random.randn()
                    gen.P[t] = p
                    self.assertEqual(gen.P[t],p)
                    q = np.random.randn()
                    gen.Q[t] = q
                    self.assertEqual(gen.Q[t],q)

            # Indexing
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power',
                           'reactive power'])

            index = 0
            for gen in net.generators:
                self.assertTrue(np.all(gen.index_P == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(gen.index_Q == range(index,index+self.T)))
                index += self.T

    def test_branches(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertTrue(net.num_branches > 0)

            self.assertEqual(net.num_branches,len(net.branches))

            for i in range(net.num_branches):

                branch = net.get_branch(i)
                same_branch = net.get_branch(i)

                # Comparisons
                self.assertFalse(branch == net.get_bus(0))
                self.assertFalse(branch is same_branch)
                self.assertTrue(branch == same_branch)
                self.assertFalse(branch != same_branch)
                if i > 0:
                    j = 0
                else:
                    j = net.num_branches-1
                if i != j:
                    other_branch = net.get_branch(j)
                    self.assertFalse(branch is other_branch)
                    self.assertFalse(branch == other_branch)
                    self.assertTrue(branch != other_branch)

                # obj type
                self.assertEqual(branch.obj_type,'branch')
                self.assertNotEqual(branch.obj_type,'unknown')

                self.assertEqual(branch.index,net.branches[i].index)

                self.assertEqual(branch.index,i)

                self.assertTrue(branch.bus_k)
                self.assertTrue(branch.bus_m)
                self.assertGreater(branch.ratio,0)

                self.assertGreaterEqual(branch.ratingA,0.)
                self.assertGreaterEqual(branch.ratingB,0.)
                self.assertGreaterEqual(branch.ratingC,0.)
                r = np.random.rand()
                branch.ratingA = r
                self.assertEqual(r,branch.ratingA)
                r = np.random.rand()
                branch.ratingB = r
                self.assertEqual(r,branch.ratingB)
                r = np.random.rand()
                branch.ratingC = r
                self.assertEqual(r,branch.ratingC)

                # tap changer v
                if branch.is_tap_changer_v():
                    self.assertTrue(branch.is_tap_changer())
                    self.assertTrue(branch.reg_bus)
                    self.assertTrue(branch.index in [b.index for b in branch.reg_bus.reg_trans])
                    self.assertTrue(branch.reg_bus.is_regulated_by_tran())
                else:
                    self.assertRaises(pf.BusError,lambda : branch.reg_bus)

                # sens
                self.assertEqual(branch.sens_P_u_bound,0.)
                self.assertEqual(branch.sens_P_l_bound,0.)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for branch in net.branches:

                self.assertEqual(branch.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(branch.ratio[t],branch.ratio[0])
                    self.assertEqual(branch.phase[t],branch.phase[0])

            # Indexing
            net.set_flags('branch',
                          'variable',
                          'any',
                          ['tap ratio',
                           'phase shift',
                           'tap ratio deviation'])

            index = 0
            for branch in net.branches:
                self.assertTrue(np.all(branch.index_ratio == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(branch.index_phase == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(branch.index_ratio_y == range(index,index+2*self.T,2)))
                self.assertTrue(np.all(branch.index_ratio_z == range(index+1,index+2*self.T,2)))
                index += 2*self.T


    def test_branch_flows(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            # set some variables
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude',
                           'voltage angle',
                           'voltage magnitude deviation',
                           'voltage magnitude violation'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'variable',
                          'any',
                          ['tap ratio','tap ratio deviation','phase shift'])
            net.set_flags('shunt',
                          'variable',
                          'any',
                          ['susceptance','susceptance deviation'])

            x0 = net.get_var_values()
            xR = x0 + np.random.random(x0.size)/10.
            
            for branch in net.branches:

                # compute branch flows
                flows =  compute_branch_flows({'ratio' : branch.ratio,
                                               'phase' : branch.phase,
                                               'bus_k.v_mag' : branch.bus_k.v_mag,
                                               'bus_k.v_ang' : branch.bus_k.v_ang,
                                               'bus_m.v_mag' : branch.bus_m.v_mag,
                                               'bus_m.v_ang' : branch.bus_m.v_ang,
                                               'g' : branch.g,
                                               'g_k' : branch.g_k,
                                               'g_m' : branch.g_m,
                                               'b' : branch.b,
                                               'b_k' : branch.b_k,
                                               'b_m' : branch.b_m})

                self.assertAlmostEqual(branch.P_k_shunt, flows['P_k_sh'])
                self.assertAlmostEqual(branch.Q_k_shunt, flows['Q_k_sh'])
                self.assertAlmostEqual(branch.P_m_shunt, flows['P_m_sh'])
                self.assertAlmostEqual(branch.Q_m_shunt, flows['Q_m_sh'])
                self.assertAlmostEqual(branch.P_km_series, flows['P_km_ser'])
                self.assertAlmostEqual(branch.Q_km_series, flows['Q_km_ser'])
                self.assertAlmostEqual(branch.P_mk_series, flows['P_mk_ser'])
                self.assertAlmostEqual(branch.Q_mk_series, flows['Q_mk_ser'])
                self.assertAlmostEqual(branch.P_km, flows['P_km'])
                self.assertAlmostEqual(branch.Q_km, flows['Q_km'])
                self.assertAlmostEqual(branch.P_mk, flows['P_mk'])
                self.assertAlmostEqual(branch.Q_mk, flows['Q_mk'])

                # check flow at bus equal to shunt + series elements
                self.assertTrue(branch.P_km == branch.P_km_series+branch.P_k_shunt)
                self.assertTrue(branch.Q_km == branch.Q_km_series+branch.Q_k_shunt)
                self.assertTrue(branch.P_mk == branch.P_mk_series+branch.P_m_shunt)
                self.assertTrue(branch.Q_mk == branch.Q_mk_series+branch.Q_m_shunt)

                # check from-to matches k-m
                ###########################
                self.assertTrue(branch.P_km == branch.P_from_to)
                self.assertTrue(branch.Q_km == branch.Q_from_to)
                self.assertTrue(branch.P_mk == branch.P_to_from)
                self.assertTrue(branch.Q_mk == branch.Q_to_from)
                self.assertTrue(branch.P_km_series == branch.P_series_from_to)
                self.assertTrue(branch.Q_km_series == branch.Q_series_from_to)
                self.assertTrue(branch.P_mk_series == branch.P_series_to_from)
                self.assertTrue(branch.Q_mk_series == branch.Q_series_to_from)
                self.assertTrue(branch.P_k_shunt == branch.P_shunt_from)
                self.assertTrue(branch.Q_k_shunt == branch.Q_shunt_from)
                self.assertTrue(branch.P_m_shunt == branch.P_shunt_to)
                self.assertTrue(branch.Q_m_shunt == branch.Q_shunt_to)

                # check passing variables to calculate flows
                ############################################
                flowsR = compute_branch_flows({'ratio': xR[branch.index_ratio],
                                               'phase': xR[branch.index_phase],
                                               'bus_k.v_mag': xR[branch.bus_k.index_v_mag],
                                               'bus_k.v_ang': xR[branch.bus_k.index_v_ang],
                                               'bus_m.v_mag': xR[branch.bus_m.index_v_mag],
                                               'bus_m.v_ang': xR[branch.bus_m.index_v_ang],
                                               'g' : branch.g,
                                               'g_k' : branch.g_k,
                                               'g_m' : branch.g_m,
                                               'b' : branch.b,
                                               'b_k' : branch.b_k,
                                               'b_m' : branch.b_m})

                self.assertAlmostEqual(flowsR['P_km_ser'], branch.get_P_km_series(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_km_ser'], branch.get_Q_km_series(var_values=xR))
                self.assertAlmostEqual(flowsR['P_mk_ser'], branch.get_P_mk_series(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_mk_ser'], branch.get_Q_mk_series(var_values=xR))
                self.assertAlmostEqual(flowsR['P_k_sh'], branch.get_P_k_shunt(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_k_sh'], branch.get_Q_k_shunt(var_values=xR))
                self.assertAlmostEqual(flowsR['P_m_sh'], branch.get_P_m_shunt(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_m_sh'], branch.get_Q_m_shunt(var_values=xR))
                self.assertAlmostEqual(flowsR['P_km'], branch.get_P_km(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_km'], branch.get_Q_km(var_values=xR))
                self.assertAlmostEqual(flowsR['P_mk'], branch.get_P_mk(var_values=xR))
                self.assertAlmostEqual(flowsR['Q_mk'], branch.get_Q_mk(var_values=xR))

        # Multi-period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for branch in net.branches:

                self.assertEqual(branch.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(branch.P_km[t],branch.P_km[0])
                    self.assertEqual(branch.Q_km[t],branch.Q_km[0])
                    self.assertEqual(branch.P_mk[t],branch.P_mk[0])
                    self.assertEqual(branch.Q_mk[t],branch.Q_mk[0])

    def test_shunts(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertGreaterEqual(net.num_shunts,0)

            self.assertEqual(len(net.shunts),net.num_shunts)

            for i in range(net.num_shunts):

                shunt = net.get_shunt(i)

                # obj type
                self.assertEqual(shunt.obj_type,'shunt')
                self.assertNotEqual(shunt.obj_type,'unknown')

                self.assertEqual(shunt.index,i)

                self.assertEqual(shunt.index,net.shunts[i].index)

                self.assertTrue(shunt.bus)

                # switched shunt v
                if shunt.is_switched_v():
                    self.assertTrue(shunt.reg_bus)
                    self.assertGreater(shunt.b_max,shunt.b_min)
                    self.assertTrue(shunt.index in [s.index for s in shunt.reg_bus.reg_shunts])
                    self.assertTrue(shunt.reg_bus.is_regulated_by_shunt())

                # fixed
                else:
                    self.assertTrue(shunt.is_fixed())
                    self.assertRaises(pf.BusError,lambda : shunt.reg_bus)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for shunt in net.shunts:

                self.assertEqual(shunt.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(shunt.b[t],shunt.b[0])

            # Indexing
            net.set_flags('shunt',
                          'variable',
                          'any',
                          ['susceptance',
                           'susceptance deviation'])

            index = 0
            for shunt in net.shunts:
                self.assertTrue(np.all(shunt.index_b == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(shunt.index_y == range(index,index+2*self.T,2)))
                self.assertTrue(np.all(shunt.index_z == range(index+1,index+2*self.T,2)))
                index += 2*self.T

    def test_loads(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertTrue(net.num_loads > 0)

            self.assertEqual(net.num_loads,len(net.loads))

            self.assertEqual(net.num_loads,sum([len(b.loads) for b in net.buses]))

            for i in range(net.num_loads):

                load = net.get_load(i)

                # obj type
                self.assertEqual(load.obj_type,'load')
                self.assertNotEqual(load.obj_type,'unknown')

                self.assertEqual(load.index,i)

                self.assertEqual(load.index,net.loads[i].index)

                self.assertTrue(load.bus)

                # P limits
                self.assertEqual(load.P,load.P_max)
                self.assertEqual(load.P,load.P_min)
                load.P_min = -1.23
                load.P_max = 2123.
                self.assertEqual(load.P_min,-1.23)
                self.assertEqual(load.P_max,2123.)

                # P, Q
                load.P = 0.3241
                load.Q = 0.1212
                self.assertEqual(load.P,0.3241)
                self.assertEqual(load.Q,0.1212)

                # adjustable
                self.assertTrue(load.is_P_adjustable())
                load.P_min = 0.5
                load.P_max = 0.5
                self.assertFalse(load.is_P_adjustable())
                load.P_min = 1.
                load.P_max = -2.
                self.assertFalse(load.is_P_adjustable())

                # utiltiy
                if case.split('.')[-1] == 'raw':
                    self.assertEqual(load.util_coeff_Q0,0.)
                    self.assertEqual(load.util_coeff_Q1,20000.)
                    self.assertEqual(load.util_coeff_Q2,-100.)
                load.util_coeff_Q0 = 1.4
                load.util_coeff_Q1 = 42.
                load.util_coeff_Q2 = -2.5
                self.assertEqual(load.util_coeff_Q0,1.4)
                self.assertEqual(load.util_coeff_Q1,42.)
                self.assertEqual(load.util_coeff_Q2,-2.5)

                # P util
                self.assertEqual(load.P_util,
                                 (load.util_coeff_Q0+
                                  load.util_coeff_Q1*load.P+
                                  load.util_coeff_Q2*(load.P**2.)))
                self.assertEqual(load.P_util,
                                 (1.4+42.*load.P-2.5*(load.P**2.)))

                # sens
                self.assertEqual(load.sens_P_u_bound,0.)
                self.assertEqual(load.sens_P_l_bound,0.)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for load in net.loads:

                self.assertEqual(load.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(load.P[t],load.P[0])
                    self.assertEqual(load.Q[t],load.Q[0])

                # Set
                x = np.random.randn(self.T)
                load.P = x
                for t in range(self.T):
                    self.assertEqual(load.P[t],x[t])
                x = np.random.randn(self.T)
                load.Q = x
                for t in range(self.T):
                    self.assertEqual(load.Q[t],x[t])

                # Set (attribute array)
                for t in range(self.T):
                    p = np.random.randn()
                    load.P[t] = p
                    self.assertEqual(load.P[t],p)
                    q = np.random.randn()
                    load.Q[t] = q
                    self.assertEqual(load.Q[t],q)

            # Indexing
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power'])

            index = 0
            for load in net.loads:
                self.assertTrue(np.all(load.index_P == range(index,index+self.T)))
                index += self.T

    def test_var_generators(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(len(net.var_generators),net.num_var_generators)

            # existing vargens
            for vg in net.var_generators:
                self.assertTrue(isinstance(vg.name,str) or isinstance(vg.name,unicode))
                self.assertEqual(net.get_var_generator_by_name(vg.name).name,vg.name)
                vg.name = "some vargen"
                self.assertEqual(vg.name,"some vargen")
                self.assertRaises(pf.NetworkError,net.get_var_generator_by_name,vg.name)

            # add vargens
            load_buses = net.get_load_buses()
            self.assertEqual(len(load_buses),
                             len([b for b in net.buses if b.loads]))
            total_load = abs(sum([l.P for l in net.loads]))

            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)

            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                self.assertEqual(vargen.index,i)
                self.assertEqual(vargen.obj_type,'variable generator')
                self.assertNotEqual(vargen.obj_type,'unknown')
                self.assertTrue(isinstance(vargen.name,str) or isinstance(vargen.name,unicode))
                self.assertEqual(vargen.name,"VARGEN %d" %(vargen.index+1))
                self.assertEqual(net.get_var_generator_by_name(vargen.name).name,vargen.name)
                self.assertEqual(vargen.P,0.5*vargen.P_max)
                self.assertEqual(vargen.P_ava,0.5*vargen.P_max)
                self.assertEqual(vargen.P_max,0.8*total_load/net.num_var_generators)
                self.assertEqual(vargen.P_min,0.)
                self.assertEqual(vargen.P_std,0.3*vargen.P_max)
                self.assertEqual(vargen.index,i)
                self.assertEqual(vargen.Q,0.)
                self.assertEqual(vargen.Q_max,0.)
                self.assertEqual(vargen.Q_min,0.)
                self.assertEqual(len(vargen.bus.var_generators),1)
                self.assertEqual(vargen.bus.var_generators[0].index,vargen.index)
                vargen.P = np.pi
                self.assertEqual(vargen.P,np.pi)
                vargen.P_ava = 142.2123
                self.assertEqual(vargen.P_ava,142.2123)
                vargen.Q = (i+1)*12.5
                self.assertEqual(vargen.Q,(i+1)*12.5)
                vargen.Q_max = (i+1)*13.5
                self.assertEqual(vargen.Q_max,(i+1)*13.5)
                vargen.Q_min = (i+1)*11.5
                self.assertEqual(vargen.Q_min,(i+1)*11.5)
                self.assertLess(vargen.Q,vargen.Q_max)
                self.assertLess(vargen.Q_min,vargen.Q)
                vargen.Q = 0.
                vargen.Q_max = 0.
                vargen.Q_min = 0.

            # check buses
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                self.assertEqual(vargen.index,net.var_generators[i].index)
                self.assertTrue(vargen.bus)
                self.assertEqual(vargen.bus,load_buses[i])
                self.assertTrue(vargen.index in [vg.index for vg in load_buses[i].var_generators])
                
            # Set P,P_ava,P_max,P_min,P_std,Q,Q_max,Q_min
            self.assertGreater(net.num_var_generators,0)
            self.assertGreater(len(net.var_generators),0)
            for vg in net.var_generators:
                self.assertEqual(vg.P,np.pi)
                self.assertEqual(vg.P_ava,142.2123)
                self.assertEqual(vg.P_max,0.8*total_load/net.num_var_generators)
                self.assertEqual(vg.P_std,0.3*vg.P_max)
                self.assertEqual(vg.P_min,0.)
                self.assertEqual(vg.Q,0.)
                self.assertEqual(vg.Q_max,0.)
                self.assertEqual(vg.Q_min,0.)
                vg.P = 1.
                vg.P_ava = 1.234
                vg.P_min = 2.3
                vg.P_max = 2.
                vg.P_std = 3.
                vg.Q = 4.
                vg.Q_max = 5.
                vg.Q_min = 6
                self.assertEqual(vg.P_max,2.)
                self.assertEqual(vg.P_min,2.3)
                self.assertEqual(vg.P_std,3)
                self.assertEqual(vg.P,1.)
                self.assertEqual(vg.P_ava,1.234)
                self.assertEqual(vg.Q,4.)
                self.assertEqual(vg.Q_max,5.)
                self.assertEqual(vg.Q_min,6.)

        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(len(net.var_generators),net.num_var_generators)

            self.assertEqual(net.var_generators_corr_radius,1.)
            self.assertEqual(net.var_generators_corr_value,0.)

            # gen buses
            gen_buses = net.get_generator_buses()
            self.assertGreater(len(gen_buses),0)
            for b in gen_buses:
                self.assertTrue(b.generators is not None)
            self.assertEqual(len(gen_buses),
                             len([b for b in net.buses if b.generators]))

            # add vargens
            power_capacity = 80.
            power_base = 50.
            power_std = 30.
            corr_radius = 5
            corr_value = 0.05
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,-80.,50.,50.,5,0.05)
            self.assertTrue(net.has_error())
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,120.,50.,5,0.05)
            self.assertTrue(net.has_error())
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,-10,50.,5,0.05)
            self.assertTrue(net.has_error())
            net.clear_error()
            self.assertFalse(net.has_error())
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,50,-10.,5,0.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,50,50.,-1,0.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,50,50.,5,1.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,50,50.,5,-1.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_var_generators,gen_buses,80.,-10,50.,5,0.05)
            self.assertTrue(net.has_error())
            net.clear_error()
            net.add_var_generators(gen_buses,80.,50,50.,5,0.05)
            self.assertFalse(net.has_error())
            self.assertEqual(net.num_var_generators,len(gen_buses))
            net.add_var_generators([],power_capacity,power_base,power_std,corr_radius,corr_value)
            self.assertEqual(net.num_var_generators,0)
            self.assertEqual(net.var_generators,[])

            net.add_var_generators(gen_buses,power_capacity,power_base,power_std,corr_radius,corr_value)
            self.assertEqual(net.num_var_generators,len(gen_buses))
            self.assertEqual(len(net.var_generators),len(gen_buses))

            self.assertEqual(net.var_generators_corr_radius,corr_radius)
            self.assertEqual(net.var_generators_corr_value,corr_value)

            total_load = abs(sum([l.P for l in net.loads]))
            total_cap = sum([vg.P_max for vg in net.var_generators])

            self.assertLess(np.abs(0.8*total_load-total_cap),1e-10)
            self.assertLess(np.abs((power_base/100.)*power_capacity-100*sum([vg.P for vg in net.var_generators])/total_load),1e-10)

            for vg in net.var_generators:
                self.assertEqual(vg.P_min,0.)
                self.assertEqual(vg.P_max,(power_capacity/100.)*total_load/net.num_var_generators)
                self.assertEqual(vg.P,(power_base/100.)*vg.P_max)
                self.assertEqual(vg.P_ava,(power_base/100.)*vg.P_max)
                self.assertEqual(vg.P_std,(power_std/100.)*vg.P_max)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            net.add_var_generators(net.get_generator_buses(),80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)

            for vargen in net.var_generators:

                self.assertEqual(vargen.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(vargen.P[t],vargen.P[0])
                    self.assertEqual(vargen.P_std[t],vargen.P_std[0])
                    self.assertEqual(vargen.P_ava[t],vargen.P_ava[0])
                    self.assertEqual(vargen.Q[t],vargen.Q[0])

                # Set
                x = np.random.randn(self.T)
                vargen.P = x
                for t in range(self.T):
                    self.assertEqual(vargen.P[t],x[t])
                x = np.random.randn(self.T)
                vargen.P_ava = x
                for t in range(self.T):
                    self.assertEqual(vargen.P_ava[t],x[t])
                x = np.random.randn(self.T)
                vargen.Q = x
                for t in range(self.T):
                    self.assertEqual(vargen.Q[t],x[t])
                vargen.P_std = x
                for t in range(self.T):
                    self.assertEqual(vargen.P_std[t],x[t])

                # Set (attribute array)
                for t in range(self.T):
                    p = np.random.randn()
                    vargen.P[t] = p
                    self.assertEqual(vargen.P[t],p)
                    pava = np.random.randn()
                    vargen.P_ava[t] = pava
                    self.assertEqual(vargen.P_ava[t],pava)
                    pstd = np.random.randn()
                    vargen.P_std[t] = pstd
                    self.assertEqual(vargen.P_std[t],pstd)
                    q = np.random.randn()
                    vargen.Q[t] = q
                    self.assertEqual(vargen.Q[t],q)

            # Indexing
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power',
                           'reactive power'])

            index = 0
            for vargen in net.var_generators:
                self.assertTrue(np.all(vargen.index_P == range(index,index+self.T)))
                index += self.T
                self.assertTrue(np.all(vargen.index_Q == range(index,index+self.T)))
                index += self.T

    def test_batteries(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            
            # Add batteries
            gen_buses = net.get_generator_buses()
            self.assertGreater(len(gen_buses),0)
            net.add_batteries(gen_buses,20.,50.,0.9,0.8)

            self.assertEqual(net.num_batteries,len(gen_buses))
            self.assertEqual(net.num_batteries,len(net.batteries))
            self.assertEqual(net.num_batteries,sum([len(b.batteries) for b in net.buses]))
            
            for i in range(net.num_batteries):

                bat = net.get_battery(i)

                self.assertTrue(isinstance(bat,pf.Battery))

                # obj type
                self.assertEqual(bat.obj_type,'battery')
                self.assertNotEqual(bat.obj_type,'unknown')

                self.assertEqual(bat.index,i)
                self.assertEqual(bat.index,net.batteries[i].index)
                self.assertTrue(bat.bus)
                self.assertTrue(bat.index in map(lambda x: x.index,bat.bus.batteries))

                # Properties
                max_total_load = abs(net.total_load_P)/net.base_power
                self.assertAlmostEqual(bat.P_max,0.2*max_total_load/net.num_batteries)
                self.assertAlmostEqual(bat.P_min,-0.2*max_total_load/net.num_batteries)
                self.assertAlmostEqual(bat.E_max,0.5*max_total_load/net.num_batteries)
                self.assertEqual(bat.eta_c,0.9)
                self.assertEqual(bat.eta_d,0.8)
                self.assertEqual(bat.P,0.)
                self.assertEqual(bat.E_init,0.5*bat.E_max)
                self.assertEqual(bat.E_final,bat.E_init)
                self.assertEqual(bat.E,bat.E_init)

                # P
                bat.P_min = -1.23
                bat.P_max = 2123.
                bat.P = 0.3241
                self.assertEqual(bat.P_min,-1.23)
                self.assertEqual(bat.P_max,2123.)
                self.assertEqual(bat.P,0.3241)

                # E
                bat.E_max = 22.2
                bat.E_init = 0.22
                bat.E_final = 1.2321
                bat.E = 0.555
                self.assertEqual(bat.E_max,22.2)
                self.assertEqual(bat.E,0.555)
                self.assertEqual(bat.E_init,0.22)
                self.assertEqual(bat.E_final,1.2321)

                # eta
                bat.eta_c = 0.91
                bat.eta_d = 0.95
                self.assertEqual(bat.eta_c,0.91)
                self.assertEqual(bat.eta_d,0.95)

            # Add zero batteries
            self.assertEqual(net.num_batteries,len(gen_buses))
            net.add_batteries([],20.,50.,0.9,0.8)
            self.assertEqual(net.num_batteries,0)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            # Add batteries
            gen_buses = net.get_generator_buses()
            self.assertGreater(len(gen_buses),0)
            net.add_batteries(gen_buses,40.,80.,0.7,0.85)

            self.assertEqual(net.num_batteries,len(gen_buses))
            self.assertEqual(net.num_batteries,len(net.batteries))
            self.assertEqual(net.num_batteries,sum([len(b.batteries) for b in net.buses]))            

            for bat in net.batteries:

                self.assertEqual(bat.num_periods,self.T)

                # Propagation
                for t in range(1,self.T):
                    self.assertEqual(bat.P[t],bat.P[0])
                    self.assertEqual(bat.E[t],bat.E[0])

                # Properties
                max_total_load = np.max(np.abs(net.total_load_P)/net.base_power)
                self.assertAlmostEqual(bat.P_max,0.4*max_total_load/net.num_batteries)
                self.assertAlmostEqual(bat.P_min,-0.4*max_total_load/net.num_batteries)
                self.assertAlmostEqual(bat.E_max,0.8*max_total_load/net.num_batteries)
                self.assertEqual(bat.eta_c,0.7)
                self.assertEqual(bat.eta_d,0.85)
                self.assertEqual(bat.E_init,0.5*bat.E_max)
                self.assertEqual(bat.E_final,bat.E_init)
                for t in range(net.num_periods):
                    self.assertEqual(bat.P[t],0.)
                    self.assertEqual(bat.E[t],bat.E_init)

                # Set
                x = np.random.randn(self.T)
                bat.P = x
                for t in range(self.T):
                    self.assertEqual(bat.P[t],x[t])
                x = np.random.randn(self.T)
                bat.E = x
                for t in range(self.T):
                    self.assertEqual(bat.E[t],x[t])

                # Set (attribute array)
                for t in range(self.T):
                    p = np.random.randn()
                    bat.P[t] = p
                    self.assertEqual(bat.P[t],p)
                    e = np.random.randn()
                    bat.E[t] = e
                    self.assertEqual(bat.E[t],e)

            # Indexing
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power',
                           'energy level'])

            index = 0
            for bat in net.batteries:
                self.assertTrue(np.all(bat.index_Pc == range(index,index+2*self.T,2)))
                self.assertTrue(np.all(bat.index_Pd == range(index+1,index+2*self.T,2)))
                index += 2*self.T
                self.assertTrue(np.all(bat.index_E == range(index,index+self.T)))
                index += self.T

    def test_clear_flags(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'variable',
                          'any',
                          'tap ratio')
            net.set_flags('shunt',
                          'variable',
                          'any',
                          'susceptance')
            net.set_flags('battery',
                          'variable',
                          'any',
                          'energy level')

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags('bus',
                          'fixed',
                          'any',
                          'voltage magnitude')
            net.set_flags('generator',
                          'fixed',
                          'any',
                          'active power')
            net.set_flags('load',
                          'fixed',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'fixed',
                          'any',
                          'tap ratio')
            net.set_flags('shunt',
                          'fixed',
                          'any',
                          'susceptance')
            net.set_flags('battery',
                          'fixed',
                          'any',
                          'energy level')

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))
            self.assertEqual(net.num_fixed,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))
            self.assertEqual(net.num_bounded,0)

            net.set_flags('bus',
                          'bounded',
                          'any',
                          'voltage magnitude')
            net.set_flags('generator',
                          'bounded',
                          'any',
                          'active power')
            net.set_flags('load',
                          'bounded',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'bounded',
                          'any',
                          'tap ratio')
            net.set_flags('shunt',
                          'bounded',
                          'any',
                          'susceptance')
            net.set_flags('battery',
                          'bounded',
                          'any',
                          'energy level')

            self.assertGreater(net.num_vars,0)
            self.assertGreater(net.num_fixed,0)
            self.assertGreater(net.num_bounded,0)

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))
            self.assertEqual(net.num_fixed,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))
            self.assertEqual(net.num_bounded,
                             (net.num_buses+
                              net.num_generators+
                              net.num_loads*2+
                              net.num_branches+
                              net.num_shunts+
                              net.num_batteries))

            net.clear_flags()

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

    def test_properties(self):

        # Single and multi period
        for case in test_cases.CASES:

            net = pf.Network()
            netMP = pf.Network(self.T)

            net.clear_properties()
            netMP.clear_properties()

            self.assertEqual(net.bus_v_max,0.)
            self.assertEqual(net.bus_v_min,0.)
            self.assertEqual(net.bus_v_vio,0.)
            self.assertEqual(net.bus_P_mis,0.)
            self.assertEqual(net.bus_Q_mis,0.)
            self.assertEqual(net.gen_P_cost,0.)
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.gen_P_vio,0.)
            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.tran_p_vio,0.)
            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)
            self.assertEqual(net.load_P_util,0.)
            self.assertEqual(net.load_P_vio,0.)
            self.assertEqual(net.num_actions,0)

            self.assertEqual(netMP.bus_v_max.shape[0],self.T)
            self.assertEqual(netMP.bus_v_min.shape[0],self.T)
            self.assertEqual(netMP.bus_v_vio.shape[0],self.T)
            self.assertEqual(netMP.bus_P_mis.shape[0],self.T)
            self.assertEqual(netMP.bus_Q_mis.shape[0],self.T)
            self.assertEqual(netMP.gen_P_cost.shape[0],self.T)
            self.assertEqual(netMP.gen_v_dev.shape[0],self.T)
            self.assertEqual(netMP.gen_Q_vio.shape[0],self.T)
            self.assertEqual(netMP.gen_P_vio.shape[0],self.T)
            self.assertEqual(netMP.tran_v_vio.shape[0],self.T)
            self.assertEqual(netMP.tran_r_vio.shape[0],self.T)
            self.assertEqual(netMP.tran_p_vio.shape[0],self.T)
            self.assertEqual(netMP.shunt_v_vio.shape[0],self.T)
            self.assertEqual(netMP.shunt_b_vio.shape[0],self.T)
            self.assertEqual(netMP.load_P_util.shape[0],self.T)
            self.assertEqual(netMP.load_P_vio.shape[0],self.T)
            self.assertEqual(netMP.num_actions.shape[0],self.T)

            self.assertTrue(np.all(netMP.bus_v_max == 0))
            self.assertTrue(np.all(netMP.bus_v_min == 0))
            self.assertTrue(np.all(netMP.bus_v_vio == 0))
            self.assertTrue(np.all(netMP.bus_P_mis == 0))
            self.assertTrue(np.all(netMP.bus_Q_mis == 0))
            self.assertTrue(np.all(netMP.gen_P_cost == 0))
            self.assertTrue(np.all(netMP.gen_v_dev == 0))
            self.assertTrue(np.all(netMP.gen_Q_vio == 0))
            self.assertTrue(np.all(netMP.gen_P_vio == 0))
            self.assertTrue(np.all(netMP.tran_v_vio == 0))
            self.assertTrue(np.all(netMP.tran_r_vio == 0))
            self.assertTrue(np.all(netMP.tran_p_vio == 0))
            self.assertTrue(np.all(netMP.shunt_v_vio == 0))
            self.assertTrue(np.all(netMP.shunt_b_vio == 0))
            self.assertTrue(np.all(netMP.load_P_util == 0))
            self.assertTrue(np.all(netMP.load_P_vio == 0))
            self.assertTrue(np.all(netMP.num_actions == 0))

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            netMP = pf.Parser(case).parse(case,self.T)
            self.assertEqual(netMP.num_periods,self.T)

            # Add vargens
            net.add_var_generators(net.get_load_buses(),80.,50.,30.,5,0.05)
            netMP.add_var_generators(netMP.get_load_buses(),80.,50.,30.,5,0.05)
            for vargen in net.var_generators:
                self.assertTrue(isinstance(vargen.name,str) or isinstance(vargen.name,unicode))
                self.assertEqual(vargen.name,"VARGEN %d" %(vargen.index+1))
                self.assertEqual(net.get_var_generator_by_name(vargen.name).name,vargen.name)
                vargen.P = 1.
                vargen.Q = 2.
                self.assertGreater(len(vargen.bus.loads),0)
            for vargen in netMP.var_generators:
                vargen.P = np.ones(self.T)
                vargen.Q = 2*np.ones(self.T)
                self.assertGreater(len(vargen.bus.loads),0)

            net.update_properties()
            netMP.update_properties()

            self.assertEqual(net.num_vars,0)
            self.assertEqual(netMP.num_vars,0)

            bus_v_max = net.bus_v_max
            bus_P_mis = net.bus_P_mis

            net.clear_properties()
            net.update_properties()
            netMP.clear_properties()
            netMP.update_properties()

            self.assertEqual(bus_v_max,net.bus_v_max)
            self.assertEqual(bus_P_mis,net.bus_P_mis)
            self.assertTrue(np.all(bus_v_max == netMP.bus_v_max))
            self.assertTrue(np.all(bus_P_mis == netMP.bus_P_mis))

            x0 = net.get_var_values()
            self.assertEqual(x0.size,0)
            net.update_properties(x0)

            x0MP = netMP.get_var_values()
            self.assertEqual(x0MP.size,0)
            netMP.update_properties(x0MP)

            self.assertGreater(net.bus_v_max,0.)
            self.assertGreater(net.bus_v_min,0.)
            self.assertGreaterEqual(net.bus_v_vio,0.)
            self.assertGreater(net.bus_P_mis,0.)
            self.assertGreater(net.bus_Q_mis,0.)
            self.assertGreaterEqual(net.gen_P_cost,0.)
            self.assertGreaterEqual(net.gen_v_dev,0.)
            self.assertGreaterEqual(net.gen_Q_vio,0.)
            self.assertGreaterEqual(net.gen_P_vio,0.)
            self.assertGreaterEqual(net.tran_v_vio,0.)
            self.assertGreaterEqual(net.tran_r_vio,0.)
            self.assertGreaterEqual(net.tran_p_vio,0.)
            self.assertGreaterEqual(net.shunt_v_vio,0.)
            self.assertGreaterEqual(net.shunt_b_vio,0.)
            self.assertNotEqual(net.load_P_util,0.)
            self.assertGreaterEqual(net.load_P_vio,0.)
            self.assertGreaterEqual(net.num_actions,0.)

            self.assertEqual(net.bus_v_max,net.get_properties()['bus_v_max'])
            self.assertEqual(net.bus_v_min,net.get_properties()['bus_v_min'])
            self.assertEqual(net.bus_v_vio,net.get_properties()['bus_v_vio'])
            self.assertEqual(net.bus_P_mis,net.get_properties()['bus_P_mis'])
            self.assertEqual(net.bus_Q_mis,net.get_properties()['bus_Q_mis'])

            self.assertEqual(net.gen_P_cost,net.get_properties()['gen_P_cost'])
            self.assertEqual(net.gen_v_dev,net.get_properties()['gen_v_dev'])
            self.assertEqual(net.gen_Q_vio,net.get_properties()['gen_Q_vio'])
            self.assertEqual(net.gen_P_vio,net.get_properties()['gen_P_vio'])

            self.assertEqual(net.tran_v_vio,net.get_properties()['tran_v_vio'])
            self.assertEqual(net.tran_r_vio,net.get_properties()['tran_r_vio'])
            self.assertEqual(net.tran_p_vio,net.get_properties()['tran_p_vio'])

            self.assertEqual(net.shunt_v_vio,net.get_properties()['shunt_v_vio'])
            self.assertEqual(net.shunt_b_vio,net.get_properties()['shunt_b_vio'])

            self.assertEqual(net.load_P_util,net.get_properties()['load_P_util'])
            self.assertEqual(net.load_P_vio,net.get_properties()['load_P_vio'])

            self.assertEqual(net.num_actions,net.get_properties()['num_actions'])

            # Buses
            #######
            vmax = -np.inf
            vmin = np.inf
            vvio = 0.
            tvvio = 0.
            svvio = 0.
            vdev = 0.
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                vmax = np.maximum(bus.v_mag,vmax)
                vmin = np.minimum(bus.v_mag,vmin)
                dv = np.max([bus.v_mag-bus.v_max,bus.v_min-bus.v_mag,0.])
                if dv > vvio:
                    vvio = dv
                if bus.is_regulated_by_tran():
                    if dv > tvvio:
                        tvvio = dv
                if bus.is_regulated_by_shunt():
                    if dv > svvio:
                        svvio = dv
                if bus.is_regulated_by_gen():
                    if np.abs(bus.v_mag-bus.v_set) > vdev:
                        vdev = np.abs(bus.v_mag-bus.v_set)
            self.assertLess(abs(net.bus_v_max-vmax),1e-10)
            self.assertLess(abs(net.bus_v_min-vmin),1e-10)
            self.assertLess(abs(net.bus_v_vio-vvio),1e-10)
            self.assertLess(abs(net.tran_v_vio-tvvio),1e-10)
            self.assertLess(abs(net.shunt_v_vio-svvio),1e-10)
            self.assertLess(abs(net.gen_v_dev-vdev),1e-10)

            self.assertTrue(np.all(np.abs(netMP.bus_v_max-vmax) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.bus_v_min-vmin) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.bus_v_vio-vvio) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.tran_v_vio-tvvio) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.shunt_v_vio-svvio) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.gen_v_dev-vdev) < 1e-10))
            
            # Generators
            Pvio = 0
            Qvio = 0
            Pcost = 0
            for gen in net.generators:
                Pcost += gen.P_cost
                dP = np.max([gen.P-gen.P_max,gen.P_min-gen.P,0.])
                if dP > Pvio:
                    Pvio = dP
                if gen.is_regulator():
                    dQ = np.max([gen.Q-gen.Q_max,gen.Q_min-gen.Q,0.])
                    if dQ > Qvio:
                        Qvio = dQ
            self.assertLess(abs(net.gen_P_vio-Pvio*net.base_power),1e-10)
            self.assertLess(abs(net.gen_Q_vio-Qvio*net.base_power),1e-10)
            self.assertLess(abs(net.gen_P_cost-Pcost),1e-8)

            self.assertTrue(np.all(np.abs(netMP.gen_P_vio-Pvio*netMP.base_power) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.gen_Q_vio-Qvio*netMP.base_power) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.gen_P_cost-Pcost) < 1e-8))

            # Loads
            Pvio = 0
            Putil = 0
            for load in net.loads:
                Putil += load.P_util
                self.assertEqual(load.P_util,
                                 (load.util_coeff_Q0+
                                  load.util_coeff_Q1*load.P+
                                  load.util_coeff_Q2*load.P**2.))
                dP = np.max([load.P-load.P_max,load.P_min-load.P,0.])
                if dP > Pvio:
                    Pvio = dP
            self.assertLess(abs(net.load_P_vio-Pvio*net.base_power),1e-10)
            self.assertLess(abs(net.load_P_util-Putil),1e-7)

            self.assertTrue(np.all(np.abs(netMP.load_P_vio-Pvio*netMP.base_power) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.load_P_util-Putil) < 1e-7))

            # Transformers
            rvio = 0.
            pvio = 0.
            for branch in net.branches:
                if branch.is_tap_changer():
                    dr = np.max([branch.ratio-branch.ratio_max,
                                 branch.ratio_min-branch.ratio,
                                 0])
                    if dr > rvio:
                        rvio = dr
                if branch.is_phase_shifter():
                    dp = np.max([branch.phase-branch.phase_max,
                                 branch.phase_min-branch.phase,
                                 0])
                    if dp > pvio:
                        pvio = dp
            self.assertLess(np.abs(rvio-net.tran_r_vio),1e-10)
            self.assertLess(np.abs(pvio-net.tran_p_vio),1e-10)

            self.assertTrue(np.all(np.abs(netMP.tran_r_vio-rvio) < 1e-10))
            self.assertTrue(np.all(np.abs(netMP.tran_p_vio-pvio) < 1e-10))

            # Mismatches 1
            constr = pf.Constraint('AC power balance',net)
            constr.analyze()
            constr.eval(x0)
            f = constr.f
            dP_list = []
            dQ_list = []
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                dP = f[bus.index_P]
                dQ = f[bus.index_Q]
                dP_list.append(dP)
                dQ_list.append(dQ)
                self.assertLess(np.abs(dP-bus.P_mis),1e-10)
                self.assertLess(np.abs(dQ-bus.Q_mis),1e-10)
            self.assertLess(np.abs(net.bus_P_mis-np.max(np.abs(dP_list))*net.base_power),1e-10)
            self.assertLess(np.abs(net.bus_Q_mis-np.max(np.abs(dQ_list))*net.base_power),1e-10)

            constrMP = pf.Constraint('AC power balance',netMP)
            constrMP.analyze()
            constrMP.eval(x0MP)
            fMP = constrMP.f
            offset = 0
            self.assertEqual(fMP.shape[0],2*netMP.num_buses*netMP.num_periods)
            for t in range(self.T):
                dP_list = []
                dQ_list = []
                ft = fMP[offset:offset+2*netMP.num_buses]
                offset += 2*netMP.num_buses
                for i in range(netMP.num_buses):
                    bus = netMP.get_bus(i)
                    dP = ft[bus.index_P]
                    dQ = ft[bus.index_Q]
                    dP_list.append(dP)
                    dQ_list.append(dQ)
                    self.assertLess(np.abs(dP-bus.P_mis[t]),1e-10)
                    self.assertLess(np.abs(dQ-bus.Q_mis[t]),1e-10)
                self.assertLess(np.abs(netMP.bus_P_mis[t]-np.max(np.abs(dP_list))*netMP.base_power),1e-10)
                self.assertLess(np.abs(netMP.bus_Q_mis[t]-np.max(np.abs(dQ_list))*netMP.base_power),1e-10)

            # Mismatches 2
            for vargen in net.var_generators:
                vargen.P = 0.
                vargen.Q = 0.
            net.update_properties()
            fsaved = f.copy()
            constr.eval(x0)
            f = constr.f
            for vargen in net.var_generators:
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       f[vargen.bus.index_P]-1.),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       vargen.bus.P_mis-1.),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       f[vargen.bus.index_Q]-2.),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       vargen.bus.Q_mis-2.),1e-10)
            for vargen in net.var_generators:
                self.assertGreater(len(vargen.bus.loads),0)
                vargen.bus.loads[0].P = vargen.bus.loads[0].P - 1.
                vargen.bus.loads[0].Q = vargen.bus.loads[0].Q - 2.
            net.update_properties()
            constr.eval(x0)
            f = constr.f
            for vargen in net.var_generators:
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       f[vargen.bus.index_P]),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_P]-
                                       vargen.bus.P_mis),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       f[vargen.bus.index_Q]),1e-10)
                self.assertLess(np.abs(fsaved[vargen.bus.index_Q]-
                                       vargen.bus.Q_mis),1e-10)

            # Mismatches 2 (multiperiod)
            for vargen in netMP.var_generators:
                vargen.P = np.zeros(netMP.num_periods)
                vargen.Q = np.zeros(netMP.num_periods)
            netMP.update_properties()
            fsaved = fMP.copy()
            constrMP.eval(x0MP)
            f = constrMP.f
            n = netMP.num_buses
            for vargen in netMP.var_generators:
                for t in range(self.T):
                    self.assertLess(np.abs(fsaved[vargen.bus.index_P+t*2*n]-
                                           f[vargen.bus.index_P+t*2*n]-1.),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_P+t*2*n]-
                                           vargen.bus.P_mis[t]-1.),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_Q+t*2*n]-
                                           f[vargen.bus.index_Q+t*2*n]-2.),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_Q+t*2*n]-
                                           vargen.bus.Q_mis[t]-2.),1e-10)
            for vargen in netMP.var_generators:
                self.assertGreater(len(vargen.bus.loads),0)
                vargen.bus.loads[0].P = vargen.bus.loads[0].P - 1.
                vargen.bus.loads[0].Q = vargen.bus.loads[0].Q - 2.
            netMP.update_properties()
            constrMP.eval(x0MP)
            f = constrMP.f
            for vargen in netMP.var_generators:
                for t in range(self.T):
                    self.assertLess(np.abs(fsaved[vargen.bus.index_P+t*2*n]-
                                           f[vargen.bus.index_P+t*2*n]),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_P+t*2*n]-
                                           vargen.bus.P_mis[t]),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_Q+t*2*n]-
                                           f[vargen.bus.index_Q+t*2*n]),1e-10)
                    self.assertLess(np.abs(fsaved[vargen.bus.index_Q+t*2*n]-
                                           vargen.bus.Q_mis[t]),1e-10)

            net.clear_properties()
            netMP.clear_properties()

            self.assertEqual(net.bus_v_max,0.)
            self.assertEqual(net.bus_v_min,0.)
            self.assertEqual(net.bus_v_vio,0.)
            self.assertEqual(net.bus_P_mis,0.)
            self.assertEqual(net.bus_Q_mis,0.)
            self.assertEqual(net.gen_P_cost,0.)
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.gen_P_vio,0.)
            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)
            self.assertEqual(net.load_P_util,0.)
            self.assertEqual(net.load_P_vio,0.)
            self.assertEqual(net.num_actions,0.)

            self.assertTrue(np.all(netMP.bus_v_max == 0))
            self.assertTrue(np.all(netMP.bus_v_min == 0))
            self.assertTrue(np.all(netMP.bus_v_vio == 0))
            self.assertTrue(np.all(netMP.bus_P_mis == 0))
            self.assertTrue(np.all(netMP.bus_Q_mis == 0))
            self.assertTrue(np.all(netMP.gen_P_cost == 0))
            self.assertTrue(np.all(netMP.gen_v_dev == 0))
            self.assertTrue(np.all(netMP.gen_Q_vio == 0))
            self.assertTrue(np.all(netMP.gen_P_vio == 0))
            self.assertTrue(np.all(netMP.tran_v_vio == 0))
            self.assertTrue(np.all(netMP.tran_r_vio == 0))
            self.assertTrue(np.all(netMP.tran_p_vio == 0))
            self.assertTrue(np.all(netMP.shunt_v_vio == 0))
            self.assertTrue(np.all(netMP.shunt_b_vio == 0))
            self.assertTrue(np.all(netMP.load_P_util == 0))
            self.assertTrue(np.all(netMP.load_P_vio == 0))
            self.assertTrue(np.all(netMP.num_actions == 0))

    def test_bus_sorting(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_bounded,0)

            # Variables
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])

            # Bounds
            net.set_flags('bus',
                          'bounded',
                          'any',
                          'voltage magnitude')

            x0 = net.get_var_values()

            self.assertEqual(net.num_vars,2*net.num_buses)
            self.assertEqual(net.num_bounded,net.num_buses)

            # Store random bus sensitivities
            constr = [pf.Constraint('AC power balance',net),
                      pf.Constraint('variable nonlinear bounds',net),
                      pf.Constraint('voltage regulation by generators',net),
                      pf.Constraint('voltage regulation by transformers',net),
                      pf.Constraint('voltage regulation by shunts',net)]
            for c in constr:
                c.analyze()
                c.eval(x0)
                c.store_sensitivities(None,np.random.randn(c.f.size),None,None)

            # Check bus largest mis and sens
            sens_types = [pf.BUS_SENS_P_BALANCE,
                          pf.BUS_SENS_Q_BALANCE,
                          pf.BUS_SENS_V_MAG_U_BOUND,
                          pf.BUS_SENS_V_MAG_L_BOUND,
                          pf.BUS_SENS_V_REG_BY_GEN,
                          pf.BUS_SENS_V_REG_BY_TRAN,
                          pf.BUS_SENS_V_REG_BY_SHUNT]
            mis_types = [pf.BUS_MIS_ACTIVE,
                         pf.BUS_MIS_REACTIVE]
            self.assertEqual(len(sens_types),len(set(sens_types)))
            self.assertEqual(len(mis_types),len(set(mis_types)))
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                sens = [abs(bus.sens_P_balance),
                        abs(bus.sens_Q_balance),
                        abs(bus.sens_v_mag_u_bound),
                        abs(bus.sens_v_mag_l_bound),
                        abs(bus.sens_v_reg_by_gen),
                        abs(bus.sens_v_reg_by_tran),
                        abs(bus.sens_v_reg_by_shunt)]
                sensm = max(sens)
                senst = sens_types[np.argmax(sens)]
                self.assertGreater(sensm,0)
                self.assertEqual(abs(bus.get_largest_sens()),sensm)
                self.assertEqual(bus.get_largest_sens_type(),senst)
                mis = [abs(bus.P_mis),
                       abs(bus.Q_mis)]
                mism = max(mis)
                mist = mis_types[np.argmax(mis)]
                self.assertEqual(abs(bus.get_largest_mis()),mism)
                if mis[0] != mis[1]:
                    self.assertEqual(bus.get_largest_mis_type(),mist)

            # Check bus quantities
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_LARGEST),
                                 bus.get_largest_sens())
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_P_BALANCE),
                                 bus.sens_P_balance)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_Q_BALANCE),
                                 bus.sens_Q_balance)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_V_MAG_U_BOUND),
                                 bus.sens_v_mag_u_bound)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_V_MAG_L_BOUND),
                                 bus.sens_v_mag_l_bound)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_V_REG_BY_GEN),
                                 bus.sens_v_reg_by_gen)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_V_REG_BY_TRAN),
                                 bus.sens_v_reg_by_tran)
                self.assertEqual(bus.get_quantity(pf.BUS_SENS_V_REG_BY_SHUNT),
                                 bus.sens_v_reg_by_shunt)
                self.assertEqual(bus.get_quantity(pf.BUS_MIS_LARGEST),
                                 bus.get_largest_mis())
                self.assertEqual(bus.get_quantity(pf.BUS_MIS_ACTIVE),
                                 bus.P_mis)
                self.assertEqual(bus.get_quantity(pf.BUS_MIS_REACTIVE),
                                 bus.Q_mis)
                self.assertEqual(bus.get_quantity(-1),
                                 0.)

            # Sort by largest mis
            bus_list = net.create_sorted_bus_list(pf.BUS_MIS_LARGEST)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.get_largest_mis()) >= abs(bus2.get_largest_mis()))
            self.assertTrue(all(r1))

            # Sort by P mismatch
            bus_list = net.create_sorted_bus_list(pf.BUS_MIS_ACTIVE)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.P_mis) >= abs(bus2.P_mis))
            self.assertTrue(all(r1))

            # Sort by largest sensitivity
            bus_list = net.create_sorted_bus_list(pf.BUS_SENS_LARGEST)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.get_largest_sens()) >= abs(bus2.get_largest_sens()))
            self.assertTrue(all(r1))

            # Sort by v_reg_by_gen sensitivity
            bus_list = net.create_sorted_bus_list(pf.BUS_SENS_V_REG_BY_GEN)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.sens_v_reg_by_gen) >= abs(bus2.sens_v_reg_by_gen))
            self.assertTrue(all(r1))

    def test_set_points(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            net.update_set_points()

            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertEqual(bus.v_mag,bus.v_set)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            for bus in net.buses:
                bus.v_mag = np.random.randn(self.T)

            net.update_set_points()

            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertTrue(np.all(bus.v_mag == bus.v_set))
                    self.assertTupleEqual(bus.v_mag.shape,(self.T,))
                    self.assertTupleEqual(bus.v_set.shape,(self.T,))

    def test_projections(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))

            # Set load Qs
            for load in net.loads:
                load.Q = load.index*3.3

            # bus vmag and vang
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])

            # gen powers
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')

            # load powers
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # branch ratio and phase
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')

            # shunt
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')

            # vargens
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # batteries
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])

            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-1) +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers_v() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              2*net.get_num_var_generators()+
                              2*net.num_loads+
                              3*net.num_batteries))

            # set vargens
            for vargen in net.var_generators:
                vargen.P = np.pi*vargen.index
                vargen.Q = -np.pi*vargen.index

            # set batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P *= -1.

            # var values (current values)
            x = net.get_var_values()

            # bus vmag
            P = net.get_var_projection('bus','voltage magnitude')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_buses-1)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_buses-1)
            vmag = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags('variable','voltage magnitude'):
                    self.assertEqual(vmag[index],bus.v_mag)
                    index += 1

            # bus vang
            P = net.get_var_projection('bus','voltage angle')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_buses-1)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_buses-1)
            vang = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags('variable','voltage angle'):
                    self.assertEqual(vang[index],bus.v_ang)
                    index += 1

            # gen active power
            P = net.get_var_projection('generator','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_slack_gens())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_slack_gens())
            gP = P*x
            index = 0
            for i in range(net.num_generators):
                gen = net.get_generator(i)
                if gen.has_flags('variable','active power'):
                    self.assertEqual(gP[index],gen.P)
                    index += 1

            # gen reactive power
            P = net.get_var_projection('generator','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_reg_gens())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_reg_gens())
            gQ = P*x
            index = 0
            for i in range(net.num_generators):
                gen = net.get_generator(i)
                if gen.has_flags('variable','reactive power'):
                    self.assertEqual(gQ[index],gen.Q)
                    index += 1

            # load active power
            P = net.get_var_projection('load','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_loads)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_loads)
            gP = P*x
            index = 0
            for i in range(net.num_loads):
                load = net.get_load(i)
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertEqual(gP[index],load.P)
                index += 1

            # load reactive power
            P = net.get_var_projection('load','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_loads)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_loads)
            gQ = P*x
            index = 0
            for i in range(net.num_loads):
                load = net.get_load(i)
                self.assertTrue(load.has_flags('variable','reactive power'))
                self.assertEqual(gQ[index],load.Q)
                self.assertEqual(gQ[index],load.index*3.3)
                index += 1

            # tap changer ratio
            P = net.get_var_projection('branch','tap ratio')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_tap_changers_v())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_tap_changers_v())
            bR = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags('variable','tap ratio'):
                    self.assertEqual(bR[index],br.ratio)
                    index += 1

            # phase shifter
            P = net.get_var_projection('branch','phase shift')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_phase_shifters())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_phase_shifters())
            bP = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags('variable','phase shift'):
                    self.assertEqual(bP[index],br.phase)
                    index += 1

            # shunt susceptance
            P = net.get_var_projection('shunt','susceptance')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_switched_shunts())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_switched_shunts())
            sS = P*x
            index = 0
            for i in range(net.num_shunts):
                shunt = net.get_shunt(i)
                if shunt.has_flags('variable','susceptance'):
                    self.assertEqual(sS[index],shunt.b)
                    index += 1

            # vargen active power
            P = net.get_var_projection('variable generator','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_var_generators)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_var_generators)
            vgP = P*x
            index = 0
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                self.assertEqual(vargen.index_P,vargen.index_Q-1)
                if vargen.has_flags('variable','active power'):
                    self.assertEqual(vgP[index],vargen.P)
                    self.assertEqual(vgP[index],vargen.index*np.pi)
                    index += 1

            # vargen reactive power
            P = net.get_var_projection('variable generator','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_var_generators)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_var_generators)
            vgQ = P*x
            index = 0
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                self.assertEqual(vargen.index_P+1,vargen.index_Q)
                if vargen.has_flags('variable','reactive power'):
                    self.assertEqual(vgQ[index],vargen.Q)
                    self.assertEqual(vgQ[index],-vargen.index*np.pi)
                    index += 1

            # battery charging
            P = net.get_var_projection('battery','charging power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],2*net.num_batteries)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,2*net.num_batteries)
            batP = P*x
            index = 0
            for i in range(net.num_batteries):
                bat = net.get_battery(i)
                self.assertEqual(bat.index_Pc,bat.index_Pd-1)
                self.assertEqual(bat.index_Pd,bat.index_E-1)
                ac = np.where(P.col == bat.index_Pc)[0]
                ad = np.where(P.col == bat.index_Pd)[0]
                self.assertEqual(ac.size,1)
                self.assertEqual(ad.size,1)
                self.assertEqual(P.col[ac[0]],bat.index_Pc)
                self.assertEqual(P.row[ac[0]],P.row[ad[0]]-1)
                if bat.has_flags('variable','charging power'):
                    if bat.P >= 0:
                        self.assertEqual(batP[index],bat.P)
                        self.assertEqual(batP[index+1],0.)
                    else:
                        self.assertEqual(batP[index],0.)
                        self.assertEqual(batP[index+1],-bat.P)
                    index += 2

            # battery energy
            P = net.get_var_projection('battery','energy level')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_batteries)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_batteries)
            batE = P*x
            index = 0
            for i in range(net.num_batteries):
                bat = net.get_battery(i)
                self.assertEqual(bat.index_Pc,bat.index_Pd-1)
                self.assertEqual(bat.index_Pd,bat.index_E-1)
                if bat.has_flags('variable','energy level'):
                    self.assertEqual(batE[index],bat.E)
                    index += 1

            # All
            Plist = [net.get_var_projection('bus','voltage magnitude'),
                     net.get_var_projection('bus','voltage angle'),
                     net.get_var_projection('generator','active power'),
                     net.get_var_projection('generator','reactive power'),
                     net.get_var_projection('load','active power'),
                     net.get_var_projection('load','reactive power'),
                     net.get_var_projection('branch','tap ratio'),
                     net.get_var_projection('branch','phase shift'),
                     net.get_var_projection('shunt','susceptance'),
                     net.get_var_projection('variable generator','active power'),
                     net.get_var_projection('variable generator','reactive power'),
                     net.get_var_projection('battery','charging power'),
                     net.get_var_projection('battery','energy level')]
            P = bmat([[P] for P in Plist if P.shape[0] > 0])
            self.assertTupleEqual(P.shape,(net.num_vars,net.num_vars))
            for i in range(10):
                x = np.random.randn(net.num_vars)
                self.assertLess(np.linalg.norm(x-P.T*P*x),1e-12)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))

            # Set load Qs
            for load in net.loads:
                load.Q = load.index*3.3*np.ones(net.num_periods)

            # bus vmag and vang
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])

            # gen powers
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')

            # load powers
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # branch ratio and phase
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')

            # shunt
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')

            # vargens
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # batteries
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])

            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-1) +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers_v() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              2*net.get_num_var_generators()+
                              2*net.num_loads+
                              3*net.num_batteries)*self.T)

            # set vargens
            for vargen in net.var_generators:
                vargen.P = np.pi*vargen.index*np.random.rand(self.T)
                vargen.Q = -np.pi*vargen.index*np.random.rand(self.T)

            # set batteries
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    bat.P = -1.*bat.P

            # var values (current values)
            x = net.get_var_values()

            # bus vmag
            P = net.get_var_projection('bus','voltage magnitude')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],(net.num_buses-1)*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,(net.num_buses-1)*self.T)
            vmag = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags('variable','voltage magnitude'):
                    for t in range(self.T):
                        self.assertEqual(vmag[index],bus.v_mag[t])
                        index += 1

            # bus vang
            P = net.get_var_projection('bus','voltage angle')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],(net.num_buses-1)*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,(net.num_buses-1)*self.T)
            vang = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags('variable','voltage angle'):
                    for t in range(self.T):
                        self.assertEqual(vang[index],bus.v_ang[t])
                        index += 1

            # gen active power
            P = net.get_var_projection('generator','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_slack_gens()*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_slack_gens()*self.T)
            gP = P*x
            index = 0
            for i in range(net.num_generators):
                gen = net.get_generator(i)
                if gen.has_flags('variable','active power'):
                    for t in range(self.T):
                        self.assertEqual(gP[index],gen.P[t])
                        index += 1

            # gen reactive power
            P = net.get_var_projection('generator','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_reg_gens()*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_reg_gens()*self.T)
            gQ = P*x
            index = 0
            for i in range(net.num_generators):
                gen = net.get_generator(i)
                if gen.has_flags('variable','reactive power'):
                    for t in range(self.T):
                        self.assertEqual(gQ[index],gen.Q[t])
                        index += 1

            # load active power
            P = net.get_var_projection('load','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_loads*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_loads*self.T)
            gP = P*x
            index = 0
            for i in range(net.num_loads):
                load = net.get_load(i)
                self.assertTrue(load.has_flags('variable','active power'))
                for t in range(self.T):
                    self.assertEqual(gP[index],load.P[t])
                    index += 1
                    
            # load reactive power
            P = net.get_var_projection('load','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_loads*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_loads*self.T)
            gQ = P*x
            index = 0
            for i in range(net.num_loads):
                load = net.get_load(i)
                self.assertTrue(load.has_flags('variable','reactive power'))
                for t in range(self.T):
                    self.assertEqual(gQ[index],load.Q[t])
                    self.assertEqual(gQ[index],load.index*3.3)
                    index += 1

            # tap changer ratio
            P = net.get_var_projection('branch','tap ratio')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_tap_changers_v()*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_tap_changers_v()*self.T)
            bR = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags('variable','tap ratio'):
                    for t in range(self.T):
                        self.assertEqual(bR[index],br.ratio[t])
                        index += 1

            # phase shifter
            P = net.get_var_projection('branch','phase shift')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_phase_shifters()*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_phase_shifters()*self.T)
            bP = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags('variable','phase shift'):
                    for t in range(self.T):
                        self.assertEqual(bP[index],br.phase[t])
                        index += 1

            # shunt susceptance
            P = net.get_var_projection('shunt','susceptance')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_switched_shunts()*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_switched_shunts()*self.T)
            sS = P*x
            index = 0
            for i in range(net.num_shunts):
                shunt = net.get_shunt(i)
                if shunt.has_flags('variable','susceptance'):
                    for t in range(self.T):
                        self.assertEqual(sS[index],shunt.b[t])
                        index += 1

            # vargen active power
            P = net.get_var_projection('variable generator','active power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_var_generators*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_var_generators*self.T)
            vgP = P*x
            index = 0
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                for t in range(self.T):
                    self.assertEqual(vargen.index_P[t],vargen.index_Q[t]-self.T)
                    if vargen.has_flags('variable','active power'):
                        self.assertEqual(vgP[index],vargen.P[t])
                        index += 1

            # vargen reactive power
            P = net.get_var_projection('variable generator','reactive power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_var_generators*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_var_generators*self.T)
            vgQ = P*x
            index = 0
            for i in range(net.num_var_generators):
                vargen = net.get_var_generator(i)
                for t in range(self.T):
                    self.assertEqual(vargen.index_P[t]+self.T,vargen.index_Q[t])
                    if vargen.has_flags('variable','reactive power'):
                        self.assertEqual(vgQ[index],vargen.Q[t])
                        index += 1

            # battery charging
            P = net.get_var_projection('battery','charging power')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],2*net.num_batteries*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,2*net.num_batteries*self.T)
            batP = P*x
            index = 0
            for i in range(net.num_batteries):
                bat = net.get_battery(i)
                for t in range(self.T):
                    self.assertEqual(bat.index_Pc[t],bat.index_Pd[t]-1)
                    self.assertEqual(bat.index_Pd[t],bat.index_E[t]-t-2*self.T+2*t+1)
                    ac = np.where(P.col == bat.index_Pc[t])[0]
                    ad = np.where(P.col == bat.index_Pd[t])[0]
                    self.assertEqual(ac.size,1)
                    self.assertEqual(ad.size,1)
                    self.assertEqual(P.col[ac[0]],bat.index_Pc[t])
                    self.assertEqual(P.row[ac[0]],P.row[ad[0]]-1)
                    if bat.has_flags('variable','charging power'):
                        if bat.P[t] >= 0:
                            self.assertEqual(batP[index],bat.P[t])
                            self.assertEqual(batP[index+1],0.)
                        else:
                            self.assertEqual(batP[index],0.)
                            self.assertEqual(batP[index+1],-bat.P[t])
                        index += 2

            # battery energy
            P = net.get_var_projection('battery','energy level')
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_batteries*self.T)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_batteries*self.T)
            batE = P*x
            index = 0
            for i in range(net.num_batteries):
                bat = net.get_battery(i)
                for t in range(self.T):
                    if bat.has_flags('variable','energy level'):
                        self.assertEqual(batE[index],bat.E[t])
                        index += 1

            # All
            Plist = [net.get_var_projection('bus','voltage magnitude'),
                     net.get_var_projection('bus','voltage angle'),
                     net.get_var_projection('generator','active power'),
                     net.get_var_projection('generator','reactive power'),
                     net.get_var_projection('load','active power'),
                     net.get_var_projection('load','reactive power'),
                     net.get_var_projection('branch','tap ratio'),
                     net.get_var_projection('branch','phase shift'),
                     net.get_var_projection('shunt','susceptance'),
                     net.get_var_projection('variable generator','active power'),
                     net.get_var_projection('variable generator','reactive power'),
                     net.get_var_projection('battery','charging power'),
                     net.get_var_projection('battery','energy level')]
            P = bmat([[P] for P in Plist if P.shape[0] > 0])
            self.assertTupleEqual(P.shape,(net.num_vars,net.num_vars))
            for i in range(10):
                x = np.random.randn(net.num_vars)
                self.assertLess(np.linalg.norm(x-P.T*P*x),1e-12)

    def test_fancy_projections(self):

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))

            # bus vmag and vang
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])

            # gen powers
            net.set_flags('generator',
                          'variable',
                          'not slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')

            # load powers
            net.set_flags('load',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # branch ratio and phase
            net.set_flags('branch',
                          'variable',
                          'tap changer - v',
                          'tap ratio')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')

            # shunt
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')

            # vargens
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])

            # batteries
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])

            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-1) +
                              net.num_generators-net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers_v() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              2*net.get_num_var_generators()+
                              2*net.num_loads+
                              3*net.num_batteries)*self.T)

            self.assertRaises(KeyError,net.get_var_projection,'all','voltage magnitude',2,4)
            self.assertFalse(net.has_error())
            net.clear_error()
            self.assertFalse(net.has_error())

            # bus all
            P = net.get_var_projection('bus','all',2,3)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(2*(net.num_buses-1)*2,net.num_vars))

            # gen all
            P = net.get_var_projection('generator','all',2,4)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,((net.num_generators-net.get_num_slack_gens() +
                                            net.get_num_reg_gens())*3,net.num_vars))

            # load all
            P = net.get_var_projection('load','all',1,4)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(net.num_loads*4*2,net.num_vars))

            # branch all
            P = net.get_var_projection('branch','all')
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,((net.get_num_tap_changers_v() +
                                            net.get_num_phase_shifters())*self.T,net.num_vars))

            # shunt all
            P = net.get_var_projection('shunt','all',3,3)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(net.get_num_switched_shunts(),net.num_vars))

            # vargen all
            P = net.get_var_projection('variable generator','all',-1,2)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(2*net.num_var_generators*3,net.num_vars))

            # battery all
            P = net.get_var_projection('battery','all',3,6)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(net.num_batteries*3*2,net.num_vars))

            # all all
            P = net.get_var_projection('all','all',2,2)
            self.assertTrue(np.all(P.data == 1.))
            self.assertTupleEqual(P.shape,(net.num_vars/self.T,net.num_vars))
            for t in range(self.T):
                for bus in net.buses:
                    if not bus.is_slack():
                        a = np.where(P.col == bus.index_v_mag[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                        a = np.where(P.col == bus.index_v_ang[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                for gen in net.generators:
                    if not gen.is_slack():
                        a = np.where(P.col == gen.index_P[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                    if gen.is_regulator():
                        a = np.where(P.col == gen.index_Q[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                for branch in net.branches:
                    if branch.is_tap_changer_v():
                        a = np.where(P.col == branch.index_ratio[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                    if branch.is_phase_shifter():
                        a = np.where(P.col == branch.index_phase[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                for load in net.loads:
                    a = np.where(P.col == load.index_P[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                    a = np.where(P.col == load.index_Q[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                for shunt in net.shunts:
                    if shunt.is_switched_v():
                        a = np.where(P.col == shunt.index_b[t])[0]
                        if t == 2:
                            self.assertTupleEqual(a.shape,(1,))
                        else:
                            self.assertTupleEqual(a.shape,(0,))
                for vargen in net.var_generators:
                    a = np.where(P.col == vargen.index_P[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                    a = np.where(P.col == vargen.index_Q[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                for bat in net.batteries:
                    a = np.where(P.col == bat.index_Pc[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                    a = np.where(P.col == bat.index_Pd[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))
                    a = np.where(P.col == bat.index_E[t])[0]
                    if t == 2:
                        self.assertTupleEqual(a.shape,(1,))
                    else:
                        self.assertTupleEqual(a.shape,(0,))

            # Recovery
            Projs = []
            for t in range(self.T):
                Projs.append(net.get_var_projection('all','all',t,t))
            x = net.get_var_values()
            self.assertTupleEqual(x.shape,(net.num_vars,))
            self.assertEqual(net.num_periods,self.T)
            y = sum([P.T*P*x for P in Projs])
            self.assertLess(np.linalg.norm(y-x),np.linalg.norm(x)*1e-10)

    def test_variable_limits(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            self.assertGreater(net.num_buses,0)
            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))

            # Add batteries
            gen_buses = net.get_generator_buses()
            net.add_batteries(gen_buses,20.,30.,0.7,0.8)
            self.assertGreater(net.num_batteries,0)
            self.assertEqual(net.num_batteries,len(gen_buses))

            # Set batteries
            for bat in net.batteries:
                self.assertEqual(bat.P,0.)
                if bat.index % 2 == 0:
                    bat.P = 1.*(bat.index+1)
                else:
                    bat.P = -1.*(bat.index+1)

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)
                load.P_max = 3.3*(load.index+1)
                load.Q = 3.5*load.index
            self.assertEqual(net.num_loads,net.get_num_P_adjust_loads())

            # vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude',
                           'voltage angle',
                           'voltage magnitude deviation',
                           'voltage magnitude violation'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          ['active power','reactive power'])
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'variable',
                          'any',
                          ['tap ratio','tap ratio deviation','phase shift'])
            net.set_flags('shunt',
                          'variable',
                          'any',
                          ['susceptance','susceptance deviation'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,
                             (6*net.num_buses +
                              2*net.num_generators +
                              2*net.num_var_generators +
                              4*net.num_branches +
                              3*net.num_shunts +
                              2*net.get_num_P_adjust_loads()+
                              3*net.num_batteries))

            # Add some interesting vargen values
            for vargen in net.var_generators:
                vargen.P = 2*vargen.index
                vargen.P_ava = 2.5*vargen.index
                vargen.P_max = 3*vargen.index
                vargen.P_min = 9*vargen.index
                vargen.Q = 4*vargen.index
                vargen.Q_min = 1*vargen.index
                vargen.Q_max = 5*vargen.index

            # current
            x = net.get_var_values()
            self.assertEqual(x.size,net.num_vars)
            for bus in net.buses:
                self.assertEqual(x[bus.index_v_mag],bus.v_mag)
                self.assertEqual(x[bus.index_v_ang],bus.v_ang)
                self.assertEqual(x[bus.index_y],np.maximum(bus.v_mag-bus.v_set,0))
                self.assertEqual(x[bus.index_z],np.maximum(bus.v_set-bus.v_mag,0))
                self.assertEqual(x[bus.index_vl],0.)
                self.assertEqual(x[bus.index_vh],0.)
            for br in net.branches:
                self.assertEqual(x[br.index_ratio],br.ratio)
                self.assertEqual(x[br.index_ratio_y],0.)
                self.assertEqual(x[br.index_ratio_z],0.)
                self.assertEqual(x[br.index_phase],br.phase)
            for gen in net.generators:
                self.assertEqual(x[gen.index_P],gen.P)
                self.assertEqual(x[gen.index_Q],gen.Q)
            for load in net.loads:
                self.assertEqual(x[load.index_P],load.P)
                self.assertEqual(x[load.index_Q],load.Q)
                self.assertEqual(x[load.index_Q],3.5*load.index)
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],vargen.index*2)
                self.assertEqual(x[vargen.index_P],vargen.P)
                self.assertEqual(x[vargen.index_Q],vargen.index*4)
                self.assertEqual(x[vargen.index_Q],vargen.Q)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b)
                self.assertEqual(x[shunt.index_y],0.)
                self.assertEqual(x[shunt.index_z],0.)
            for bat in net.batteries:
                if bat.index % 2 == 0:
                    self.assertGreater(bat.P,0)
                    self.assertEqual(x[bat.index_Pc],bat.P)
                    self.assertEqual(x[bat.index_Pc],1.*(bat.index+1))
                    self.assertEqual(x[bat.index_Pd],0.)
                else:
                    self.assertLess(bat.P,0.)
                    self.assertEqual(x[bat.index_Pc],0.)
                    self.assertEqual(x[bat.index_Pd],-bat.P)
                    self.assertEqual(x[bat.index_Pd],1.*(bat.index+1))
                self.assertEqual(x[bat.index_E],bat.E)

            # upper limits
            x = net.get_var_values('upper limits')
            self.assertEqual(x.size,net.num_vars)
            self.assertEqual(pf.BUS_INF_V_MAG,100.)
            self.assertEqual(pf.BRANCH_INF_RATIO,100.)
            self.assertEqual(pf.SHUNT_INF_SUSC,1000.)
            for bus in net.buses:
                self.assertEqual(x[bus.index_v_mag],bus.v_max)
                self.assertEqual(x[bus.index_v_ang],pf.BUS_INF_V_ANG)
                self.assertEqual(x[bus.index_y],pf.BUS_INF_V_MAG)
                self.assertEqual(x[bus.index_z],pf.BUS_INF_V_MAG)
                self.assertEqual(x[bus.index_vl],pf.BUS_INF_V_MAG)
                self.assertEqual(x[bus.index_vh],pf.BUS_INF_V_MAG)
            for br in net.branches:
                self.assertEqual(x[br.index_ratio],br.ratio_max)
                self.assertEqual(x[br.index_ratio_y],pf.BRANCH_INF_RATIO)
                self.assertEqual(x[br.index_ratio_z],pf.BRANCH_INF_RATIO)
                self.assertEqual(x[br.index_phase],br.phase_max)
            for gen in net.generators:
                self.assertEqual(x[gen.index_P],gen.P_max)
                self.assertEqual(x[gen.index_Q],gen.Q_max)
            for load in net.loads:
                self.assertEqual(x[load.index_P],load.P_max)
                self.assertEqual(x[load.index_P],3.3*(load.index+1))
                self.assertEqual(x[load.index_Q],pf.LOAD_INF_Q)
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],2.5*vargen.index)
                self.assertEqual(x[vargen.index_P],vargen.P_ava)
                self.assertEqual(x[vargen.index_Q],5*vargen.index)
                self.assertEqual(x[vargen.index_Q],vargen.Q_max)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_max)
                self.assertEqual(x[shunt.index_y],pf.SHUNT_INF_SUSC)
                self.assertEqual(x[shunt.index_z],pf.SHUNT_INF_SUSC)
            for bat in net.batteries:
                self.assertEqual(x[bat.index_Pc],bat.P_max)
                self.assertEqual(x[bat.index_Pd],-bat.P_min)
                self.assertEqual(x[bat.index_E],bat.E_max)

            # lower limits
            x = net.get_var_values('lower limits')
            self.assertEqual(x.size,net.num_vars)
            for bus in net.buses:
                self.assertEqual(x[bus.index_v_mag],bus.v_min)
                self.assertEqual(x[bus.index_v_ang],-pf.BUS_INF_V_ANG)
                self.assertEqual(x[bus.index_y],0.)
                self.assertEqual(x[bus.index_z],0.)
                self.assertEqual(x[bus.index_vl],0.)
                self.assertEqual(x[bus.index_vh],0.)
            for br in net.branches:
                self.assertEqual(x[br.index_ratio],br.ratio_min)
                self.assertEqual(x[br.index_ratio_y],0.)
                self.assertEqual(x[br.index_ratio_z],0.)
                self.assertEqual(x[br.index_phase],br.phase_min)
            for gen in net.generators:
                self.assertEqual(x[gen.index_P],gen.P_min)
                self.assertEqual(x[gen.index_Q],gen.Q_min)
            for load in net.loads:
                self.assertEqual(x[load.index_P],load.P_min)
                self.assertEqual(x[load.index_P],-2.4*(load.index+1))
                self.assertEqual(x[load.index_Q],-pf.LOAD_INF_Q)
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],9.*vargen.index)
                self.assertEqual(x[vargen.index_P],vargen.P_min)
                self.assertEqual(x[vargen.index_Q],1.*vargen.index)
                self.assertEqual(x[vargen.index_Q],vargen.Q_min)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_min)
                self.assertEqual(x[shunt.index_y],0.)
                self.assertEqual(x[shunt.index_z],0.)
            for bat in net.batteries:
                self.assertEqual(x[bat.index_Pc],0.)
                self.assertEqual(x[bat.index_Pd],0.)
                self.assertEqual(x[bat.index_E],0.)

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            self.assertGreater(net.num_buses,0)
            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_var_generators(load_buses,80.,50.,30.,5,0.05)
            self.assertGreater(net.num_var_generators,0)
            self.assertEqual(net.num_var_generators,len(load_buses))

            # Add batteries
            gen_buses = net.get_generator_buses()
            net.add_batteries(gen_buses,20.,40.,0.7,0.8)
            self.assertGreater(net.num_batteries,0)
            self.assertEqual(net.num_batteries,len(gen_buses))

            # Set batteries
            for bat in net.batteries:
                self.assertTrue(np.all(bat.P == 0.))
                x = np.random.rand(net.num_periods)
                for t in range(net.num_periods):
                    if t % 2 == 0:
                        x[t] *= -1.
                self.assertTrue(np.any(x < 0))
                self.assertTrue(np.any(x > 0))
                bat.P = x
                for t in range(net.num_periods):
                    self.assertEqual(x[t],bat.P[t])
                self.assertTrue(np.any(bat.P < 0))
                self.assertTrue(np.any(bat.P > 0))

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)
                load.P_max = 3.3*(load.index+1)
                load.Q = 3.5*load.index*np.array(range(net.num_periods))
                for t in range(net.num_periods):
                    self.assertEqual(load.Q[t],3.5*load.index*t)
            self.assertEqual(net.num_loads,net.get_num_P_adjust_loads())

            # vars
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude',
                           'voltage angle',
                           'voltage magnitude deviation',
                           'voltage magnitude violation'])
            net.set_flags('generator',
                          'variable',
                          'any',
                          ['active power',
                           'reactive power'])
            net.set_flags('load',
                          'variable',
                          'adjustable active power',
                          ['active power','reactive power'])
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          ['active power','reactive power'])
            net.set_flags('branch',
                          'variable',
                          'any',
                          ['tap ratio',
                           'tap ratio deviation',
                           'phase shift'])
            net.set_flags('shunt',
                          'variable',
                          'any',
                          ['susceptance','susceptance deviation'])
            net.set_flags('battery',
                          'variable',
                          'any',
                          ['charging power','energy level'])
            self.assertEqual(net.num_vars,
                             (6*net.num_buses +
                              2*net.num_generators +
                              2*net.num_var_generators +
                              4*net.num_branches +
                              3*net.num_shunts +
                              2*net.get_num_P_adjust_loads()+
                              3*net.num_batteries)*net.num_periods)

            # Add some interesting vargen values
            for vargen in net.var_generators:
                vargen.P = 2.*vargen.index*np.array(range(self.T))
                vargen.P_ava = 2.5*vargen.index*np.array(range(self.T))
                vargen.P_max = 3*vargen.index
                vargen.P_min = 9*vargen.index
                vargen.Q = np.ones(self.T)*4*vargen.index
                vargen.Q_min = 1*vargen.index
                vargen.Q_max = 5*vargen.index

            # current
            x = net.get_var_values()
            self.assertEqual(x.size,net.num_vars)
            for t in range(self.T):
                for bus in net.buses:
                    self.assertEqual(x[bus.index_v_mag[t]],bus.v_mag[t])
                    self.assertEqual(x[bus.index_v_ang[t]],bus.v_ang[t])
                    self.assertEqual(x[bus.index_y[t]],np.maximum(bus.v_mag[t]-bus.v_set[t],0))
                    self.assertEqual(x[bus.index_z[t]],np.maximum(bus.v_set[t]-bus.v_mag[t],0))
                    self.assertEqual(x[bus.index_vl[t]],0.)
                    self.assertEqual(x[bus.index_vh[t]],0.)
                for br in net.branches:
                    self.assertEqual(x[br.index_ratio[t]],br.ratio[t])
                    self.assertEqual(x[br.index_ratio_y[t]],0.)
                    self.assertEqual(x[br.index_ratio_z[t]],0.)
                    self.assertEqual(x[br.index_phase[t]],br.phase[t])
                for gen in net.generators:
                    self.assertEqual(x[gen.index_P[t]],gen.P[t])
                    self.assertEqual(x[gen.index_Q[t]],gen.Q[t])
                for load in net.loads:
                    self.assertEqual(x[load.index_P[t]],load.P[t])
                    self.assertEqual(x[load.index_Q[t]],load.Q[t])
                    self.assertEqual(x[load.index_Q[t]],load.index*3.5*t)
                for vargen in net.var_generators:
                    self.assertEqual(x[vargen.index_P[t]],vargen.index*2.*t)
                    self.assertEqual(x[vargen.index_P[t]],vargen.P[t])
                    self.assertEqual(x[vargen.index_Q[t]],vargen.index*4)
                    self.assertEqual(x[vargen.index_Q[t]],vargen.Q[t])
                for shunt in net.shunts:
                    self.assertEqual(x[shunt.index_b[t]],shunt.b[t])
                    self.assertEqual(x[shunt.index_y[t]],0.)
                    self.assertEqual(x[shunt.index_z[t]],0.)
                for bat in net.batteries:
                    if bat.P[t] >= 0:
                        self.assertEqual(x[bat.index_Pc[t]],bat.P[t])
                        self.assertEqual(x[bat.index_Pd[t]],0.)
                    else:
                        self.assertEqual(x[bat.index_Pc[t]],0.)
                        self.assertEqual(x[bat.index_Pd[t]],-bat.P[t])
                    self.assertEqual(x[bat.index_E[t]],bat.E[t])

            # upper limits
            x = net.get_var_values('upper limits')
            self.assertEqual(x.size,net.num_vars)
            self.assertEqual(pf.BUS_INF_V_MAG,100.)
            self.assertEqual(pf.BRANCH_INF_RATIO,100.)
            self.assertEqual(pf.SHUNT_INF_SUSC,1000.)
            for t in range(self.T):
                for bus in net.buses:
                    self.assertEqual(x[bus.index_v_mag[t]],bus.v_max)
                    self.assertEqual(x[bus.index_v_ang[t]],pf.BUS_INF_V_ANG)
                    self.assertEqual(x[bus.index_y[t]],pf.BUS_INF_V_MAG)
                    self.assertEqual(x[bus.index_z[t]],pf.BUS_INF_V_MAG)
                    self.assertEqual(x[bus.index_vl[t]],pf.BUS_INF_V_MAG)
                    self.assertEqual(x[bus.index_vh[t]],pf.BUS_INF_V_MAG)
                for br in net.branches:
                    self.assertEqual(x[br.index_ratio[t]],br.ratio_max)
                    self.assertEqual(x[br.index_ratio_y[t]],pf.BRANCH_INF_RATIO)
                    self.assertEqual(x[br.index_ratio_z[t]],pf.BRANCH_INF_RATIO)
                    self.assertEqual(x[br.index_phase[t]],br.phase_max)
                for gen in net.generators:
                    self.assertEqual(x[gen.index_P[t]],gen.P_max)
                    self.assertEqual(x[gen.index_Q[t]],gen.Q_max)
                for load in net.loads:
                    self.assertEqual(x[load.index_P[t]],load.P_max)
                    self.assertEqual(x[load.index_P[t]],3.3*(load.index+1))
                    self.assertEqual(x[load.index_Q[t]],pf.LOAD_INF_Q)
                for vargen in net.var_generators:
                    self.assertEqual(x[vargen.index_P[t]],2.5*vargen.index*t)
                    self.assertEqual(x[vargen.index_P[t]],vargen.P_ava[t])
                    self.assertEqual(x[vargen.index_Q[t]],5*vargen.index)
                    self.assertEqual(x[vargen.index_Q[t]],vargen.Q_max)
                for shunt in net.shunts:
                    self.assertEqual(x[shunt.index_b[t]],shunt.b_max)
                    self.assertEqual(x[shunt.index_y[t]],pf.SHUNT_INF_SUSC)
                    self.assertEqual(x[shunt.index_z[t]],pf.SHUNT_INF_SUSC)
                for bat in net.batteries:
                    self.assertEqual(x[bat.index_Pc[t]],bat.P_max)
                    self.assertEqual(x[bat.index_Pd[t]],-bat.P_min)
                    self.assertEqual(x[bat.index_E[t]],bat.E_max)

            # lower limits
            x = net.get_var_values('lower limits')
            self.assertEqual(x.size,net.num_vars)
            for t in range(self.T):
                for bus in net.buses:
                    self.assertEqual(x[bus.index_v_mag[t]],bus.v_min)
                    self.assertEqual(x[bus.index_v_ang[t]],-pf.BUS_INF_V_ANG)
                    self.assertEqual(x[bus.index_y[t]],0.)
                    self.assertEqual(x[bus.index_z[t]],0.)
                    self.assertEqual(x[bus.index_vl[t]],0.)
                    self.assertEqual(x[bus.index_vh[t]],0.)
                for br in net.branches:
                    self.assertEqual(x[br.index_ratio[t]],br.ratio_min)
                    self.assertEqual(x[br.index_ratio_y[t]],0.)
                    self.assertEqual(x[br.index_ratio_z[t]],0.)
                    self.assertEqual(x[br.index_phase[t]],br.phase_min)
                for gen in net.generators:
                    self.assertEqual(x[gen.index_P[t]],gen.P_min)
                    self.assertEqual(x[gen.index_Q[t]],gen.Q_min)
                for load in net.loads:
                    self.assertEqual(x[load.index_P[t]],load.P_min)
                    self.assertEqual(x[load.index_P[t]],-2.4*(load.index+1))
                    self.assertEqual(x[load.index_Q[t]],-pf.LOAD_INF_Q)
                for vargen in net.var_generators:
                    self.assertEqual(x[vargen.index_P[t]],9.*vargen.index)
                    self.assertEqual(x[vargen.index_P[t]],vargen.P_min)
                    self.assertEqual(x[vargen.index_Q[t]],1.*vargen.index)
                    self.assertEqual(x[vargen.index_Q[t]],vargen.Q_min)
                for shunt in net.shunts:
                    self.assertEqual(x[shunt.index_b[t]],shunt.b_min)
                    self.assertEqual(x[shunt.index_y[t]],0.)
                    self.assertEqual(x[shunt.index_z[t]],0.)
                for bat in net.batteries:
                    self.assertEqual(x[bat.index_Pc[t]],0.)
                    self.assertEqual(x[bat.index_Pd[t]],0.)
                    self.assertEqual(x[bat.index_E[t]],0.)
                    
    def test_var_generators_P_sigma(self):

        # Single period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)

            lmin = np.min([l.P for l in net.loads])
            for l in net.loads:
                l.P = l.P + np.abs(lmin)

            spread = 2
            corr = 0.1
            P_MIN = 1e-5

            # Add renewable sources
            gen_buses = net.get_generator_buses()
            total_load = abs(sum([l.P for l in net.loads]))
            net.add_var_generators(gen_buses,80.,50.,30.,5,0.05)
            self.assertEqual(net.num_var_generators,len([b for b in net.buses if b.generators]))
            for vg in net.var_generators:
                self.assertTrue(isinstance(vg.name,str) or isinstance(vg.name,unicode))
                self.assertEqual(vg.name,"VARGEN %d" %(vg.index+1))
                self.assertEqual(net.get_var_generator_by_name(vg.name).name,vg.name)
                self.assertEqual(vg.P,0.5*vg.P_max)
                self.assertEqual(vg.P_ava,0.5*vg.P_max)
                self.assertEqual(vg.P_min,0)
                self.assertEqual(vg.P_max,0.8*total_load/net.num_var_generators)
                self.assertEqual(vg.P_std,0.3*vg.P_max)
                self.assertGreater(len(vg.bus.generators),0)
                self.assertNotEqual(vg.P_max,0)
                self.assertNotEqual(vg.P_std,0)
            self.assertLess(np.abs(sum([vg.P_max for vg in net.var_generators])-0.8*total_load),1e-10)

            # Variables
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            self.assertEqual(net.num_vars,net.num_var_generators)

            # Correlation
            sigma = net.create_var_generators_P_sigma(spread,corr)

            # Check
            self.assertTrue(np.all(sigma.row >= sigma.col))
            indexP2vargen = {}
            for vg in net.var_generators:
                indexP2vargen[vg.index_P] = vg
            for k in range(sigma.nnz):
                i = sigma.row[k]
                j = sigma.col[k]
                d = sigma.data[k]
                vg1 = indexP2vargen[i]
                vg2 = indexP2vargen[j]
                if i == j:
                    self.assertLess(np.abs(d - vg1.P_std**2.),1e-12)
                else:
                    self.assertLess(np.abs(d - corr*vg1.P_std*vg2.P_std),1e-12)

            # Posdef
            try:
                from scikits.sparse.cholmod import cholesky
                sigma = (sigma + sigma.T - triu(sigma)).tocoo()
                factor = cholesky(sigma.tocsc())
            except ImportError:
                pass

        # Multi period
        for case in test_cases.CASES:

            net = pf.Parser(case).parse(case,self.T)
            self.assertEqual(net.num_periods,self.T)

            lmin = np.min([l.P[0] for l in net.loads])
            for l in net.loads:
                l.P = l.P + np.abs(lmin)

            spread = 2
            corr = 0.1
            P_MIN = 1e-5

            # Add renewable sources
            gen_buses = net.get_generator_buses()
            total_load = abs(sum([l.P[0] for l in net.loads]))
            net.add_var_generators(gen_buses,80.,50.,30.,5,0.05)
            self.assertEqual(net.num_var_generators,len([b for b in net.buses if b.generators]))
            for vg in net.var_generators:
                self.assertEqual(net.get_var_generator_by_name(vg.name).name,vg.name)
                self.assertEqual(vg.P_min,0)
                self.assertEqual(vg.P_max,0.8*total_load/net.num_var_generators)
                self.assertGreater(len(vg.bus.generators),0)
                self.assertNotEqual(vg.P_max,0)
                for t in range(self.T):
                    self.assertEqual(vg.P[t],0.5*vg.P_max)
                    self.assertEqual(vg.P_ava[t],0.5*vg.P_max)
                    self.assertEqual(vg.P_std[t],0.3*vg.P_max)
                    self.assertNotEqual(vg.P_std[t],0)
            self.assertLess(np.abs(sum([vg.P_max for vg in net.var_generators])-0.8*total_load),1e-10)

            # Variables
            net.set_flags('variable generator',
                          'variable',
                          'any',
                          'active power')
            self.assertEqual(net.num_vars,net.num_var_generators*self.T)

            # Correlation
            sigma = net.create_var_generators_P_sigma(spread,corr)

            # Check
            self.assertTrue(np.all(sigma.row >= sigma.col))
            indexP2vargen = {}
            for vg in net.var_generators:
                for t in range(self.T):
                    indexP2vargen[vg.index_P[t]] = vg
            for k in range(sigma.nnz):
                i = sigma.row[k]
                j = sigma.col[k]
                d = sigma.data[k]
                vg1 = indexP2vargen[i]
                vg2 = indexP2vargen[j]
                for t in range(self.T):
                    if i == j:
                        self.assertLess(np.abs(d - vg1.P_std[t]**2.),1e-12)
                    else:
                        self.assertLess(np.abs(d - corr*vg1.P_std[t]*vg2.P_std[t]),1e-12)

    def tearDown(self):

        pass

def compute_branch_flows(parameters):
    """
    Compute branch flows for the given

    Parameters
    ----------
    parameters : dict

    Returns
    -------
    flows : dict
    """
    
    # Transformer tap ratios
    a_km = parameters['ratio']
    a_mk = 1.
    
    # Transformer phase shift
    phi = parameters['phase']
    
    # Voltage magnitude and angles
    v_k = parameters['bus_k.v_mag']
    w_k = parameters['bus_k.v_ang']
    v_m = parameters['bus_m.v_mag']
    w_m = parameters['bus_m.v_ang']
    
    # Conductances
    g_km = parameters['g']
    g_k_sh = parameters['g_k']
    g_mk = parameters['g']
    g_m_sh = parameters['g_m']
    
    # Susceptances
    b_km = parameters['b']
    b_k_sh = parameters['b_k']
    b_mk = parameters['b']
    b_m_sh = parameters['b_m']

    # Intermediate calculations
    v_k_tap_squared = math.pow(v_k,2)*math.pow(a_km,2)
    v_m_tap_squared = math.pow(v_m,2)*math.pow(a_mk,2)
    v_k_v_m_tap = v_k*v_m*a_km*a_mk
    cos_km = math.cos(w_k-w_m-phi)
    sin_km = math.sin(w_k-w_m-phi)
    cos_mk = math.cos(w_m-w_k+phi)
    sin_mk = math.sin(w_m-w_k+phi)

    flows = {}
    
    # Flows in shunt elements of pi model
    
    # P_k_sh = v_k^2*a_km^2*g_k_sh
    flows['P_k_sh'] = v_k_tap_squared*g_k_sh
    
    # Q_k_sh = -v_k^2*a_km^2*b_k_sh
    flows['Q_k_sh'] = -v_k_tap_squared*b_k_sh
    
    # P_m_sh = v_m^2*a_mk^2*g_m_sh
    flows['P_m_sh'] = v_m_tap_squared*g_m_sh
    
    # Q_m_sh = -v_m^2*a_mk^2*b_m_sh
    flows['Q_m_sh'] = -v_m_tap_squared*b_m_sh

    # Flows in series elements of pi model
    flows['P_km_ser'] = v_k_tap_squared*g_km - v_k_v_m_tap*(g_km*cos_km + b_km*sin_km)
    flows['Q_km_ser'] = -v_k_tap_squared*b_km - v_k_v_m_tap*(g_km*sin_km - b_km*cos_km)
    flows['P_mk_ser'] = v_m_tap_squared*g_mk - v_k_v_m_tap*(g_mk*cos_mk + b_mk*sin_mk)
    flows['Q_mk_ser'] = -v_m_tap_squared*b_mk - v_k_v_m_tap*(g_mk*sin_mk - b_mk*cos_mk)

    # Flows as measured from the bus
    flows['P_km'] = flows['P_km_ser'] + flows['P_k_sh']
    flows['Q_km'] = flows['Q_km_ser'] + flows['Q_k_sh']
    flows['P_mk'] = flows['P_mk_ser'] + flows['P_m_sh']
    flows['Q_mk'] = flows['Q_mk_ser'] + flows['Q_m_sh']

    return flows
