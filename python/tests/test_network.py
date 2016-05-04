#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
from . import test_cases
import numpy as np
from scipy.sparse import coo_matrix, bmat, triu

class TestNetwork(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_network(self):

        net = self.net

        for case in test_cases.CASES:

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

            net.load(case)
            
            self.assertGreater(net.num_buses,0)
            self.assertGreater(net.num_gens,0)
            self.assertGreater(net.num_branches,0)
            self.assertGreater(net.num_loads,0)
            self.assertGreaterEqual(net.num_shunts,0)
            self.assertGreaterEqual(net.num_vargens,0)
            self.assertGreaterEqual(net.num_bats,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            self.assertRaises(pf.NetworkError,net.get_bus,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_bus,net.num_buses)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_gen,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_gen,net.num_gens)
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
            self.assertRaises(pf.NetworkError,net.get_vargen,-1)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_vargen,net.num_vargens)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_bat,net.num_bats)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.get_bat,-1)
            net.clear_error()

            # Counters
            self.assertEqual(net.get_num_P_adjust_gens(),
                             len([g for g in net.generators if g.P_min < g.P_max]))
            self.assertEqual(net.get_num_P_adjust_loads(),
                             len([l for l in net.loads if l.P_min < l.P_max]))

    def test_variables(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            num_so_far = 0

            # P adjust gens
            num_adj = 0
            for gen in net.generators:
                if gen.P_min < gen.P_max:
                    num_adj += 1

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.get_num_P_adjust_gens(),num_adj)

            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_P_ADJUST,
                          pf.GEN_VAR_P)
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

            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_P_ADJUST,
                          pf.LOAD_VAR_P)
            num_so_far += net.get_num_P_adjust_loads()
            
            self.assertEqual(net.num_vars,
                             num_so_far)
           
            for load in net.loads:
                if load.index % 2 == 0:
                    self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                else:
                    self.assertFalse(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))

            # batter charging
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P)
            num_so_far += net.num_bats

            self.assertEqual(net.num_vars,
                             num_so_far)

            for bat in net.batteries:
                self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P))
                self.assertFalse(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_E))

            # batter energy
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_E)
            num_so_far += net.num_bats

            self.assertEqual(net.num_vars,
                             num_so_far)

            for bat in net.batteries:
                self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P))
                self.assertTrue(bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_E))

    def test_buses(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
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
                self.assertEqual(bus.obj_type,pf.OBJ_BUS)
                self.assertNotEqual(bus.obj_type,pf.OBJ_UNKNOWN)

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
                self.assertTrue(isinstance(bus.reg_gens,list))
                self.assertTrue(isinstance(bus.reg_trans,list))
                self.assertTrue(isinstance(bus.reg_shunts,list))
                self.assertTrue(isinstance(bus.gens,list))

                # name
                self.assertTrue(isinstance(bus.name,str) or isinstance(bus.name,unicode))
                bus.name = "some bus"
                self.assertTrue(isinstance(bus.name,str) or isinstance(bus.name,unicode))
                self.assertEqual(bus.name,"some bus")
                self.assertEqual(bus.name,'some bus')

                # generators
                for g in bus.gens:
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
                    self.assertTrue(bus == b.bus_from or bus == b.bus_to)

                # total injections
                self.assertEqual(bus.get_total_gen_P(),sum([g.P for g in bus.gens],0))
                self.assertEqual(bus.get_total_gen_Q(),sum([g.Q for g in bus.gens],0))
                self.assertEqual(bus.get_total_load_P(),sum([l.P for l in bus.loads],0))
                self.assertEqual(bus.get_total_load_Q(),sum([l.Q for l in bus.loads],0))

                # regulated by gen
                if bus.is_regulated_by_gen():
                    self.assertGreater(len(bus.reg_gens),0)
                    for gen in bus.reg_gens:
                        self.assertTrue(gen.is_regulator())
                        self.assertEqual(gen.reg_bus.number,bus.number)
                        
                else:
                    self.assertTrue(len(bus.reg_gens) == 0)

                # slack
                if bus.is_slack():
                    self.assertTrue(bus.is_regulated_by_gen())
                    self.assertGreater(len(bus.gens),0)
                    self.assertGreater(len(bus.reg_gens),0)
                    self.assertEqual(len(bus.gens),len(bus.reg_gens))
                    self.assertEqual(set([g.index for g in bus.gens]),set([g.index for g in bus.reg_gens]))
                    for gen in bus.gens:
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
                        if bus.number == tran.bus_from.number: # reg bus in from side -> neg sensitivity
                            self.assertFalse(tran.has_pos_ratio_v_sens()) 
                        elif bus.number == tran.bus_to.number: # reg bus in to side -> pos sensitivity
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
                self.assertTrue(isinstance(bus.branches_from,list))
                for br in bus.branches_from:
                    self.assertEqual(bus.number,br.bus_from.number)
                self.assertTrue(isinstance(bus.branches_to,list))
                for br in bus.branches_to:
                    self.assertEqual(bus.number,br.bus_to.number)
                self.assertGreater(len(bus.branches),0)
                self.assertEqual(len(bus.branches),len(bus.branches_from)+len(bus.branches_to))
                self.assertEqual(len(bus.branches),bus.degree)

                # vargens
                for vg in bus.vargens:
                    self.assertEqual(bus,vg.bus)

                # batters
                for b in bus.bats:
                    self.assertTrue(isinstance(b,pf.Battery))
                    self.assertEqual(bus,b.bus)

            # sum of degrees
            sum_deg = 0
            for i in range(net.num_buses):
                sum_deg += net.get_bus(i).degree
            self.assertEqual(sum_deg,2*net.num_branches)                
                            
    def test_gens(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertTrue(net.num_gens > 0)

            self.assertEqual(len(net.generators),net.num_gens)

            for i in range(net.num_gens):

                gen = net.get_gen(i)
                same_gen = net.get_gen(i)

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
                    j = net.num_gens-1
                if i != j:
                    other_gen = net.get_gen(j)
                    self.assertFalse(gen is other_gen)
                    self.assertFalse(gen == other_gen)
                    self.assertTrue(gen != other_gen)

                # obj type
                self.assertEqual(gen.obj_type,pf.OBJ_GEN)
                self.assertNotEqual(gen.obj_type,pf.OBJ_UNKNOWN)
 
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

    def test_branches(self):
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
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
                self.assertEqual(branch.obj_type,pf.OBJ_BRANCH)
                self.assertNotEqual(branch.obj_type,pf.OBJ_UNKNOWN)

                self.assertEqual(branch.index,net.branches[i].index)

                self.assertEqual(branch.index,i)
                
                self.assertTrue(branch.bus_from)
                self.assertTrue(branch.bus_to)
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

    def test_shunts(self):
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertGreaterEqual(net.num_shunts,0)

            self.assertEqual(len(net.shunts),net.num_shunts)

            for i in range(net.num_shunts):

                shunt = net.get_shunt(i)

                # obj type
                self.assertEqual(shunt.obj_type,pf.OBJ_SHUNT)
                self.assertNotEqual(shunt.obj_type,pf.OBJ_UNKNOWN)

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

    def test_loads(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertTrue(net.num_loads > 0)

            self.assertEqual(net.num_loads,len(net.loads))
            
            self.assertEqual(net.num_loads,sum([len(b.loads) for b in net.buses]))

            for i in range(net.num_loads):

                load = net.get_load(i)

                # obj type
                self.assertEqual(load.obj_type,pf.OBJ_LOAD)
                self.assertNotEqual(load.obj_type,pf.OBJ_UNKNOWN)

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

    def test_vargens(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(len(net.var_generators),net.num_vargens)

            # existing vargens
            for vg in net.var_generators:
                self.assertTrue(isinstance(vg.name,str) or isinstance(vg.name,unicode))
                self.assertEqual(net.get_vargen_by_name(vg.name).name,vg.name)
                vg.name = "some vargen"
                self.assertEqual(vg.name,"some vargen")
                self.assertRaises(pf.NetworkError,net.get_vargen_by_name,vg.name)

            # add vargens
            load_buses = net.get_load_buses()
            self.assertEqual(len(load_buses),
                             len([b for b in net.buses if b.loads]))
            total_load = sum([l.P for l in net.loads])
            net.add_vargens(load_buses,50.,30.,5,0.05)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len(load_buses))
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.index,i)
                self.assertEqual(vargen.obj_type,pf.OBJ_VARGEN)
                self.assertNotEqual(vargen.obj_type,pf.OBJ_UNKNOWN)
                self.assertTrue(isinstance(vargen.name,str) or isinstance(vargen.name,unicode))
                self.assertEqual(vargen.name,"VARGEN %d" %(vargen.index+1))
                self.assertEqual(net.get_vargen_by_name(vargen.name).name,vargen.name)
                self.assertEqual(vargen.P,0.5*vargen.P_max)
                self.assertEqual(vargen.P_max,total_load/net.num_vargens)
                self.assertEqual(vargen.P_min,0.)
                self.assertEqual(vargen.P_std,0.3*vargen.P_max)
                self.assertEqual(vargen.index,i)
                self.assertEqual(vargen.Q,0.)
                self.assertEqual(vargen.Q_max,0.)
                self.assertEqual(vargen.Q_min,0.)
                self.assertEqual(len(vargen.bus.vargens),1)
                self.assertEqual(vargen.bus.vargens[0].index,vargen.index)
                vargen.P = np.pi
                self.assertEqual(vargen.P,np.pi)
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
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.index,net.var_generators[i].index)
                self.assertTrue(vargen.bus)
                self.assertEqual(vargen.bus,load_buses[i])
                self.assertTrue(vargen.index in [vg.index for vg in load_buses[i].vargens])

            # Set P,P_max,P_min,P_std,Q,Q_max,Q_min
            self.assertGreater(net.num_vargens,0)
            self.assertGreater(len(net.var_generators),0)
            for vg in net.var_generators:
                self.assertEqual(vg.P,np.pi)
                self.assertEqual(vg.P_max,total_load/net.num_vargens)
                self.assertEqual(vg.P_std,0.3*vg.P_max)
                self.assertEqual(vg.P_min,0.)
                self.assertEqual(vg.Q,0.)
                self.assertEqual(vg.Q_max,0.)
                self.assertEqual(vg.Q_min,0.)
                vg.P = 1.
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
                self.assertEqual(vg.Q,4.)
                self.assertEqual(vg.Q_max,5.)
                self.assertEqual(vg.Q_min,6.)

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(len(net.var_generators),net.num_vargens)

            self.assertEqual(net.vargen_corr_radius,1.)
            self.assertEqual(net.vargen_corr_value,0.)

            # gen buses
            gen_buses = net.get_gen_buses()
            self.assertGreater(len(gen_buses),0)
            for b in gen_buses:
                self.assertTrue(b.gens is not None)
            self.assertEqual(len(gen_buses),
                             len([b for b in net.buses if b.gens]))

            # add vargens
            penetration = 50.
            uncertainty = 30.
            corr_radius = 5
            corr_value = 0.05
            self.assertRaises(pf.NetworkError,net.add_vargens,gen_buses,-10,50.,5,0.05)
            self.assertTrue(net.has_error())
            net.clear_error()
            self.assertFalse(net.has_error())
            self.assertRaises(pf.NetworkError,net.add_vargens,gen_buses,50,-10.,5,0.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_vargens,gen_buses,50,50.,-1,0.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_vargens,gen_buses,50,50.,5,1.05)
            net.clear_error()
            self.assertRaises(pf.NetworkError,net.add_vargens,gen_buses,50,50.,5,-1.05)
            net.clear_error()
            net.add_vargens([],penetration,uncertainty,corr_radius,corr_value)
            self.assertEqual(net.num_vargens,0)
            self.assertEqual(net.var_generators,[])

            net.add_vargens(gen_buses,penetration,uncertainty,corr_radius,corr_value)
            self.assertEqual(net.num_vargens,len(gen_buses))
            self.assertEqual(len(net.var_generators),len(gen_buses))

            self.assertEqual(net.vargen_corr_radius,corr_radius)
            self.assertEqual(net.vargen_corr_value,corr_value)
            
            total_load = sum([l.P for l in net.loads])            
            total_cap = sum([vg.P_max for vg in net.var_generators])

            self.assertLess(np.abs(total_load-total_cap),1e-10)
            self.assertLess(np.abs(penetration-100*sum([vg.P for vg in net.var_generators])/total_load),1e-10)
                        
            for vg in net.var_generators:
                self.assertEqual(vg.P_min,0.)
                self.assertEqual(vg.P_max,total_load/net.num_vargens)
                self.assertEqual(vg.P,(penetration/100.)*vg.P_max)
                self.assertEqual(vg.P_std,(uncertainty/100.)*vg.P_max)

    def test_batteries(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertGreaterEqual(net.num_bats,0)

            self.assertEqual(net.num_bats,len(net.batteries))
            
            self.assertEqual(net.num_bats,sum([len(b.bats) for b in net.buses]))

            for i in range(net.num_bats):

                bat = net.get_bat(i)

                self.assertTrue(isinstance(bat,pf.Battery))

                # obj type
                self.assertEqual(bat.obj_type,pf.OBJ_BAT)
                self.assertNotEqual(bat.obj_type,pf.OBJ_UNKNOWN)

                self.assertEqual(bat.index,i)

                self.assertEqual(bat.index,net.batteries[i].index)

                self.assertTrue(bat.bus)

                self.assertTrue(bat.index in map(lambda x: x.index,bat.bus.bats))

                # P
                bat.P_min = -1.23
                bat.P_max = 2123.
                bat.P = 0.3241
                self.assertEqual(bat.P_min,-1.23)
                self.assertEqual(bat.P_max,2123.)
                self.assertEqual(bat.P,0.3241)

                # E
                bat.E_max = 22.2
                bat.E = 0.555
                self.assertEqual(bat.E_max,22.2)
                self.assertEqual(bat.E,0.555)

                # eta
                bat.eta_c = 0.91
                bat.eta_d = 0.95
                self.assertEqual(bat.eta_c,0.91)
                self.assertEqual(bat.eta_d,0.95)
                                
    def test_clear_flags(self):
       
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_E)

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_FIXED,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_FIXED,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_FIXED,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_FIXED,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_E)

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))
            self.assertEqual(net.num_fixed,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))
            self.assertEqual(net.num_bounded,0)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_BOUNDED,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_BOUNDED,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_BOUNDED,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_BOUNDED,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_E)

            self.assertGreater(net.num_vars,0)
            self.assertGreater(net.num_fixed,0)
            self.assertGreater(net.num_bounded,0)

            self.assertEqual(net.num_vars,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))
            self.assertEqual(net.num_fixed,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))
            self.assertEqual(net.num_bounded,
                             (net.num_buses+
                              net.num_gens+
                              net.num_loads+
                              net.num_branches+
                              net.num_shunts+
                              net.num_bats))

            net.clear_flags()

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

    def test_properties(self):

        net = self.net
        
        for case in test_cases.CASES:

            net.clear_properties()

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
            
            net.load(case)

            # add vargens
            net.add_vargens(net.get_load_buses(),50.,30.,5,0.05)
            for vargen in net.var_generators:
                self.assertTrue(isinstance(vargen.name,str) or isinstance(vargen.name,unicode))
                self.assertEqual(vargen.name,"VARGEN %d" %(vargen.index+1))
                self.assertEqual(net.get_vargen_by_name(vargen.name).name,vargen.name)
                vargen.P = 1.
                vargen.Q = 2.
                self.assertGreater(len(vargen.bus.loads),0)

            net.update_properties()

            self.assertEqual(net.num_vars,0)

            bus_v_max = net.bus_v_max
            bus_P_mis = net.bus_P_mis

            net.clear_properties()
            net.update_properties()

            self.assertEqual(bus_v_max,net.bus_v_max)
            self.assertEqual(bus_P_mis,net.bus_P_mis)
            
            x0 = net.get_var_values()

            self.assertEqual(x0.size,0)

            net.update_properties(x0)

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
                
            # Mismatches
            constr = pf.Constraint(pf.CONSTR_TYPE_PF,net)
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

            net.clear_properties()
            
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

    def test_bus_sorting(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_bounded,0)

            # Variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)

            # Bounds
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)

            x0 = net.get_var_values()
            
            self.assertEqual(net.num_vars,2*net.num_buses)
            self.assertEqual(net.num_bounded,net.num_buses)
            
            # Store random bus sensitivities
            constr = [pf.Constraint(pf.CONSTR_TYPE_PF,net),
                      pf.Constraint(pf.CONSTR_TYPE_BOUND,net),
                      pf.Constraint(pf.CONSTR_TYPE_REG_GEN,net),
                      pf.Constraint(pf.CONSTR_TYPE_REG_TRAN,net),
                      pf.Constraint(pf.CONSTR_TYPE_REG_SHUNT,net)]
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

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            net.update_set_points()
            
            for bus in net.buses:
                if bus.is_regulated_by_gen():
                    self.assertEqual(bus.v_mag,bus.v_set)

    def test_projections(self):
       
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertEqual(net.num_vars,0)

            # Add vargens
            load_buses = net.get_load_buses()
            net.add_vargens(load_buses,50.,30.,5,0.05)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len(load_buses))

            # bus vmag and vang
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            
            # gen powers
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_SLACK,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)

            # load active powers
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_ANY,
                          pf.LOAD_VAR_P)
            
            # branch ratio and phase
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER_V,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            
            # shunt
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)

            # vargens
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)

            # batteries
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)

            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-1) +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers_v() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              2*net.get_num_vargens()+
                              net.num_loads+
                              2*net.num_bats))

            # set vargens
            for vargen in net.var_generators:
                vargen.P = np.pi*vargen.index
                vargen.Q = -np.pi*vargen.index

            # var values (current values)
            x = net.get_var_values()

            # bus vmag
            P = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VMAG)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_buses-1)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_buses-1)
            vmag = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VMAG):
                    self.assertEqual(vmag[index],bus.v_mag)
                    index += 1

            # bus vang
            P = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_buses-1)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_buses-1)
            vang = P*x
            index = 0
            for i in range(net.num_buses):
                bus = net.get_bus(i)
                if bus.has_flags(pf.FLAG_VARS,pf.BUS_VAR_VANG):
                    self.assertEqual(vang[index],bus.v_ang)
                    index += 1

            # gen active power
            P = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_slack_gens())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_slack_gens())
            gP = P*x
            index = 0
            for i in range(net.num_gens):
                gen = net.get_gen(i)
                if gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P):
                    self.assertEqual(gP[index],gen.P)
                    index += 1

            # gen reactive power
            P = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_Q)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_reg_gens())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_reg_gens())
            gQ = P*x
            index = 0
            for i in range(net.num_gens):
                gen = net.get_gen(i)
                if gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_Q):
                    self.assertEqual(gQ[index],gen.Q)
                    index += 1

            # load active power
            P = net.get_var_projection(pf.OBJ_LOAD,pf.LOAD_VAR_P)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_loads)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_loads)
            gP = P*x
            index = 0
            for i in range(net.num_loads):
                load = net.get_load(i)
                self.assertTrue(load.has_flags(pf.FLAG_VARS,pf.LOAD_VAR_P))
                self.assertEqual(gP[index],load.P)
                index += 1
                    
            # tap changer ratio
            P = net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_RATIO)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_tap_changers_v())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_tap_changers_v())
            bR = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO):
                    self.assertEqual(bR[index],br.ratio)
                    index += 1
                    
            # phase shifter
            P = net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_PHASE)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_phase_shifters())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_phase_shifters())
            bP = P*x
            index = 0
            for i in range(net.num_branches):
                br = net.get_branch(i)
                if br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_PHASE):
                    self.assertEqual(bP[index],br.phase)
                    index += 1
                    
            # shunt susceptance
            P = net.get_var_projection(pf.OBJ_SHUNT,pf.SHUNT_VAR_SUSC)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.get_num_switched_shunts())
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.get_num_switched_shunts())
            sS = P*x
            index = 0
            for i in range(net.num_shunts):
                shunt = net.get_shunt(i)
                if shunt.has_flags(pf.FLAG_VARS,pf.SHUNT_VAR_SUSC):
                    self.assertEqual(sS[index],shunt.b)
                    index += 1
                
            # vargen active power
            P = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_vargens)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_vargens)
            vgP = P*x
            index = 0
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.index_P,vargen.index_Q-1)
                if vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P):
                    self.assertEqual(vgP[index],vargen.P)
                    self.assertEqual(vgP[index],vargen.index*np.pi)
                    index += 1

            # vargen reactive power
            P = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_Q)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_vargens)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_vargens)
            vgQ = P*x
            index = 0
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.index_P+1,vargen.index_Q)
                if vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_Q):
                    self.assertEqual(vgQ[index],vargen.Q)
                    self.assertEqual(vgQ[index],-vargen.index*np.pi)
                    index += 1

            # battery charging
            P = net.get_var_projection(pf.OBJ_BAT,pf.BAT_VAR_P)
            self.assertTrue(isinstance(P,coo_matrix))
            self.assertEqual(P.shape[0],net.num_bats)
            self.assertEqual(P.shape[1],net.num_vars)
            self.assertEqual(P.nnz,net.num_bats)
            batP = P*x
            index = 0
            for i in range(net.num_bats):
                bat = net.get_bat(i)
                self.assertEqual(bat.index_P,bat.index_E-1)
                if bat.has_flags(pf.FLAG_VARS,pf.BAT_VAR_P):
                    self.assertEqual(batP[index],bat.P)
                    index += 1
            print 'prin'

            # All
            Plist = [net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VMAG),
                     net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG),
                     net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P),
                     net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_Q),
                     net.get_var_projection(pf.OBJ_LOAD,pf.LOAD_VAR_P),           
                     net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_RATIO),
                     net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_PHASE),
                     net.get_var_projection(pf.OBJ_SHUNT,pf.SHUNT_VAR_SUSC),
                     net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P),
                     net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_Q),
                     net.get_var_projection(pf.OBJ_BAT,pf.BAT_VAR_P),
                     net.get_var_projection(pf.OBJ_BAT,pf.BAT_VAR_E)]
            P = bmat([[P] for P in Plist if P.shape[0] > 0])
            self.assertTupleEqual(P.shape,(net.num_vars,net.num_vars))
            for i in range(10):
                x = np.random.randn(net.num_vars)
                self.assertLess(np.linalg.norm(x-P.T*P*x),1e-12)

    def test_variable_limits(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertGreater(net.num_buses,0)            
            self.assertEqual(net.num_vars,0)
            
            # vargens
            load_buses = net.get_load_buses()
            net.add_vargens(load_buses,50.,30.,5,0.05)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len(load_buses))

            # loads
            for load in net.loads:
                load.P_min = -2.4*(load.index+1)
                load.P_max = 3.3*(load.index+1)
            self.assertEqual(net.num_loads,net.get_num_P_adjust_loads())

            # vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_LOAD,
                          pf.FLAG_VARS,
                          pf.LOAD_PROP_P_ADJUST,
                          pf.LOAD_VAR_P)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P|pf.VARGEN_VAR_Q)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV|pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV)
            net.set_flags(pf.OBJ_BAT,
                          pf.FLAG_VARS,
                          pf.BAT_PROP_ANY,
                          pf.BAT_VAR_P|pf.BAT_VAR_E)
            self.assertEqual(net.num_vars,
                             (6*net.num_buses +
                              2*net.num_gens + 
                              2*net.num_vargens +
                              4*net.num_branches + 
                              3*net.num_shunts +
                              net.get_num_P_adjust_loads()+
                              2*net.num_bats))
            
            # Add some interesting vargen values
            for vargen in net.var_generators:
                vargen.P = 2*vargen.index
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
                self.assertEqual(x[bat.index_P],bat.P)
                self.assertEqual(x[bat.index_E],bat.E)
                
            # upper limits
            x = net.get_var_values(pf.UPPER_LIMITS)
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
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],3*vargen.index)
                self.assertEqual(x[vargen.index_P],vargen.P_max)
                self.assertEqual(x[vargen.index_Q],5*vargen.index)
                self.assertEqual(x[vargen.index_Q],vargen.Q_max)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_max)
                self.assertEqual(x[shunt.index_y],pf.SHUNT_INF_SUSC)
                self.assertEqual(x[shunt.index_z],pf.SHUNT_INF_SUSC)
            for bat in net.batteries:
                self.assertEqual(x[bat.index_P],bat.P_max)
                self.assertEqual(x[bat.index_E],bat.E_max)

            # lower limits
            x = net.get_var_values(pf.LOWER_LIMITS)
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
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],9.*vargen.index)
                self.assertEqual(x[vargen.index_Q],1.*vargen.index)
                self.assertEqual(x[vargen.index_Q],vargen.Q_min)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_min)
                self.assertEqual(x[shunt.index_y],0.)
                self.assertEqual(x[shunt.index_z],0.)
            for bat in net.batteries:
                self.assertEqual(x[bat.index_P],bat.P_min)
                self.assertEqual(x[bat.index_E],0.)

    def test_vargen_P_sigma(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            spread = 2
            corr = 0.1
            P_MIN = 1e-5
        
            # Add renewable sources (at load buses)
            gen_buses = net.get_gen_buses()
            total_load = sum([l.P for l in net.loads])
            net.add_vargens(gen_buses,50.,30.,5,0.05)
            self.assertEqual(net.num_vargens,len([b for b in net.buses if b.gens]))
            for vg in net.var_generators:
                self.assertTrue(isinstance(vg.name,str) or isinstance(vg.name,unicode))
                self.assertEqual(vg.name,"VARGEN %d" %(vg.index+1))
                self.assertEqual(net.get_vargen_by_name(vg.name).name,vg.name)
                self.assertEqual(vg.P,0.5*vg.P_max)
                self.assertEqual(vg.P_min,0)
                self.assertEqual(vg.P_max,total_load/net.num_vargens)
                self.assertEqual(vg.P_std,0.3*vg.P_max)
                self.assertGreater(len(vg.bus.gens),0)
                self.assertNotEqual(vg.P_max,0)
                self.assertNotEqual(vg.P_std,0)
            self.assertLess(np.abs(sum([vg.P_max for vg in net.var_generators])-total_load),1e-10)

            # Variables
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            self.assertEqual(net.num_vars,net.num_vargens)
            
            # Correlation
            sigma = net.create_vargen_P_sigma(spread,corr)

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

    def tearDown(self):
        
        pass
                
                
