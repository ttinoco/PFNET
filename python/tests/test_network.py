#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
import test_cases
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
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)
            
            net.load(case)
            
            self.assertGreater(net.num_buses,0)
            self.assertGreater(net.num_gens,0)
            self.assertGreater(net.num_branches,0)
            self.assertGreater(net.num_shunts,0)
            self.assertEqual(net.num_vargens,0)
            self.assertEqual(net.num_vars,0)
            self.assertEqual(net.num_fixed,0)
            self.assertEqual(net.num_bounded,0)

            self.assertRaises(pf.NetworkError,net.get_bus,-1)
            self.assertRaises(pf.NetworkError,net.get_bus,net.num_buses)
            self.assertRaises(pf.NetworkError,net.get_gen,-1)
            self.assertRaises(pf.NetworkError,net.get_gen,net.num_gens)
            self.assertRaises(pf.NetworkError,net.get_branch,-1)
            self.assertRaises(pf.NetworkError,net.get_branch,net.num_branches)
            self.assertRaises(pf.NetworkError,net.get_shunt,-1)
            self.assertRaises(pf.NetworkError,net.get_shunt,net.num_shunts)
            self.assertRaises(pf.NetworkError,net.get_load,-1)
            self.assertRaises(pf.NetworkError,net.get_load,net.num_loads)
            self.assertRaises(pf.NetworkError,net.get_vargen,-1)
            self.assertRaises(pf.NetworkError,net.get_vargen,net.num_vargens)            

            # Counters
            self.assertEqual(net.get_num_P_adjust_gens(),
                             len([g for g in net.generators if g.P_min < g.P_max]))

    def test_variables(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            self.assertEqual(net.num_vars,0)
            self.assertGreater(net.get_num_P_adjust_gens(),0)

            # P adjust gens
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_P_ADJUST,
                          pf.GEN_VAR_P)
            self.assertEqual(net.num_vars,
                             net.get_num_P_adjust_gens())
            
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
                
                # hash
                self.assertEqual(bus.number,net.get_bus_by_number(bus.number).number)

                # values
                self.assertGreater(bus.number,0)
                self.assertGreater(bus.v_mag,0)
                self.assertGreater(bus.v_set,0)
                self.assertEqual(bus.sens_P_balance,0.)
                self.assertEqual(bus.sens_Q_balance,0.)
                self.assertEqual(bus.sens_v_mag_u_bound,0.)
                self.assertEqual(bus.sens_v_mag_l_bound,0.)
                self.assertEqual(bus.sens_v_reg_by_gen,0.)
                self.assertEqual(bus.sens_v_reg_by_tran,0.)
                self.assertEqual(bus.sens_v_reg_by_shunt,0.)
                self.assertTrue(isinstance(bus.reg_gens,list))
                self.assertTrue(isinstance(bus.reg_trans,list))
                self.assertTrue(isinstance(bus.reg_shunts,list))
                self.assertTrue(isinstance(bus.gens,list))

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
                self.assertEqual(bus.vargens,[])
                
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

                self.assertEqual(gen.index,i)

                self.assertEqual(gen.index,net.generators[i].index)

                self.assertTrue(gen.bus)
            
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

                # setting P_min and P_max
                self.assertNotEqual(gen.P_max,2*np.pi)
                self.assertNotEqual(gen.P_min,np.pi)
                gen.P_min = np.pi
                gen.P_max = 2*np.pi
                self.assertEqual(gen.P_min,np.pi)
                self.assertEqual(gen.P_max,2*np.pi)

    def test_branches(self):
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertTrue(net.num_branches > 0)

            self.assertEqual(net.num_branches,len(net.branches))

            for i in range(net.num_branches):

                branch = net.get_branch(i)

                self.assertEqual(branch.index,net.branches[i].index)

                self.assertEqual(branch.index,i)
                
                self.assertTrue(branch.bus_from)
                self.assertTrue(branch.bus_to)
                self.assertGreater(branch.ratio,0)

                # tap changer v
                if branch.is_tap_changer_v():
                    self.assertTrue(branch.is_tap_changer())
                    self.assertTrue(branch.reg_bus)
                    self.assertTrue(branch.index in [b.index for b in branch.reg_bus.reg_trans])
                    self.assertTrue(branch.reg_bus.is_regulated_by_tran())
                else:
                    self.assertRaises(pf.BusError,lambda : branch.reg_bus)

    def test_shunts(self):
        
        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            self.assertTrue(net.num_shunts > 0)

            self.assertEqual(len(net.shunts),net.num_shunts)

            for i in range(net.num_shunts):

                shunt = net.get_shunt(i)

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

                self.assertEqual(load.index,i)

                self.assertEqual(load.index,net.loads[i].index)

                self.assertTrue(load.bus)

                load.P = 0.3241
                load.Q = 0.1212
                self.assertEqual(load.P,0.3241)
                self.assertEqual(load.Q,0.1212)

    def test_vargens(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            self.assertEqual(net.num_vargens,0)
            self.assertEqual(len(net.var_generators),net.num_vargens)

            # Allocate
            load_buses = net.get_load_buses()
            self.assertEqual(len(load_buses),
                             len([b for b in net.buses if b.loads]))
            vargen_array = pf.VarGeneratorArray(len(load_buses))
            self.assertEqual(vargen_array.size,len(load_buses))

            # Set array
            net.set_vargen_array(vargen_array)
            self.assertGreater(net.num_vargens,0)
            self.assertEqual(net.num_vargens,len(load_buses))
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.P,0.)
                self.assertEqual(vargen.P_max,0.)
                self.assertEqual(vargen.P_std,0.)
                self.assertEqual(vargen.index,i)
                self.assertRaises(pf.BusError,lambda : vargen.bus)
                vargen.P = np.pi
                self.assertEqual(vargen.P,np.pi)
                vargen.P = 0.

            # Set buses
            net.set_vargen_buses(load_buses)
            for i in range(net.num_vargens):
                vargen = net.get_vargen(i)
                self.assertEqual(vargen.index,net.var_generators[i].index)
                self.assertTrue(vargen.bus)
                self.assertEqual(vargen.bus,load_buses[i])
                self.assertTrue(vargen.index in [vg.index for vg in load_buses[i].vargens])

            # Set P,P_max,P_std
            self.assertGreater(net.num_vargens,0)
            self.assertGreater(len(net.var_generators),0)
            for vg in net.var_generators:
                self.assertEqual(vg.P,0.)
                self.assertEqual(vg.P_max,0.)
                self.assertEqual(vg.P_std,0.)
                vg.P = 1.
                vg.P_max = 2.
                vg.P_std = 3.
                self.assertEqual(vg.P,1.)
                self.assertEqual(vg.P_max,2.)
                self.assertEqual(vg.P_std,3)                
                
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
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_FIXED,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_FIXED,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_FIXED,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_FIXED,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)

            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_BOUNDED,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_BOUNDED,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_BOUNDED,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_BOUNDED,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC)

            self.assertGreater(net.num_vars,0)
            self.assertGreater(net.num_fixed,0)
            self.assertGreater(net.num_bounded,0)

            self.assertEqual(net.num_vars,
                             net.num_buses+net.num_gens+net.num_branches+net.num_shunts)
            self.assertEqual(net.num_fixed,
                             net.num_buses+net.num_gens+net.num_branches+net.num_shunts)
            self.assertEqual(net.num_bounded,
                             net.num_buses+net.num_gens+net.num_branches+net.num_shunts)

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
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.gen_P_vio,0.)
            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.tran_p_vio,0.)
            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)
            self.assertEqual(net.num_actions,0)
            
            net.load(case)

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
            self.assertGreaterEqual(net.gen_v_dev,0.)
            self.assertGreaterEqual(net.gen_Q_vio,0.)
            self.assertGreaterEqual(net.gen_P_vio,0.)
            self.assertGreaterEqual(net.tran_v_vio,0.)
            self.assertGreaterEqual(net.tran_r_vio,0.)
            self.assertGreaterEqual(net.tran_p_vio,0.)
            self.assertGreaterEqual(net.shunt_v_vio,0.)
            self.assertGreaterEqual(net.shunt_b_vio,0.)
            self.assertGreaterEqual(net.num_actions,0.)

            self.assertEqual(net.bus_v_max,net.get_properties()['bus_v_max'])
            self.assertEqual(net.bus_v_min,net.get_properties()['bus_v_min'])
            self.assertEqual(net.bus_v_vio,net.get_properties()['bus_v_vio'])
            self.assertEqual(net.bus_P_mis,net.get_properties()['bus_P_mis'])
            self.assertEqual(net.bus_Q_mis,net.get_properties()['bus_Q_mis'])
            self.assertEqual(net.gen_v_dev,net.get_properties()['gen_v_dev'])
            self.assertEqual(net.gen_Q_vio,net.get_properties()['gen_Q_vio'])
            self.assertEqual(net.gen_P_vio,net.get_properties()['gen_P_vio'])
            self.assertEqual(net.tran_v_vio,net.get_properties()['tran_v_vio'])
            self.assertEqual(net.tran_r_vio,net.get_properties()['tran_r_vio'])
            self.assertEqual(net.tran_p_vio,net.get_properties()['tran_p_vio'])
            self.assertEqual(net.shunt_v_vio,net.get_properties()['shunt_v_vio'])
            self.assertEqual(net.shunt_b_vio,net.get_properties()['shunt_b_vio'])            
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
            for gen in net.generators:
                dP = np.max([gen.P-gen.P_max,gen.P_min-gen.P,0.])
                if dP > Pvio:
                    Pvio = dP
                if gen.is_regulator():
                    dQ = np.max([gen.Q-gen.Q_max,gen.Q_min-gen.Q,0.])
                    if dQ > Qvio:
                        Qvio = dQ                    
            self.assertLess(abs(net.gen_P_vio-Pvio*net.base_power),1e-10)
            self.assertLess(abs(net.gen_Q_vio-Qvio*net.base_power),1e-10)
            
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

            net.clear_properties()
            
            self.assertEqual(net.bus_v_max,0.)
            self.assertEqual(net.bus_v_min,0.)
            self.assertEqual(net.bus_v_vio,0.)
            self.assertEqual(net.bus_P_mis,0.)
            self.assertEqual(net.bus_Q_mis,0.)
            self.assertEqual(net.gen_v_dev,0.)
            self.assertEqual(net.gen_Q_vio,0.)
            self.assertEqual(net.gen_P_vio,0.)
            self.assertEqual(net.tran_v_vio,0.)
            self.assertEqual(net.tran_r_vio,0.)
            self.assertEqual(net.shunt_v_vio,0.)
            self.assertEqual(net.shunt_b_vio,0.)
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
                c.store_sensitivities(np.random.randn(c.f.size))
                
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
            r2 = []
            r3 = []
            r4 = []
            r5 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.get_largest_mis()) >= abs(bus2.get_largest_mis()))
                    r2.append(abs(bus2.get_largest_mis()) >= abs(bus1.get_largest_mis()))
                    r3.append(abs(bus1.get_largest_sens()) >= abs(bus2.get_largest_sens()))
                    r4.append(bus1.get_largest_mis() >= bus2.get_largest_mis())
                    r5.append(abs(bus1.P_mis) >= abs(bus2.P_mis))
            self.assertTrue(all(r1))
            self.assertFalse(all(r2))
            self.assertFalse(all(r3))
            self.assertFalse(all(r4))
            self.assertFalse(all(r5))

            # Sort by P mismatch
            bus_list = net.create_sorted_bus_list(pf.BUS_MIS_ACTIVE)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            r2 = []
            r3 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.P_mis) >= abs(bus2.P_mis))
                    r2.append(abs(bus2.P_mis) >= abs(bus1.P_mis))
                    r3.append(abs(bus1.get_largest_mis()) >= abs(bus2.get_largest_mis()))
            self.assertTrue(all(r1))
            self.assertFalse(all(r2))
            self.assertFalse(all(r3))

            # Sort by largest sensitivity
            bus_list = net.create_sorted_bus_list(pf.BUS_SENS_LARGEST)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            r2 = []
            r3 = []
            r4 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus1.get_largest_mis()) >= abs(bus2.get_largest_mis()))
                    r2.append(abs(bus2.get_largest_mis()) >= abs(bus1.get_largest_mis()))
                    r3.append(abs(bus1.get_largest_sens()) >= abs(bus2.get_largest_sens()))
                    r4.append(bus1.get_largest_sens() >= bus2.get_largest_sens())
            self.assertFalse(all(r1))
            self.assertFalse(all(r2))
            self.assertTrue(all(r3))
            self.assertFalse(all(r4))

            # Sort by v_reg_by_gen sensitivity
            bus_list = net.create_sorted_bus_list(pf.BUS_SENS_V_REG_BY_GEN)
            self.assertTrue(isinstance(bus_list,list))
            self.assertEqual(len(bus_list),net.num_buses)
            r1 = []
            r2 = []
            r3 = []
            for i in range(len(bus_list)):
                if i > 0:
                    bus1 = bus_list[i-1]
                    bus2 = bus_list[i]
                    self.assertTrue(isinstance(bus1,pf.Bus))
                    self.assertTrue(isinstance(bus2,pf.Bus))
                    r1.append(abs(bus2.sens_v_reg_by_gen) >= abs(bus1.sens_v_reg_by_gen))
                    r2.append(abs(bus1.sens_v_reg_by_gen) >= abs(bus2.sens_v_reg_by_gen))
                    r3.append(abs(bus1.get_largest_sens()) >= abs(bus2.get_largest_sens()))
            self.assertFalse(all(r1))
            self.assertTrue(all(r2))
            self.assertFalse(all(r3))

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
            self.assertEqual(net.num_vargens,0)

            # Add vargens
            load_buses = net.get_load_buses()
            vargen_array = pf.VarGeneratorArray(len(load_buses))
            net.set_vargen_array(vargen_array)
            net.set_vargen_buses(load_buses)
            self.assertGreater(net.num_vargens,0)

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
                          pf.VARGEN_VAR_P)

            self.assertEqual(net.num_vars,
                             (2*(net.num_buses-1) +
                              net.get_num_slack_gens() +
                              net.get_num_reg_gens() +
                              net.get_num_tap_changers_v() +
                              net.get_num_phase_shifters() +
                              net.get_num_switched_shunts() +
                              net.get_num_vargens()))

            # set vargens
            for vargen in net.var_generators:
                vargen.P = np.pi*vargen.index

            # var values
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
                if vargen.has_flags(pf.FLAG_VARS,pf.VARGEN_VAR_P):
                    self.assertEqual(vgP[index],vargen.P)
                    self.assertEqual(vgP[index],vargen.index*np.pi)
                    index += 1

            # All
            Plist = [net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VMAG),
                     net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG),
                     net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P),
                     net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_Q),           
                     net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_RATIO),
                     net.get_var_projection(pf.OBJ_BRANCH,pf.BRANCH_VAR_PHASE),
                     net.get_var_projection(pf.OBJ_SHUNT,pf.SHUNT_VAR_SUSC),
                     net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)]
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
            vargen_array = pf.VarGeneratorArray(len(load_buses))
            net.set_vargen_array(vargen_array)
            net.set_vargen_buses(load_buses)
            self.assertGreater(net.num_vargens,0)

            # vars
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG|pf.BUS_VAR_VDEV|pf.BUS_VAR_VVIO)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P|pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_VARGEN,
                          pf.FLAG_VARS,
                          pf.VARGEN_PROP_ANY,
                          pf.VARGEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_ANY,
                          pf.BRANCH_VAR_RATIO|pf.BRANCH_VAR_RATIO_DEV|pf.BRANCH_VAR_PHASE)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_ANY,
                          pf.SHUNT_VAR_SUSC|pf.SHUNT_VAR_SUSC_DEV)
            self.assertEqual(net.num_vars,
                             (6*net.num_buses +
                              2*net.num_gens + 
                              net.num_vargens +
                              4*net.num_branches + 
                              3*net.num_shunts))
            
            # Add some interesting vargen values
            for vargen in net.var_generators:
                vargen.P = 2*vargen.index
                vargen.P_max = 3*vargen.index

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
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],vargen.P)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b)
                self.assertEqual(x[shunt.index_y],0.)
                self.assertEqual(x[shunt.index_z],0.)
                
            # upper limits
            x = net.get_var_values(pf.UPPER_LIMITS)
            self.assertEqual(x.size,net.num_vars)
            for bus in net.buses:
                self.assertEqual(x[bus.index_v_mag],bus.v_max)
                self.assertEqual(x[bus.index_v_ang],pf.PI)
                self.assertEqual(x[bus.index_y],pf.INF)
                self.assertEqual(x[bus.index_z],pf.INF)
                self.assertEqual(x[bus.index_vl],pf.INF)
                self.assertEqual(x[bus.index_vh],pf.INF)
            for br in net.branches:
                self.assertEqual(x[br.index_ratio],br.ratio_max)
                self.assertEqual(x[br.index_ratio_y],pf.INF)
                self.assertEqual(x[br.index_ratio_z],pf.INF)
                self.assertEqual(x[br.index_phase],br.phase_max)
            for gen in net.generators:
                self.assertEqual(x[gen.index_P],gen.P_max)
                self.assertEqual(x[gen.index_Q],gen.Q_max)
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],vargen.P_max)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_max)
                self.assertEqual(x[shunt.index_y],pf.INF)
                self.assertEqual(x[shunt.index_z],pf.INF)

            # lower limits
            x = net.get_var_values(pf.LOWER_LIMITS)
            self.assertEqual(x.size,net.num_vars)
            for bus in net.buses:
                self.assertEqual(x[bus.index_v_mag],bus.v_min)
                self.assertEqual(x[bus.index_v_ang],-pf.PI)
                self.assertEqual(x[bus.index_y],-pf.INF)
                self.assertEqual(x[bus.index_z],-pf.INF)
                self.assertEqual(x[bus.index_vl],-pf.INF)
                self.assertEqual(x[bus.index_vh],-pf.INF)
            for br in net.branches:
                self.assertEqual(x[br.index_ratio],br.ratio_min)
                self.assertEqual(x[br.index_ratio_y],-pf.INF)
                self.assertEqual(x[br.index_ratio_z],-pf.INF)
                self.assertEqual(x[br.index_phase],br.phase_min)
            for gen in net.generators:
                self.assertEqual(x[gen.index_P],gen.P_min)
                self.assertEqual(x[gen.index_Q],gen.Q_min)
            for vargen in net.var_generators:
                self.assertEqual(x[vargen.index_P],0.)
            for shunt in net.shunts:
                self.assertEqual(x[shunt.index_b],shunt.b_min)
                self.assertEqual(x[shunt.index_y],-pf.INF)
                self.assertEqual(x[shunt.index_z],-pf.INF)

    def test_vargen_P_sigma(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)

            spread = 2
            corr = 0.1
            P_MIN = 1e-5
        
            # Add renewable sources (at load buses)
            self.assertEqual(net.num_vargens,0)
            load_buses = net.get_load_buses()
            vargen_array = pf.VarGeneratorArray(len(load_buses))
            net.set_vargen_array(vargen_array)
            net.set_vargen_buses(load_buses)
            self.assertEqual(net.num_vargens,len([b for b in net.buses if b.loads]))
            for vg in net.var_generators:
                self.assertEqual(vg.P,0)
                self.assertEqual(vg.P_max,0)
                self.assertEqual(vg.P_std,0)
                self.assertGreater(len(vg.bus.loads),0)
                vg.P_max = np.maximum(vg.bus.get_total_load_P(),P_MIN)
                vg.P_std = 0.25*vg.P_max
                self.assertGreater(vg.P_max,0)
                self.assertGreater(vg.P_std,0)
            self.assertGreaterEqual(sum([vg.P_max for vg in net.var_generators]),
                                    sum([l.P for l in net.loads]))

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
                
                
