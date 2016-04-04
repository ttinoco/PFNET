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
from scipy.sparse import coo_matrix, bmat, triu

TEST_BRANCHES = 100
TEST_GENS = 100
TEST_BUSES = 20

class TestContingency(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_construction(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            # outage init
            for gen in net.generators:
                self.assertFalse(gen.is_on_outage())
                self.assertFalse(gen.outage)
            for branch in net.branches:
                self.assertFalse(branch.is_on_outage())
                self.assertFalse(branch.outage)            

            # outage set
            gen = net.get_gen(0)
            branch = net.get_branch(0)
            def s1():
                gen.outage = True
            def s2():
                branch.outage = True
            self.assertRaises(AttributeError,s1)
            self.assertRaises(AttributeError,s2)

            # outages at construction
            c0 = pf.Contingency([net.get_gen(0)],
                                [net.get_branch(1)])
            self.assertTrue(c0.has_gen_outage(net.get_gen(0)))
            self.assertTrue(c0.has_branch_outage(net.get_branch(1)))
            c1 = pf.Contingency(gens=[net.get_gen(0)],
                                branches=[net.get_branch(2)])
            self.assertTrue(c1.has_gen_outage(net.get_gen(0)))
            self.assertTrue(c1.has_branch_outage(net.get_branch(2)))
            
            # contingency
            g0 = net.get_gen(0)
            bus0 = g0.bus
            reg_bus0 = g0.reg_bus
            br7 = net.get_branch(0)
            br3 = net.get_branch(1)
            bus_from7 = br7.bus_from
            bus_from3 = br3.bus_from
            bus_to7 = br7.bus_to
            bus_to3 = br3.bus_to
            bus_from7_degree = br7.bus_from.degree
            bus_to7_degree = br7.bus_to.degree
            bus_from3_degree = br3.bus_from.degree
            bus_to3_degree = br3.bus_to.degree
            if net.num_gens > 5:
                g5 = net.get_gen(5)
                bus5 = g5.bus
                reg_bus5 = g5.reg_bus
            cont = pf.Contingency()
            self.assertEqual(cont.num_gen_outages,0)
            self.assertEqual(cont.num_branch_outages,0)
            cont.add_gen_outage(g0)
            cont.add_branch_outage(br7)
            self.assertEqual(cont.num_gen_outages,1)
            self.assertEqual(cont.num_branch_outages,1)
            cont.add_gen_outage(g0)
            cont.add_branch_outage(br7)
            self.assertEqual(cont.num_gen_outages,1)
            self.assertEqual(cont.num_branch_outages,1)
            if net.num_gens > 5:
                cont.add_gen_outage(g5)
                self.assertEqual(cont.num_gen_outages,2)
                self.assertTrue(cont.has_gen_outage(g5))
            self.assertTrue(cont.has_gen_outage(g0))
            cont.add_branch_outage(br3)
            self.assertEqual(cont.num_branch_outages,2)
            self.assertTrue(cont.has_branch_outage(br3))
            self.assertTrue(cont.has_branch_outage(br7))

            # apply
            self.assertEqual(len([g for g in net.generators if not g.outage]),net.num_gens)
            self.assertEqual(len([b for b in net.branches if not b.outage]),net.num_branches)
            self.assertEqual(len([g for g in net.generators if g.outage]),0)
            self.assertEqual(len([b for b in net.branches if b.outage]),0)
            cont.apply()
            if net.num_gens > 5:
                self.assertEqual(net.get_num_gens_not_on_outage(),net.num_gens-2)
                self.assertEqual(len([g for g in net.generators if g.outage]),2)
            else:
                self.assertEqual(net.get_num_gens_not_on_outage(),net.num_gens-1)
                self.assertEqual(len([g for g in net.generators if g.outage]),1)
            self.assertEqual(len([b for b in net.branches if b.outage]),2)
            self.assertEqual(net.get_num_branches_not_on_outage(),net.num_branches-2)
            for g in net.generators:
                if g.index == 0 or g.index == 5:
                    self.assertFalse(g.is_regulator())
                    self.assertTrue(g.is_on_outage())
                    self.assertTrue(g.outage)
                    self.assertRaises(pf.BusError,lambda x: x.bus,g)
                    self.assertRaises(pf.BusError,lambda x: x.reg_bus,g)
                    if g.index == 0:
                        self.assertFalse(g.index in [y.index for y in bus0.gens])
                        self.assertFalse(g.index in [y.index for y in reg_bus0.reg_gens])
                    elif g.index == 5:
                        self.assertFalse(g.index in [y.index for y in bus5.gens])
                        self.assertFalse(g.index in [y.index for y in reg_bus5.reg_gens])
                else:
                    self.assertFalse(g.is_on_outage())
                    self.assertFalse(g.outage)
            if bus_from7 == bus_from3 or bus_from7 == bus_to3:
                self.assertEqual(bus_from7_degree-2,bus_from7.degree)
            else:
                self.assertEqual(bus_from7_degree-1,bus_from7.degree)
            if bus_to7 == bus_from3 or bus_to7 == bus_to3:
                self.assertEqual(bus_to7_degree-2,bus_to7.degree)
            else:
                self.assertEqual(bus_to7_degree-1,bus_to7.degree)
            if bus_from3 == bus_from7 or bus_from3 == bus_to7:
                self.assertEqual(bus_from3_degree-2,bus_from3.degree)
            else:
                self.assertEqual(bus_from3_degree-1,bus_from3.degree)
            if bus_to3 == bus_from7 or bus_to3 == bus_to7:
                self.assertEqual(bus_to3_degree-2,bus_to3.degree)
            else:
                self.assertEqual(bus_to3_degree-1,bus_to3.degree)
            for b in net.branches:
                if b.index == 0 or b.index == 1:
                    self.assertTrue(b.is_line() or b.is_fixed_tran())
                    self.assertTrue(b.is_on_outage())
                    self.assertTrue(b.outage)
                    self.assertRaises(pf.BusError,lambda x: x.bus_from,b)
                    self.assertRaises(pf.BusError,lambda x: x.bus_to,b)
                    self.assertRaises(pf.BusError,lambda x: x.reg_bus,b)
                    if b.index == 0:
                        self.assertFalse(b.index in [y.index for y in bus_from7.branches_from])
                        self.assertFalse(b.index in [y.index for y in bus_from7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_from7.branches_to])
                        self.assertFalse(b.index in [y.index for y in bus_to7.branches_from])
                        self.assertFalse(b.index in [y.index for y in bus_to7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_to7.branches_to])
                    elif b.index == 1:
                        self.assertFalse(b.index in [y.index for y in bus_from3.branches_from])
                        self.assertFalse(b.index in [y.index for y in bus_from3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_from3.branches_to])
                        self.assertFalse(b.index in [y.index for y in bus_to3.branches_from])
                        self.assertFalse(b.index in [y.index for y in bus_to3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_to3.branches_to])
                else:
                    self.assertFalse(b.is_on_outage())
                    self.assertFalse(b.outage)
            cont2 = pf.Contingency()
            cont2.add_branch_outage(net.get_branch(2))
            self.assertFalse(net.get_branch(2).outage)
            cont2.apply()
            self.assertTrue(net.get_branch(2).outage)
            self.assertTrue(net.get_branch(0).outage)
            self.assertTrue(net.get_branch(1).outage)
            self.assertTrue(net.get_gen(0).outage)
            if net.num_gens > 5:
                self.assertTrue(net.get_gen(5).outage)
            self.assertEqual(cont2.num_branch_outages,1)
            self.assertEqual(cont2.num_gen_outages,0)

            # clear
            cont.clear()
            self.assertTrue(net.get_branch(2).outage)
            self.assertFalse(net.get_branch(1).outage)
            self.assertFalse(net.get_branch(0).outage)
            self.assertFalse(net.get_gen(0).outage)
            if net.num_gens > 5:
                self.assertFalse(net.get_gen(5).outage)
            self.assertEqual(len([b for b in net.branches if b.outage]),1)
            cont2.clear()
            self.assertEqual(len([b for b in net.branches if b.outage]),0)
            for g in net.generators:
                if g.index == 0 or g.index == 5:
                    self.assertFalse(g.is_on_outage())
                    self.assertFalse(g.outage)
                    if g.index == 0:
                        self.assertEqual(g.bus.index,bus0.index)
                        self.assertEqual(g.reg_bus.index,reg_bus0.index)
                        self.assertTrue(g.index in [y.index for y in bus0.gens])
                        self.assertTrue(g.index in [y.index for y in reg_bus0.reg_gens])
                    elif g.index == 5:
                        self.assertEqual(g.bus.index,bus5.index)
                        self.assertEqual(g.reg_bus.index,reg_bus5.index)
                        self.assertTrue(g.index in [y.index for y in bus5.gens])
                        self.assertTrue(g.index in [y.index for y in reg_bus5.reg_gens])
                else:
                    self.assertFalse(g.is_on_outage())
                    self.assertFalse(g.outage)
            self.assertEqual(bus_from7_degree,bus_from7.degree)
            self.assertEqual(bus_to7_degree,bus_to7.degree)
            self.assertEqual(bus_from3_degree,bus_from3.degree)
            self.assertEqual(bus_to3_degree,bus_to3.degree)
            for b in net.branches:
                if b.index == 0 or b.index == 1:
                    self.assertFalse(b.is_on_outage())
                    self.assertFalse(b.outage)
                    if b.index == 0:
                        self.assertEqual(b.bus_from.index,bus_from7.index)
                        self.assertEqual(b.bus_to.index,bus_to7.index)
                        self.assertTrue(b.index in [y.index for y in bus_from7.branches_from])
                        self.assertTrue(b.index in [y.index for y in bus_from7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_from7.branches_to])
                        self.assertFalse(b.index in [y.index for y in bus_to7.branches_from])
                        self.assertTrue(b.index in [y.index for y in bus_to7.branches])
                        self.assertTrue(b.index in [y.index for y in bus_to7.branches_to])
                    elif b.index == 1:
                        self.assertEqual(b.bus_from.index,bus_from3.index)
                        self.assertEqual(b.bus_to.index,bus_to3.index)
                        self.assertTrue(b.index in [y.index for y in bus_from3.branches_from])
                        self.assertTrue(b.index in [y.index for y in bus_from3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_from3.branches_to])
                        self.assertFalse(b.index in [y.index for y in bus_to3.branches_from])
                        self.assertTrue(b.index in [y.index for y in bus_to3.branches])
                        self.assertTrue(b.index in [y.index for y in bus_to3.branches_to])
                
            # do it again
            net.clear_properties()
            net.load(case)
            net.clear_flags()
            
            # generator single contingencies
            for gen in net.generators:

                self.assertFalse(gen.is_on_outage())
                cont = pf.Contingency()
                cont.add_gen_outage(gen)                
                self.assertFalse(gen.is_on_outage())

                bus = gen.bus
                reg_bus = gen.reg_bus if gen.is_regulator() else None

                gens = bus.gens
                if reg_bus is not None:
                    reg_gens = reg_bus.reg_gens
                    
                cont.apply()
                cont.apply()
                cont.apply()
                
                self.assertTrue(gen.is_on_outage())
                self.assertFalse(gen.is_regulator())
                self.assertRaises(pf.BusError,lambda x: x.bus,gen)
                self.assertRaises(pf.BusError,lambda x: x.reg_bus,gen)
                self.assertFalse(gen.index in [x.index for x in bus.gens])
                self.assertTrue(gen.index in [x.index for x in gens])
                if reg_bus is not None:
                    self.assertFalse(gen.is_regulator())
                    self.assertFalse(gen.index in [x.index for x in reg_bus.reg_gens])
                    self.assertTrue(gen.index in [x.index for x in reg_gens])

                cont.clear()
                cont.clear()
                cont.clear()

                self.assertFalse(gen.is_on_outage())
                self.assertEqual(gen.bus.index,bus.index)
                if reg_bus is not None:
                    self.assertTrue(gen.is_regulator())
                    self.assertEqual(gen.reg_bus.index,reg_bus.index)
                self.assertTrue(gen.index in [x.index for x in gens])
                self.assertTrue(gen.index in [x.index for x in bus.gens])
                self.assertEqual(len(gens),len(bus.gens))
                self.assertTrue(set(map(lambda x: x.index,bus.gens)) == 
                                set(map(lambda x: x.index,gens)))
                if reg_bus is not None:
                    self.assertTrue(gen.index in [x.index for x in reg_gens])
                    self.assertTrue(gen.index in [x.index for x in reg_bus.reg_gens])
                    self.assertEqual(len(reg_gens),len(reg_bus.reg_gens))
                    self.assertTrue(set(map(lambda x: x.index,reg_bus.reg_gens)) == 
                                    set(map(lambda x: x.index,reg_gens)))

            # branch single contingencies
            for br in net.branches:

                self.assertFalse(br.is_on_outage())
                cont = pf.Contingency()
                cont.add_branch_outage(br)                
                self.assertFalse(br.is_on_outage())

                bus_from = br.bus_from
                bus_to = br.bus_to
                reg_bus = br.reg_bus if br.is_tap_changer_v() else None

                bus_from_branches_from = bus_from.branches_from
                bus_from_branches_to = bus_from.branches_to
                bus_from_branches = bus_from.branches
                bus_to_branches_from = bus_to.branches_from
                bus_to_branches_to = bus_to.branches_to
                bus_to_branches = bus_to.branches
                
                if reg_bus is not None:
                    reg_trans = reg_bus.reg_trans

                br_types = [br.is_line(),
                            br.is_fixed_tran(),
                            br.is_phase_shifter(),
                            br.is_tap_changer(),
                            br.is_tap_changer_v(),
                            br.is_tap_changer_Q()]

                cont.apply()
                cont.apply()
                cont.apply()
                
                self.assertTrue(br.is_on_outage())
                self.assertFalse(br.is_tap_changer())
                if br_types[0]:
                    self.assertTrue(br.is_line())
                else:
                    self.assertTrue(br.is_fixed_tran())
                self.assertRaises(pf.BusError,lambda x: x.bus_from,br)
                self.assertRaises(pf.BusError,lambda x: x.bus_to,br)
                self.assertRaises(pf.BusError,lambda x: x.reg_bus,br)
                self.assertFalse(br.index in [x.index for x in bus_from.branches_from])
                self.assertFalse(br.index in [x.index for x in bus_from.branches_to])
                self.assertFalse(br.index in [x.index for x in bus_from.branches])
                self.assertFalse(br.index in [x.index for x in bus_to.branches_from])
                self.assertFalse(br.index in [x.index for x in bus_to.branches_to])
                self.assertFalse(br.index in [x.index for x in bus_to.branches])

                self.assertTrue(br.index in [x.index for x in bus_from_branches_from])
                self.assertFalse(br.index in [x.index for x in bus_from_branches_to])
                self.assertTrue(br.index in [x.index for x in bus_from_branches])
                self.assertFalse(br.index in [x.index for x in bus_to_branches_from])
                self.assertTrue(br.index in [x.index for x in bus_to_branches_to])
                self.assertTrue(br.index in [x.index for x in bus_to_branches])

                if reg_bus is not None:
                    self.assertFalse(br.is_tap_changer())
                    self.assertTrue(br.is_fixed_tran())
                    self.assertFalse(br.index in [x.index for x in reg_bus.reg_trans])
                    self.assertTrue(br.index in [x.index for x in reg_trans])

                cont.clear()
                cont.clear()
                cont.clear()

                self.assertFalse(br.is_on_outage())
                self.assertEqual(br.bus_from.index,bus_from.index)
                self.assertEqual(br.bus_to.index,bus_to.index)
                if reg_bus is not None:
                    self.assertTrue(br.is_tap_changer_v())
                    self.assertTrue(br.is_tap_changer())
                    self.assertEqual(br.reg_bus.index,reg_bus.index)

                self.assertTrue(br.index in [x.index for x in bus_from_branches_from])
                self.assertTrue(br.index in [x.index for x in bus_from_branches])
                self.assertFalse(br.index in [x.index for x in bus_from_branches_to])
                self.assertTrue(br.index in [x.index for x in bus_from.branches_from])
                self.assertTrue(br.index in [x.index for x in bus_from.branches])
                self.assertFalse(br.index in [x.index for x in bus_from.branches_to])
                self.assertEqual(len(bus_from_branches_from),len(bus_from.branches_from))
                self.assertEqual(len(bus_from_branches_to),len(bus_from.branches_to))
                self.assertEqual(len(bus_from_branches),len(bus_from.branches))

                self.assertFalse(br.index in [x.index for x in bus_to_branches_from])
                self.assertTrue(br.index in [x.index for x in bus_to_branches])
                self.assertTrue(br.index in [x.index for x in bus_to_branches_to])
                self.assertFalse(br.index in [x.index for x in bus_to.branches_from])
                self.assertTrue(br.index in [x.index for x in bus_to.branches])
                self.assertTrue(br.index in [x.index for x in bus_to.branches_to])
                self.assertEqual(len(bus_to_branches_from),len(bus_to.branches_from))
                self.assertEqual(len(bus_to_branches_to),len(bus_to.branches_to))
                self.assertEqual(len(bus_to_branches),len(bus_to.branches))

                self.assertTrue(set(map(lambda x: x.index,bus_from_branches_from)) == 
                                set(map(lambda x: x.index,bus_from.branches_from)))
                self.assertTrue(set(map(lambda x: x.index,bus_from_branches_to)) == 
                                set(map(lambda x: x.index,bus_from.branches_to)))
                self.assertTrue(set(map(lambda x: x.index,bus_from_branches)) == 
                                set(map(lambda x: x.index,bus_from.branches)))

                self.assertTrue(set(map(lambda x: x.index,bus_to_branches_from)) == 
                                set(map(lambda x: x.index,bus_to.branches_from)))
                self.assertTrue(set(map(lambda x: x.index,bus_to_branches_to)) == 
                                set(map(lambda x: x.index,bus_to.branches_to)))
                self.assertTrue(set(map(lambda x: x.index,bus_to_branches)) == 
                                set(map(lambda x: x.index,bus_to.branches)))

                if reg_bus is not None:
                    self.assertTrue(br.index in [x.index for x in reg_trans])
                    self.assertTrue(br.index in [x.index for x in reg_bus.reg_trans])
                    self.assertEqual(len(reg_trans),len(reg_bus.reg_trans))
                    self.assertTrue(set(map(lambda x: x.index,reg_bus.reg_trans)) == 
                                    set(map(lambda x: x.index,reg_trans)))
                    
                self.assertTupleEqual(tuple(br_types),
                                      (br.is_line(),
                                       br.is_fixed_tran(),
                                       br.is_phase_shifter(),
                                       br.is_tap_changer(),
                                       br.is_tap_changer_v(),
                                       br.is_tap_changer_Q()))

    def test_gen_cost(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            # variables
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            self.assertEqual(net.num_vars,net.num_gens)

            # pre contingency
            net.update_properties()
            gen_cost_base = net.gen_P_cost
            func = pf.Function(pf.FUNC_TYPE_GEN_COST,1.,net)
            func.analyze()
            func.eval(net.get_var_values())
            phi_base = func.phi
            gphi_base = func.gphi.copy()
            Hphi_base = func.Hphi.copy()
            self.assertEqual(phi_base,gen_cost_base)
            self.assertTupleEqual(gphi_base.shape,(net.num_vars,))
            self.assertEqual(Hphi_base.nnz,net.num_vars)

            # gen outages
            counter = 0
            for gen in net.generators:
                
                cont = pf.Contingency()
                cont.add_gen_outage(gen)
                cont.apply()

                func.del_matvec()
                func.analyze()                
                func.eval(net.get_var_values())

                net.update_properties()

                # value
                self.assertLess(np.abs(phi_base-gen.P_cost-func.phi),1e-8)
                self.assertLess(np.abs(gen_cost_base-gen.P_cost-net.gen_P_cost),1e-8)

                # grad
                gphi = func.gphi
                self.assertEqual(gphi[gen.index_P],0.)
                gphi[gen.index_P] = gphi_base[gen.index_P]
                self.assertLess(np.linalg.norm(gphi-gphi_base,np.inf),1e-8)

                # Hessian
                Hphi = func.Hphi
                self.assertTrue(np.all(Hphi.row != gen.index_P))
                self.assertTrue(np.all(Hphi.col != gen.index_P))

                cont.clear()
                counter += 1
                if counter > TEST_GENS:
                    break

            # branch outages
            counter = 0
            for br in net.branches:
                
                if br.bus_from.degree == 1 or br.bus_to.degree == 1:
                    continue
                
                cont = pf.Contingency()
                cont.add_branch_outage(br)
                cont.apply()

                func.del_matvec()
                func.analyze()                
                func.eval(net.get_var_values())

                net.update_properties()

                # value
                self.assertLess(np.abs(phi_base-func.phi),1e-8)
                self.assertLess(np.abs(gen_cost_base-net.gen_P_cost),1e-8)

                # grad
                self.assertLess(np.linalg.norm(func.gphi-gphi_base,np.inf),1e-8)

                # Hessian
                E = func.Hphi-Hphi_base
                self.assertEqual(E.nnz,0)

                cont.clear()
                counter += 1
                if counter > TEST_BRANCHES:
                    break

    def test_pf(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            # variables
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_REG,
                          pf.GEN_VAR_Q)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_ANY,
                          pf.BUS_VAR_VMAG|pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_TAP_CHANGER,
                          pf.BRANCH_VAR_RATIO)
            net.set_flags(pf.OBJ_SHUNT,
                          pf.FLAG_VARS,
                          pf.SHUNT_PROP_SWITCHED_V,
                          pf.SHUNT_VAR_SUSC)
            self.assertEqual(net.num_vars,
                             (net.num_gens +
                              net.get_num_reg_gens() +
                              2*net.num_buses +
                              net.get_num_tap_changers() +
                              net.get_num_switched_shunts()))

            # pre contingency
            net.update_properties()
            mismatches = np.zeros(2*net.num_buses)
            for bus in net.buses:
                mismatches[bus.index_P] = bus.P_mis
                mismatches[bus.index_Q] = bus.Q_mis
            constr = pf.Constraint(pf.CONSTR_TYPE_PF,net)
            constr.analyze()
            constr.eval(net.get_var_values())
            f = constr.f.copy()
            self.assertLess(np.linalg.norm(f-mismatches),1e-8)

            # gen outages
            counter = 0
            for gen in net.generators:
 
                bus = gen.bus

                cont = pf.Contingency()
                cont.add_gen_outage(gen)
                cont.apply()

                constr.del_matvec()
                constr.analyze()                
                constr.eval(net.get_var_values())
                
                net.update_properties()

                self.assertLess(np.abs(constr.f[bus.index_P]-bus.P_mis),1e-8)
                self.assertLess(np.abs(constr.f[bus.index_Q]-bus.Q_mis),1e-8)
                self.assertLess(np.abs(f[bus.index_P]-constr.f[bus.index_P]-gen.P),1e-8)
                self.assertLess(np.abs(f[bus.index_Q]-constr.f[bus.index_Q]-gen.Q),1e-8)
                self.assertLess(np.abs(mismatches[bus.index_P]-bus.P_mis-gen.P),1e-8)
                self.assertLess(np.abs(mismatches[bus.index_Q]-bus.Q_mis-gen.Q),1e-8)

                counter1 = 0
                for bus1 in net.buses:
                    if bus != bus1:
                        self.assertLess(np.abs(constr.f[bus1.index_P]-f[bus1.index_P]),1e-8)
                        self.assertLess(np.abs(constr.f[bus1.index_Q]-f[bus1.index_Q]),1e-8)
                        self.assertLess(np.abs(bus1.P_mis-mismatches[bus1.index_P]),1e-8)
                        self.assertLess(np.abs(bus1.Q_mis-mismatches[bus1.index_Q]),1e-8)
                        counter1 += 1
                        if counter1 > TEST_BUSES:
                            break

                cont.clear()
                counter += 1
                if counter > TEST_GENS:
                    break

            # branch outages
            counter = 0
            for br in net.branches:
                
                if br.bus_from.degree == 1 or br.bus_to.degree == 1:
                    continue
                
                cont = pf.Contingency()
                cont.add_branch_outage(br)
                cont.apply()

                constr.del_matvec()
                constr.analyze()                
                constr.eval(net.get_var_values())

                net.update_properties()
               
                # NEED TO TEST
                #*************
                
                cont.clear()
                counter += 1
                if counter > TEST_BRANCHES:
                    break

    def test_dcpf(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            # variables
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_ANY,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            self.assertEqual(net.num_vars,
                             (net.num_gens +
                              net.num_buses-net.get_num_slack_buses() +
                              net.get_num_phase_shifters()))

            # pre contingency
            constr = pf.Constraint(pf.CONSTR_TYPE_DCPF,net)
            constr.analyze()
            constr.eval(net.get_var_values())
            A = constr.A.copy()
            b = constr.b.copy()
            x = net.get_var_values()
           
            # gen outages
            counter = 0
            for gen in net.generators:
 
                bus = gen.bus

                cont = pf.Contingency()
                cont.add_gen_outage(gen)
                cont.apply()
            
                constr.del_matvec()
                constr.analyze()                
                constr.eval(x)
                
                self.assertFalse(np.all(A.col != gen.index_P))
                self.assertTrue(np.all(constr.A.col != gen.index_P))

                self.assertLess(np.abs((A*x-b)[bus.index]-
                                       (constr.A*x-constr.b)[bus.index]-gen.P),1e-8)

                counter1 = 0
                for bus1 in net.buses:
                    if bus != bus1:
                        self.assertLess(np.abs((A*x-b)[bus1.index]-
                                               (constr.A*x-constr.b)[bus1.index]),1e-8)
                        counter1 += 1
                        if counter1 > TEST_BUSES:
                            break

                cont.clear()
                counter += 1
                if counter > TEST_GENS:
                    break

            # branch outages
            counter = 0
            for br in net.branches:

                bus_from = br.bus_from
                bus_to = br.bus_to
                
                if br.bus_from.degree == 1 or br.bus_to.degree == 1:
                    continue

                Pflow = br.P_flow_DC
                
                cont = pf.Contingency()
                cont.add_branch_outage(br)
                cont.apply()

                constr.del_matvec()
                constr.analyze()                
                constr.eval(net.get_var_values())
                
                self.assertLess(np.abs((A*x-b)[bus_from.index]-
                                       ((constr.A*x-constr.b)[bus_from.index]-Pflow)),1e-8)
                self.assertLess(np.abs((A*x-b)[bus_to.index]-
                                       ((constr.A*x-constr.b)[bus_to.index]+Pflow)),1e-8)
 
                cont.clear()
                counter += 1
                if counter > TEST_BRANCHES:
                    break

    def test_dc_flow_lim(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            # variables
            net.set_flags(pf.OBJ_BUS,
                          pf.FLAG_VARS,
                          pf.BUS_PROP_NOT_SLACK,
                          pf.BUS_VAR_VANG)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_PHASE_SHIFTER,
                          pf.BRANCH_VAR_PHASE)
            self.assertEqual(net.num_vars,
                             (net.num_buses-net.get_num_slack_buses() +
                              net.get_num_phase_shifters()))

            # pre contingency
            constr = pf.Constraint(pf.CONSTR_TYPE_DC_FLOW_LIM,net)
            constr.analyze()
            constr.eval(net.get_var_values())
            G = constr.G.copy()
            l = constr.l.copy()
            u = constr.u.copy()
            x = net.get_var_values()
           
            # gen outages
            counter = 0
            for gen in net.generators:
 
                cont = pf.Contingency()
                cont.add_gen_outage(gen)
                cont.apply()
            
                constr.del_matvec()
                constr.analyze()                
                constr.eval(x)
                
                self.assertEqual((G-constr.G).nnz,0)
                self.assertEqual(np.linalg.norm(l-constr.l),0.)
                self.assertEqual(np.linalg.norm(u-constr.u),0.)

                cont.clear()
                counter += 1
                if counter > TEST_GENS:
                    break

            # branch outages
            counter = 0
            for br in net.branches:
                
                if br.bus_from.degree == 1 or br.bus_to.degree == 1:
                    continue
                
                cont = pf.Contingency()
                cont.add_branch_outage(br)
                cont.apply()

                constr.del_matvec()
                constr.analyze()                
                constr.eval(net.get_var_values())

                lnew = constr.l.copy()
                unew = constr.u.copy()
                Gnew = constr.G.copy()

                self.assertTupleEqual(lnew.shape,(l.size-1,))
                self.assertTupleEqual(unew.shape,(u.size-1,))
                self.assertTupleEqual(Gnew.shape,(G.shape[0]-1,G.shape[1]))

                self.assertLess(np.linalg.norm(lnew-np.hstack((l[:br.index],l[br.index+1:])),np.inf),1e-8)
                self.assertLess(np.linalg.norm(unew-np.hstack((u[:br.index],u[br.index+1:])),np.inf),1e-8)
                indices = G.row != br.index
                row = G.row[indices]
                row = row - 1*(row>br.index)
                col = G.col[indices]
                data = G.data[indices] 
                Gcut = coo_matrix((data,(row,col)),shape=(net.num_branches-1,net.num_vars))
                self.assertEqual((Gnew-Gcut).nnz,0)

                cont.clear()
                counter += 1
                if counter > TEST_BRANCHES:
                    break

    def test_variables(self):

        net = self.net

        for case in test_cases.CASES:

            net.clear_properties()
            net.load(case)
            net.clear_flags()

            cont = pf.Contingency()
            cont.add_gen_outage(net.get_gen(0))
            cont.add_branch_outage(net.get_branch(0))
            cont.add_branch_outage(net.get_branch(1))

            cont.apply()

            # variables
            net.set_flags(pf.OBJ_GEN,
                          pf.FLAG_VARS,
                          pf.GEN_PROP_NOT_OUT,
                          pf.GEN_VAR_P)
            net.set_flags(pf.OBJ_BRANCH,
                          pf.FLAG_VARS,
                          pf.BRANCH_PROP_NOT_OUT,
                          pf.BRANCH_VAR_RATIO)
            self.assertEqual(net.num_vars,
                             (net.num_gens-1+net.num_branches-2))
            for gen in net.generators:
                if gen.index != 0:
                    self.assertTrue(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
                else:
                    self.assertFalse(gen.has_flags(pf.FLAG_VARS,pf.GEN_VAR_P))
            for br in net.branches:
                if br.index not in [0,1]:
                    self.assertTrue(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))
                else:
                    self.assertFalse(br.has_flags(pf.FLAG_VARS,pf.BRANCH_VAR_RATIO))

 
    def tearDown(self):
        
        pass
