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
            bus_k7 = br7.bus_k
            bus_k3 = br3.bus_k
            bus_m7 = br7.bus_m
            bus_m3 = br3.bus_m
            bus_k7_degree = br7.bus_k.degree
            bus_m7_degree = br7.bus_m.degree
            bus_k3_degree = br3.bus_k.degree
            bus_m3_degree = br3.bus_m.degree
            if net.num_generators > 5:
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
            if net.num_generators > 5:
                cont.add_gen_outage(g5)
                self.assertEqual(cont.num_gen_outages,2)
                self.assertTrue(cont.has_gen_outage(g5))
            self.assertTrue(cont.has_gen_outage(g0))
            cont.add_branch_outage(br3)
            self.assertEqual(cont.num_branch_outages,2)
            self.assertTrue(cont.has_branch_outage(br3))
            self.assertTrue(cont.has_branch_outage(br7))

            # apply
            self.assertEqual(len([g for g in net.generators if not g.outage]),net.num_generators)
            self.assertEqual(len([b for b in net.branches if not b.outage]),net.num_branches)
            self.assertEqual(len([g for g in net.generators if g.outage]),0)
            self.assertEqual(len([b for b in net.branches if b.outage]),0)
            cont.apply()
            if net.num_generators > 5:
                self.assertEqual(net.get_num_gens_not_on_outage(),net.num_generators-2)
                self.assertEqual(len([g for g in net.generators if g.outage]),2)
            else:
                self.assertEqual(net.get_num_gens_not_on_outage(),net.num_generators-1)
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
                        self.assertFalse(g.index in [y.index for y in bus0.generators])
                        self.assertFalse(g.index in [y.index for y in reg_bus0.reg_generators])
                    elif g.index == 5:
                        self.assertFalse(g.index in [y.index for y in bus5.generators])
                        self.assertFalse(g.index in [y.index for y in reg_bus5.reg_generators])
                else:
                    self.assertFalse(g.is_on_outage())
                    self.assertFalse(g.outage)
            if bus_k7 == bus_k3 or bus_k7 == bus_m3:
                self.assertEqual(bus_k7_degree-2,bus_k7.degree)
            else:
                self.assertEqual(bus_k7_degree-1,bus_k7.degree)
            if bus_m7 == bus_k3 or bus_m7 == bus_m3:
                self.assertEqual(bus_m7_degree-2,bus_m7.degree)
            else:
                self.assertEqual(bus_m7_degree-1,bus_m7.degree)
            if bus_k3 == bus_k7 or bus_k3 == bus_m7:
                self.assertEqual(bus_k3_degree-2,bus_k3.degree)
            else:
                self.assertEqual(bus_k3_degree-1,bus_k3.degree)
            if bus_m3 == bus_k7 or bus_m3 == bus_m7:
                self.assertEqual(bus_m3_degree-2,bus_m3.degree)
            else:
                self.assertEqual(bus_m3_degree-1,bus_m3.degree)
            for b in net.branches:
                if b.index == 0 or b.index == 1:
                    self.assertTrue(b.is_line() or b.is_fixed_tran())
                    self.assertTrue(b.is_on_outage())
                    self.assertTrue(b.outage)
                    self.assertRaises(pf.BusError,lambda x: x.bus_k,b)
                    self.assertRaises(pf.BusError,lambda x: x.bus_m,b)
                    self.assertRaises(pf.BusError,lambda x: x.reg_bus,b)
                    if b.index == 0:
                        self.assertFalse(b.index in [y.index for y in bus_k7.branches_k])
                        self.assertFalse(b.index in [y.index for y in bus_k7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_k7.branches_m])
                        self.assertFalse(b.index in [y.index for y in bus_m7.branches_k])
                        self.assertFalse(b.index in [y.index for y in bus_m7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_m7.branches_m])
                    elif b.index == 1:
                        self.assertFalse(b.index in [y.index for y in bus_k3.branches_k])
                        self.assertFalse(b.index in [y.index for y in bus_k3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_k3.branches_m])
                        self.assertFalse(b.index in [y.index for y in bus_m3.branches_k])
                        self.assertFalse(b.index in [y.index for y in bus_m3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_m3.branches_m])
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
            if net.num_generators > 5:
                self.assertTrue(net.get_gen(5).outage)
            self.assertEqual(cont2.num_branch_outages,1)
            self.assertEqual(cont2.num_gen_outages,0)

            # clear
            cont.clear()
            self.assertTrue(net.get_branch(2).outage)
            self.assertFalse(net.get_branch(1).outage)
            self.assertFalse(net.get_branch(0).outage)
            self.assertFalse(net.get_gen(0).outage)
            if net.num_generators > 5:
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
                        self.assertTrue(g.index in [y.index for y in bus0.generators])
                        self.assertTrue(g.index in [y.index for y in reg_bus0.reg_generators])
                    elif g.index == 5:
                        self.assertEqual(g.bus.index,bus5.index)
                        self.assertEqual(g.reg_bus.index,reg_bus5.index)
                        self.assertTrue(g.index in [y.index for y in bus5.generators])
                        self.assertTrue(g.index in [y.index for y in reg_bus5.reg_generators])
                else:
                    self.assertFalse(g.is_on_outage())
                    self.assertFalse(g.outage)
            self.assertEqual(bus_k7_degree,bus_k7.degree)
            self.assertEqual(bus_m7_degree,bus_m7.degree)
            self.assertEqual(bus_k3_degree,bus_k3.degree)
            self.assertEqual(bus_m3_degree,bus_m3.degree)
            for b in net.branches:
                if b.index == 0 or b.index == 1:
                    self.assertFalse(b.is_on_outage())
                    self.assertFalse(b.outage)
                    if b.index == 0:
                        self.assertEqual(b.bus_k.index,bus_k7.index)
                        self.assertEqual(b.bus_m.index,bus_m7.index)
                        self.assertTrue(b.index in [y.index for y in bus_k7.branches_k])
                        self.assertTrue(b.index in [y.index for y in bus_k7.branches])
                        self.assertFalse(b.index in [y.index for y in bus_k7.branches_m])
                        self.assertFalse(b.index in [y.index for y in bus_m7.branches_k])
                        self.assertTrue(b.index in [y.index for y in bus_m7.branches])
                        self.assertTrue(b.index in [y.index for y in bus_m7.branches_m])
                    elif b.index == 1:
                        self.assertEqual(b.bus_k.index,bus_k3.index)
                        self.assertEqual(b.bus_m.index,bus_m3.index)
                        self.assertTrue(b.index in [y.index for y in bus_k3.branches_k])
                        self.assertTrue(b.index in [y.index for y in bus_k3.branches])
                        self.assertFalse(b.index in [y.index for y in bus_k3.branches_m])
                        self.assertFalse(b.index in [y.index for y in bus_m3.branches_k])
                        self.assertTrue(b.index in [y.index for y in bus_m3.branches])
                        self.assertTrue(b.index in [y.index for y in bus_m3.branches_m])

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

                gens = bus.generators
                if reg_bus is not None:
                    reg_gens = reg_bus.reg_generators

                cont.apply()
                cont.apply()
                cont.apply()

                self.assertTrue(gen.is_on_outage())
                self.assertFalse(gen.is_regulator())
                self.assertRaises(pf.BusError,lambda x: x.bus,gen)
                self.assertRaises(pf.BusError,lambda x: x.reg_bus,gen)
                self.assertFalse(gen.index in [x.index for x in bus.generators])
                self.assertTrue(gen.index in [x.index for x in gens])
                if reg_bus is not None:
                    self.assertFalse(gen.is_regulator())
                    self.assertFalse(gen.index in [x.index for x in reg_bus.reg_generators])
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
                self.assertTrue(gen.index in [x.index for x in bus.generators])
                self.assertEqual(len(gens),len(bus.generators))
                self.assertTrue(set(map(lambda x: x.index,bus.generators)) ==
                                set(map(lambda x: x.index,gens)))
                if reg_bus is not None:
                    self.assertTrue(gen.index in [x.index for x in reg_gens])
                    self.assertTrue(gen.index in [x.index for x in reg_bus.reg_generators])
                    self.assertEqual(len(reg_gens),len(reg_bus.reg_generators))
                    self.assertTrue(set(map(lambda x: x.index,reg_bus.reg_generators)) ==
                                    set(map(lambda x: x.index,reg_gens)))

            # branch single contingencies
            for br in net.branches:

                self.assertFalse(br.is_on_outage())
                cont = pf.Contingency()
                cont.add_branch_outage(br)
                self.assertFalse(br.is_on_outage())

                bus_k = br.bus_k
                bus_m = br.bus_m
                reg_bus = br.reg_bus if br.is_tap_changer_v() else None

                bus_k_branches_k = bus_k.branches_k
                bus_k_branches_m = bus_k.branches_m
                bus_k_branches = bus_k.branches
                bus_m_branches_k = bus_m.branches_k
                bus_m_branches_m = bus_m.branches_m
                bus_m_branches = bus_m.branches

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
                self.assertRaises(pf.BusError,lambda x: x.bus_k,br)
                self.assertRaises(pf.BusError,lambda x: x.bus_m,br)
                self.assertRaises(pf.BusError,lambda x: x.reg_bus,br)
                self.assertFalse(br.index in [x.index for x in bus_k.branches_k])
                self.assertFalse(br.index in [x.index for x in bus_k.branches_m])
                self.assertFalse(br.index in [x.index for x in bus_k.branches])
                self.assertFalse(br.index in [x.index for x in bus_m.branches_k])
                self.assertFalse(br.index in [x.index for x in bus_m.branches_m])
                self.assertFalse(br.index in [x.index for x in bus_m.branches])

                self.assertTrue(br.index in [x.index for x in bus_k_branches_k])
                self.assertFalse(br.index in [x.index for x in bus_k_branches_m])
                self.assertTrue(br.index in [x.index for x in bus_k_branches])
                self.assertFalse(br.index in [x.index for x in bus_m_branches_k])
                self.assertTrue(br.index in [x.index for x in bus_m_branches_m])
                self.assertTrue(br.index in [x.index for x in bus_m_branches])

                if reg_bus is not None:
                    self.assertFalse(br.is_tap_changer())
                    self.assertTrue(br.is_fixed_tran())
                    self.assertFalse(br.index in [x.index for x in reg_bus.reg_trans])
                    self.assertTrue(br.index in [x.index for x in reg_trans])

                cont.clear()
                cont.clear()
                cont.clear()

                self.assertFalse(br.is_on_outage())
                self.assertEqual(br.bus_k.index,bus_k.index)
                self.assertEqual(br.bus_m.index,bus_m.index)
                if reg_bus is not None:
                    self.assertTrue(br.is_tap_changer_v())
                    self.assertTrue(br.is_tap_changer())
                    self.assertEqual(br.reg_bus.index,reg_bus.index)

                self.assertTrue(br.index in [x.index for x in bus_k_branches_k])
                self.assertTrue(br.index in [x.index for x in bus_k_branches])
                self.assertFalse(br.index in [x.index for x in bus_k_branches_m])
                self.assertTrue(br.index in [x.index for x in bus_k.branches_k])
                self.assertTrue(br.index in [x.index for x in bus_k.branches])
                self.assertFalse(br.index in [x.index for x in bus_k.branches_m])
                self.assertEqual(len(bus_k_branches_k),len(bus_k.branches_k))
                self.assertEqual(len(bus_k_branches_m),len(bus_k.branches_m))
                self.assertEqual(len(bus_k_branches),len(bus_k.branches))

                self.assertFalse(br.index in [x.index for x in bus_m_branches_k])
                self.assertTrue(br.index in [x.index for x in bus_m_branches])
                self.assertTrue(br.index in [x.index for x in bus_m_branches_m])
                self.assertFalse(br.index in [x.index for x in bus_m.branches_k])
                self.assertTrue(br.index in [x.index for x in bus_m.branches])
                self.assertTrue(br.index in [x.index for x in bus_m.branches_m])
                self.assertEqual(len(bus_m_branches_k),len(bus_m.branches_k))
                self.assertEqual(len(bus_m_branches_m),len(bus_m.branches_m))
                self.assertEqual(len(bus_m_branches),len(bus_m.branches))

                self.assertTrue(set(map(lambda x: x.index,bus_k_branches_k)) ==
                                set(map(lambda x: x.index,bus_k.branches_k)))
                self.assertTrue(set(map(lambda x: x.index,bus_k_branches_m)) ==
                                set(map(lambda x: x.index,bus_k.branches_m)))
                self.assertTrue(set(map(lambda x: x.index,bus_k_branches)) ==
                                set(map(lambda x: x.index,bus_k.branches)))

                self.assertTrue(set(map(lambda x: x.index,bus_m_branches_k)) ==
                                set(map(lambda x: x.index,bus_m.branches_k)))
                self.assertTrue(set(map(lambda x: x.index,bus_m_branches_m)) ==
                                set(map(lambda x: x.index,bus_m.branches_m)))
                self.assertTrue(set(map(lambda x: x.index,bus_m_branches)) ==
                                set(map(lambda x: x.index,bus_m.branches)))

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
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            self.assertEqual(net.num_vars,net.num_generators)

            # pre contingency
            net.update_properties()
            gen_cost_base = net.gen_P_cost
            func = pf.Function('generation cost',1.,net)
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

                if br.bus_k.degree == 1 or br.bus_m.degree == 1:
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
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            net.set_flags('bus',
                          'variable',
                          'any',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('branch',
                          'variable',
                          'tap changer',
                          'tap ratio')
            net.set_flags('shunt',
                          'variable',
                          'switching - v',
                          'susceptance')
            self.assertEqual(net.num_vars,
                             (net.num_generators +
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
            constr = pf.Constraint('AC power balance',net)
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

                if br.bus_k.degree == 1 or br.bus_m.degree == 1:
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
            net.set_flags('generator',
                          'variable',
                          'any',
                          'active power')
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            self.assertEqual(net.num_vars,
                             (net.num_generators +
                              net.num_buses-net.get_num_slack_buses() +
                              net.get_num_phase_shifters()))

            # pre contingency
            constr = pf.Constraint('DC power balance',net)
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

                bus_k = br.bus_k
                bus_m = br.bus_m

                if br.bus_k.degree == 1 or br.bus_m.degree == 1:
                    continue

                Pkm = br.P_km_DC

                cont = pf.Contingency()
                cont.add_branch_outage(br)
                cont.apply()

                constr.del_matvec()
                constr.analyze()
                constr.eval(net.get_var_values())

                self.assertLess(np.abs((A*x-b)[bus_k.index]-
                                       ((constr.A*x-constr.b)[bus_k.index]-Pkm)),1e-8)
                self.assertLess(np.abs((A*x-b)[bus_m.index]-
                                       ((constr.A*x-constr.b)[bus_m.index]+Pkm)),1e-8)

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

            for branch in net.branches:
                if branch.ratingA == 0.:
                    branch.ratingA = 100.

            # variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('branch',
                          'variable',
                          'phase shifter',
                          'phase shift')
            self.assertEqual(net.num_vars,
                             (net.num_buses-net.get_num_slack_buses() +
                              net.get_num_phase_shifters()))

            # pre contingency
            constr = pf.Constraint('DC branch flow limits',net)
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

                if br.bus_k.degree == 1 or br.bus_m.degree == 1:
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
            net.set_flags('generator',
                          'variable',
                          'not on outage',
                          'active power')
            net.set_flags('branch',
                          'variable',
                          'not on outage',
                          'tap ratio')
            self.assertEqual(net.num_vars,
                             (net.num_generators-1+net.num_branches-2))
            for gen in net.generators:
                if gen.index != 0:
                    self.assertTrue(gen.has_flags('variable','active power'))
                else:
                    self.assertFalse(gen.has_flags('variable','active power'))
            for br in net.branches:
                if br.index not in [0,1]:
                    self.assertTrue(br.has_flags('variable','tap ratio'))
                else:
                    self.assertFalse(br.has_flags('variable','tap ratio'))


    def tearDown(self):

        pass
