#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import pfnet as pf
import unittest
from . import test_cases
import numpy as np

class TestParser(unittest.TestCase):

    def setUp(self):

        pass

    def test_dummy_parser(self):

        p = pf.parsers.DummyParser()
        
        self.assertEqual(p.some_init_data,3)
        
        self.assertRaises(pf.ParserError,p.parse,"foo.mat")
        
        net = p.parse("foo.dummy")
        self.assertTrue(isinstance(net,pf.Network))
        self.assertEqual(net.num_buses,0)

    def test_sys_problem2(self):

        for case in test_cases.CASES:
            if case == '../data/sys_problem2.mat':

                net = pf.ParserMAT().parse(case)

                self.assertEqual(net.base_power,100.)

                self.assertEqual(net.num_buses,3)
                self.assertEqual(net.num_generators,4)
                self.assertEqual(net.num_loads,3)
                self.assertEqual(net.num_branches,3)

                bus1 = net.get_bus_by_number(1)
                bus2 = net.get_bus_by_number(2)
                bus3 = net.get_bus_by_number(3)

                self.assertEqual(bus1.number,1)
                self.assertEqual(bus2.number,2)
                self.assertEqual(bus3.number,3)

                branch13 = net.get_branch(0)
                branch23 = net.get_branch(1)
                branch12 = net.get_branch(2)
                
                self.assertEqual(len(bus1.generators),2)
                self.assertEqual(len(bus2.generators),1)
                self.assertEqual(len(bus3.generators),1)

                self.assertEqual(len(bus1.loads),1)
                self.assertEqual(len(bus2.loads),1)
                self.assertEqual(len(bus3.loads),1)

                gen0 = net.get_generator(3)
                gen1 = net.get_generator(2)
                gen2 = net.get_generator(1) 
                gen3 = net.get_generator(0)

                load0 = net.get_load(2)
                load1 = net.get_load(1)
                load2 = net.get_load(0)

                self.assertEqual(load0.bus,bus1)
                self.assertEqual(load1.bus,bus2)
                self.assertEqual(load2.bus,bus3)
                
                self.assertEqual(gen0.P,50/100.)
                self.assertEqual(gen1.P,40/100.)
                self.assertEqual(gen2.P,30/100.)
                self.assertEqual(gen3.P,20/100.)
                self.assertEqual(gen0.P_max,50)
                self.assertEqual(gen1.P_max,50)
                self.assertEqual(gen2.P_max,25)
                self.assertEqual(gen3.P_max,19)
                self.assertEqual(gen0.bus,bus1)
                self.assertEqual(gen1.bus,bus1)
                self.assertEqual(gen2.bus,bus2)
                self.assertEqual(gen3.bus,bus3)
                for gen in net.generators:
                    self.assertEqual(gen.P_min,0.)

                self.assertEqual(branch13.bus_k.number,bus1.number)
                self.assertEqual(branch13.bus_m.number,bus3.number)

                self.assertEqual(branch23.bus_k.number,bus2.number)
                self.assertEqual(branch23.bus_m.number,bus3.number)

                self.assertEqual(branch12.bus_k.number,bus1.number)
                self.assertEqual(branch12.bus_m.number,bus2.number)

                self.assertEqual(branch13.g,0)
                self.assertLess(abs(branch13.b + 1./0.1),1e-10)
                self.assertEqual(branch13.ratingA,30.95)
                self.assertEqual(branch13.ratingB,30.95)
                self.assertEqual(branch13.ratingC,30.95)

                self.assertEqual(branch23.g,0)
                self.assertLess(abs(branch23.b + 1./0.2),1e-10)
                self.assertEqual(branch23.ratingA,13)
                self.assertEqual(branch23.ratingB,13)
                self.assertEqual(branch23.ratingC,13)

                self.assertEqual(branch12.g,0)
                self.assertLess(abs(branch12.b + 1./0.2),1e-10)
                self.assertEqual(branch12.ratingA,15)
                self.assertEqual(branch12.ratingB,15)
                self.assertEqual(branch12.ratingC,15)

                self.assertEqual(gen0.cost_coeff_Q0,0)
                self.assertEqual(gen0.cost_coeff_Q1,6.*net.base_power)
                self.assertEqual(gen0.cost_coeff_Q2,0.03*(net.base_power**2.))

                self.assertEqual(gen1.cost_coeff_Q0,0)
                self.assertEqual(gen1.cost_coeff_Q1,5.*net.base_power)
                self.assertEqual(gen1.cost_coeff_Q2,0.02*(net.base_power**2.))

                self.assertEqual(gen2.cost_coeff_Q0,0)
                self.assertEqual(gen2.cost_coeff_Q1,12.*net.base_power)
                self.assertEqual(gen2.cost_coeff_Q2,0.06*(net.base_power**2.))

                self.assertEqual(gen3.cost_coeff_Q0,0)
                self.assertEqual(gen3.cost_coeff_Q1,10.*net.base_power)
                self.assertEqual(gen3.cost_coeff_Q2,0.08*(net.base_power**2.))

                # Load utility
                self.assertEqual(load0.util_coeff_Q0,0)
                self.assertEqual(load0.util_coeff_Q1,400.*net.base_power)
                self.assertEqual(load0.util_coeff_Q2,-0.03*(net.base_power**2.))

                self.assertEqual(load1.util_coeff_Q0,0)
                self.assertEqual(load1.util_coeff_Q1,450.*net.base_power)
                self.assertEqual(load1.util_coeff_Q2,-0.02*(net.base_power**2.))

                self.assertEqual(load2.util_coeff_Q0,0)
                self.assertEqual(load2.util_coeff_Q1,300.*net.base_power)
                self.assertEqual(load2.util_coeff_Q2,-0.03*(net.base_power**2.))

    def test_sys_problem3(self):

        for case in test_cases.CASES:
            if case == '../data/sys_problem3.mat':

                net = pf.ParserMAT().parse(case)

                self.assertEqual(net.base_power,100.)

                # numbers
                self.assertEqual(net.num_buses,10)
                self.assertEqual(net.num_generators,5)
                self.assertEqual(net.num_loads,10)
                self.assertEqual(net.num_shunts,0)
                self.assertEqual(net.num_var_generators,0)
                self.assertEqual(net.num_branches,13)

                # buses
                bus1 = net.get_bus_by_number(1)
                bus2 = net.get_bus_by_number(2)
                bus3 = net.get_bus_by_number(3)
                bus4 = net.get_bus_by_number(4)
                bus5 = net.get_bus_by_number(5)
                bus6 = net.get_bus_by_number(6)
                bus7 = net.get_bus_by_number(7)
                bus8 = net.get_bus_by_number(8)
                bus9 = net.get_bus_by_number(9)
                bus10 = net.get_bus_by_number(10)

                # loads
                for bus in net.buses:
                    self.assertEqual(len(bus.loads),1)
                load1 = bus1.loads[0]
                load2 = bus2.loads[0]
                load3 = bus3.loads[0]
                load4 = bus4.loads[0]
                load5 = bus5.loads[0]
                load6 = bus6.loads[0]
                load7 = bus7.loads[0]
                load8 = bus8.loads[0]
                load9 = bus9.loads[0]
                load10 = bus10.loads[0]

                self.assertEqual(load1.bus,bus1)
                self.assertEqual(load2.bus,bus2)
                self.assertEqual(load3.bus,bus3)
                self.assertEqual(load4.bus,bus4)
                self.assertEqual(load5.bus,bus5)
                self.assertEqual(load6.bus,bus6)
                self.assertEqual(load7.bus,bus7)
                self.assertEqual(load8.bus,bus8)
                self.assertEqual(load9.bus,bus9)
                self.assertEqual(load10.bus,bus10)

                self.assertEqual(load1.P,55./100.)
                self.assertEqual(load2.P,55/100.)
                self.assertEqual(load3.P,1300/100.)
                self.assertEqual(load4.P,650/100.)
                self.assertEqual(load5.P,650/100.)
                self.assertEqual(load6.P,200/100.)
                self.assertEqual(load7.P,2600/100.)
                self.assertEqual(load8.P,3600/100.)
                self.assertEqual(load9.P,1100/100.)
                self.assertEqual(load10.P,1900/100.)
                for load in net.loads:
                    self.assertEqual(load.P_max,load.P)
                    self.assertEqual(load.P_min,load.P)

                # generators
                self.assertEqual(len(bus1.generators),0)
                self.assertEqual(len(bus2.generators),1)
                self.assertEqual(len(bus3.generators),1)
                self.assertEqual(len(bus4.generators),0)
                self.assertEqual(len(bus5.generators),1)
                self.assertEqual(len(bus6.generators),0)
                self.assertEqual(len(bus7.generators),1)
                self.assertEqual(len(bus8.generators),1)
                self.assertEqual(len(bus9.generators),0)
                self.assertEqual(len(bus10.generators),0)
                gen1 = bus2.generators[0]
                gen2 = bus3.generators[0]
                gen3 = bus5.generators[0]
                gen4 = bus7.generators[0]
                gen5 = bus8.generators[0]

                self.assertEqual(gen1.bus,bus2)
                self.assertEqual(gen2.bus,bus3)
                self.assertEqual(gen3.bus,bus5)
                self.assertEqual(gen4.bus,bus7)
                self.assertEqual(gen5.bus,bus8)

                self.assertEqual(gen1.P,50./100.)
                self.assertEqual(gen1.P_min,0)
                self.assertEqual(gen1.P_max,1200./100.)
                self.assertEqual(gen1.cost_coeff_Q0,0)
                self.assertEqual(gen1.cost_coeff_Q1,6.9*100)
                self.assertEqual(gen1.cost_coeff_Q2,0.00067*(100**2.))

                self.assertEqual(gen2.P,40./100.)
                self.assertEqual(gen2.P_min,0)
                self.assertEqual(gen2.P_max,8000./100.)
                self.assertEqual(gen2.cost_coeff_Q0,0)
                self.assertEqual(gen2.cost_coeff_Q1,24.3*100)
                self.assertEqual(gen2.cost_coeff_Q2,0.00040*(100**2.))

                self.assertEqual(gen3.P,30./100.)
                self.assertEqual(gen3.P_min,0)
                self.assertEqual(gen3.P_max,3000./100.)
                self.assertEqual(gen3.cost_coeff_Q0,0)
                self.assertEqual(gen3.cost_coeff_Q1,29.1*100)
                self.assertEqual(gen3.cost_coeff_Q2,0.00006*(100**2.))

                self.assertEqual(gen4.P,20./100.)
                self.assertEqual(gen4.P_min,0)
                self.assertEqual(gen4.P_max,800./100.)
                self.assertEqual(gen4.cost_coeff_Q0,0)
                self.assertEqual(gen4.cost_coeff_Q1,6.9*100)
                self.assertEqual(gen4.cost_coeff_Q2,0.00026*(100**2.))

                self.assertEqual(gen5.P,10./100.)
                self.assertEqual(gen5.P_min,0)
                self.assertEqual(gen5.P_max,2000./100.)
                self.assertEqual(gen5.cost_coeff_Q0,0)
                self.assertEqual(gen5.cost_coeff_Q1,50.*100)
                self.assertEqual(gen5.cost_coeff_Q2,0.0015*(100**2.))

                # branches
                branch1 = net.get_branch(12)
                branch2 = net.get_branch(11)
                branch3 = net.get_branch(10)
                branch4 = net.get_branch(9)
                branch5 = net.get_branch(8)
                branch6 = net.get_branch(7)
                branch7 = net.get_branch(6)
                branch8 = net.get_branch(5)
                branch9 = net.get_branch(4)
                branch10 = net.get_branch(3)
                branch11 = net.get_branch(2)
                branch12 = net.get_branch(1)
                branch13 = net.get_branch(0)

                self.assertEqual(branch1.bus_k,bus1)
                self.assertEqual(branch1.bus_m,bus3)
                self.assertLess(abs(branch1.b + 1./0.1),1e-10)
                self.assertEqual(branch1.ratingA,3000./100.)
                self.assertEqual(branch1.ratingB,3000./100.)
                self.assertEqual(branch1.ratingC,3000./100.)

                self.assertEqual(branch2.bus_k,bus1)
                self.assertEqual(branch2.bus_m,bus10)
                self.assertLess(abs(branch2.b + 1./0.27),1e-10)
                self.assertEqual(branch2.ratingA,2000./100.)
                self.assertEqual(branch2.ratingB,2000./100.)
                self.assertEqual(branch2.ratingC,2000./100.)

                self.assertEqual(branch3.bus_k,bus2)
                self.assertEqual(branch3.bus_m,bus3)
                self.assertLess(abs(branch3.b + 1./0.12),1e-10)
                self.assertEqual(branch3.ratingA,6500./100.)
                self.assertEqual(branch3.ratingB,6500./100.)
                self.assertEqual(branch3.ratingC,6500./100.)

                self.assertEqual(branch4.bus_k,bus2)
                self.assertEqual(branch4.bus_m,bus9)
                self.assertLess(abs(branch4.b + 1./0.07),1e-10)
                self.assertEqual(branch4.ratingA,5500./100.)
                self.assertEqual(branch4.ratingB,5500./100.)
                self.assertEqual(branch4.ratingC,5500./100.)

                self.assertEqual(branch5.bus_k,bus2)
                self.assertEqual(branch5.bus_m,bus10)
                self.assertLess(abs(branch5.b + 1./0.14),1e-10)
                self.assertEqual(branch5.ratingA,5500./100.)
                self.assertEqual(branch5.ratingB,5500./100.)
                self.assertEqual(branch5.ratingC,5500./100.)

                self.assertEqual(branch6.bus_k,bus3)
                self.assertEqual(branch6.bus_m,bus4)
                self.assertLess(abs(branch6.b + 1./0.1),1e-10)
                self.assertEqual(branch6.ratingA,3000./100.)
                self.assertEqual(branch6.ratingB,3000./100.)
                self.assertEqual(branch6.ratingC,3000./100.)

                self.assertEqual(branch7.bus_k,bus3)
                self.assertEqual(branch7.bus_m,bus5)
                self.assertLess(abs(branch7.b + 1./0.17),1e-10)
                self.assertEqual(branch7.ratingA,4000./100.)
                self.assertEqual(branch7.ratingB,4000./100.)
                self.assertEqual(branch7.ratingC,4000./100.)

                self.assertEqual(branch8.bus_k,bus4)
                self.assertEqual(branch8.bus_m,bus5)
                self.assertLess(abs(branch8.b + 1./0.17),1e-10)
                self.assertEqual(branch8.ratingA,4000./100.)
                self.assertEqual(branch8.ratingB,4000./100.)
                self.assertEqual(branch8.ratingC,4000./100.)

                self.assertEqual(branch9.bus_k,bus5)
                self.assertEqual(branch9.bus_m,bus6)
                self.assertLess(abs(branch9.b + 1./0.17),1e-10)
                self.assertEqual(branch9.ratingA,5000./100.)
                self.assertEqual(branch9.ratingB,5000./100.)
                self.assertEqual(branch9.ratingC,5000./100.)

                self.assertEqual(branch10.bus_k,bus6)
                self.assertEqual(branch10.bus_m,bus7)
                self.assertLess(abs(branch10.b + 1./0.16),1e-10)
                self.assertEqual(branch10.ratingA,2000./100.)
                self.assertEqual(branch10.ratingB,2000./100.)
                self.assertEqual(branch10.ratingC,2000./100.)

                self.assertEqual(branch11.bus_k,bus7)
                self.assertEqual(branch11.bus_m,bus8)
                self.assertLess(abs(branch11.b + 1./0.25),1e-10)
                self.assertEqual(branch11.ratingA,3000./100.)
                self.assertEqual(branch11.ratingB,3000./100.)
                self.assertEqual(branch11.ratingC,3000./100.)

                self.assertEqual(branch12.bus_k,bus8)
                self.assertEqual(branch12.bus_m,bus9)
                self.assertLess(abs(branch12.b + 1./0.25),1e-10)
                self.assertEqual(branch12.ratingA,2500./100.)
                self.assertEqual(branch12.ratingB,2500./100.)
                self.assertEqual(branch12.ratingC,2500./100.)

                self.assertEqual(branch13.bus_k,bus8)
                self.assertEqual(branch13.bus_m,bus10)
                self.assertLess(abs(branch13.b + 1./0.07),1e-10)
                self.assertEqual(branch13.ratingA,4000./100.)
                self.assertEqual(branch13.ratingB,4000./100.)
                self.assertEqual(branch13.ratingC,4000./100.)

    def test_cas32art(self):

        for case in test_cases.CASES:
            if case == '../data/case32.art':

                net = pf.ParserART().parse(case)

                self.assertEqual(net.num_buses,31)
                self.assertEqual(net.num_batteries,3)
                b1 = net.get_battery(0)
                b2 = net.get_battery(1)
                b3 = net.get_battery(2)
                self.assertRaises(pf.NetworkError,net.get_battery,3)

                self.assertEqual(b1.bus.name,"N18")
                self.assertEqual(b1.index,0)
                self.assertTrue(all(list(map(lambda y: isinstance(y,pf.Battery),b1.bus.batteries))))
                self.assertEqual(len(b1.bus.batteries),2)
                self.assertEqual(list(map(lambda y: y.index,b1.bus.batteries)),[b1.index,b2.index])
                self.assertEqual(b1.P,6./net.base_power)
                self.assertEqual(b1.P_min,-7./net.base_power)
                self.assertEqual(b1.P_max,8./net.base_power)
                self.assertEqual(b1.E,14./net.base_power)
                self.assertEqual(b1.E_max,22./net.base_power)
                self.assertEqual(b1.eta_c,0.93)
                self.assertEqual(b1.eta_d,0.97)

                self.assertEqual(b2.bus.name,"N18")
                self.assertEqual(b2.index,1)
                self.assertTrue(all(list(map(lambda y: isinstance(y,pf.Battery),b2.bus.batteries))))
                self.assertEqual(len(b2.bus.batteries),2)
                self.assertEqual(list(map(lambda y: y.index,b2.bus.batteries)),[b1.index,b2.index])
                self.assertEqual(b2.P,3./net.base_power)
                self.assertEqual(b2.P_min,-3./net.base_power)
                self.assertEqual(b2.P_max,9./net.base_power)
                self.assertEqual(b2.E,12./net.base_power)
                self.assertEqual(b2.E_max,21./net.base_power)
                self.assertEqual(b2.eta_c,0.94)
                self.assertEqual(b2.eta_d,0.92)

                self.assertEqual(b3.bus.name,"N15")
                self.assertEqual(b3.index,2)
                self.assertTrue(all(list(map(lambda y: isinstance(y,pf.Battery),b3.bus.batteries))))
                self.assertEqual(len(b3.bus.batteries),1)
                self.assertEqual(list(map(lambda y: y.index,b3.bus.batteries)),[b3.index])
                self.assertEqual(b3.P,2./net.base_power)
                self.assertEqual(b3.P_min,-4./net.base_power)
                self.assertEqual(b3.P_max,5./net.base_power)
                self.assertEqual(b3.E,10./net.base_power)
                self.assertEqual(b3.E_max,20./net.base_power)
                self.assertEqual(b3.eta_c,0.95)
                self.assertEqual(b3.eta_d,0.93)

    def test_ieee14_gen_cost(self):

        for case in test_cases.CASES:
            if case == '../data/ieee14.mat':

                net = pf.ParserMAT().parse(case)

                self.assertEqual(net.base_power,100.)
                self.assertEqual(net.num_buses,14)
                self.assertEqual(net.num_generators,5)

                gen0 = net.get_generator(net.num_generators-1)
                gen1 = net.get_generator(net.num_generators-2)

                self.assertEqual(gen0.P,232.4/100.)
                self.assertEqual(gen0.cost_coeff_Q2,(4.3029259900e-02)*(net.base_power**2.))
                self.assertEqual(gen0.cost_coeff_Q1,20.*net.base_power)

                self.assertEqual(gen1.P,40./100.)
                self.assertEqual(gen1.cost_coeff_Q2,0.25*(net.base_power**2.))
                self.assertEqual(gen1.cost_coeff_Q1,20.*net.base_power)
            
    def test_type_parsers(self):

        for case in test_cases.CASES:

            if case.split('.')[-1] == 'mat':
                self.assertRaises(pf.ParserError,pf.ParserART().parse,case)
                if pf.info['raw_parser']:
                    self.assertRaises(pf.ParserError,pf.ParserRAW().parse,case)
                net = pf.ParserMAT().parse(case)
                self.assertGreater(net.num_buses,0)
            elif case.split('.')[-1] == 'art':
                self.assertRaises(pf.ParserError,pf.ParserMAT().parse,case)
                if pf.info['raw_parser']:
                    self.assertRaises(pf.ParserError,pf.ParserRAW().parse,case)
                net = pf.ParserART().parse(case)
                self.assertGreater(net.num_buses,0)
            elif case.split('.')[-1] == 'raw':
                self.assertRaises(pf.ParserError,pf.ParserMAT().parse,case)
                self.assertRaises(pf.ParserError,pf.ParserART().parse,case)
                if pf.info['raw_parser']:
                    net = pf.ParserRAW().parse(case)
                    self.assertGreater(net.num_buses,0)

    def test_json_parser(self):

        import os
        from numpy.linalg import norm
        eps = 1e-10

        norminf = lambda x: norm(x,np.inf) if not np.isscalar(x) else np.abs(x)

        for case in test_cases.CASES:

            T = 4
                
            net = pf.Parser(case).parse(case,T)
            self.assertEqual(net.num_periods,T)
            
            # Set flags
            net.set_flags('bus','variable','any','voltage magnitude')
            self.assertEqual(net.num_vars,net.num_buses*T)
            
            # Add vargens and betteries
            net.add_var_generators(net.get_load_buses(),100.,50.,30.,5,0.05)
            net.add_batteries(net.get_generator_buses(),20.,50.)
            
            # Some perturbations to reduce luck
            for bus in net.buses:
                bus.price = np.random.randn(T)

            try:
                
                json_parser = pf.ParserJSON()
                
                json_parser.write(net,"temp_json.json")
                
                new_net = json_parser.parse("temp_json.json")
                self.assertEqual(new_net.num_periods,T)

                # Network
                self.assertTrue(net is not new_net)
                self.assertEqual(net.num_periods,new_net.num_periods)
                self.assertEqual(net.base_power,new_net.base_power)
                self.assertNotEqual(net.num_vars,new_net.num_vars)
                self.assertEqual(new_net.num_vars,0)
                
                # Buses
                self.assertEqual(net.num_buses,new_net.num_buses)
                for i in range(net.num_buses):
                    bus = net.buses[i]
                    new_bus = new_net.buses[i]
                    self.assertEqual(bus.num_periods,T)
                    self.assertTrue(bus is not new_bus)
                    self.assertEqual(bus.number,new_bus.number)
                    self.assertEqual(bus.name,new_bus.name)
                    self.assertLess(norminf(bus.v_mag-new_bus.v_mag),eps)
                    self.assertLess(norminf(bus.v_ang-new_bus.v_ang),eps)
                    self.assertLess(norminf(bus.v_set-new_bus.v_set),eps)
                    self.assertLess(norminf(bus.v_mag-new_bus.v_mag),eps)
                    self.assertLess(norminf(bus.v_max_reg-new_bus.v_max_reg),eps)
                    self.assertLess(norminf(bus.v_min_reg-new_bus.v_min_reg),eps)
                    self.assertLess(norminf(bus.v_max_norm-new_bus.v_max_norm),eps)
                    self.assertLess(norminf(bus.v_min_norm-new_bus.v_min_norm),eps)
                    self.assertLess(norminf(bus.v_max_emer-new_bus.v_max_emer),eps)
                    self.assertLess(norminf(bus.v_min_emer-new_bus.v_min_emer),eps)
                    self.assertEqual(bus.is_slack(),new_bus.is_slack())
                    self.assertEqual(bus.is_regulated_by_gen(),new_bus.is_regulated_by_gen())
                    self.assertEqual(bus.is_regulated_by_tran(),new_bus.is_regulated_by_tran())
                    self.assertEqual(bus.is_regulated_by_shunt(),new_bus.is_regulated_by_shunt())
                    self.assertLess(norminf(bus.price-new_bus.price),eps)
                    self.assertEqual(set([o.index for o in bus.generators]),
                                     set([o.index for o in new_bus.generators]))
                    self.assertEqual(set([o.index for o in bus.reg_generators]),
                                     set([o.index for o in new_bus.reg_generators]))
                    self.assertEqual(set([o.index for o in bus.loads]),
                                     set([o.index for o in new_bus.loads]))
                    self.assertEqual(set([o.index for o in bus.shunts]),
                                     set([o.index for o in new_bus.shunts]))
                    self.assertEqual(set([o.index for o in bus.reg_shunts]),
                                     set([o.index for o in new_bus.reg_shunts]))
                    self.assertEqual(set([o.index for o in bus.branches_k]),
                                     set([o.index for o in new_bus.branches_k]))
                    self.assertEqual(set([o.index for o in bus.branches_m]),
                                     set([o.index for o in new_bus.branches_m]))
                    self.assertEqual(set([o.index for o in bus.reg_trans]),
                                     set([o.index for o in new_bus.reg_trans]))
                    self.assertEqual(set([o.index for o in bus.var_generators]),
                                     set([o.index for o in new_bus.var_generators]))
                    self.assertEqual(set([o.index for o in bus.batteries]),
                                     set([o.index for o in new_bus.batteries]))

                # Branches
                self.assertEqual(net.num_branches,new_net.num_branches)
                for i in range(net.num_branches):
                    branch = net.branches[i]
                    new_branch = new_net.branches[i]
                    self.assertTrue(branch is not new_branch)
                    self.assertEqual(branch.num_periods,T)
                    self.assertEqual(branch.num_periods,new_branch.num_periods)
                    self.assertEqual(branch.bus_k.index,new_branch.bus_k.index)
                    self.assertEqual(branch.bus_m.index,new_branch.bus_m.index)
                    self.assertEqual(branch.is_fixed_tran(),new_branch.is_fixed_tran())
                    self.assertEqual(branch.is_line(),new_branch.is_line())
                    self.assertEqual(branch.is_phase_shifter(),new_branch.is_phase_shifter())
                    self.assertEqual(branch.is_tap_changer(),new_branch.is_tap_changer())
                    self.assertEqual(branch.is_tap_changer_v(),new_branch.is_tap_changer_v())
                    self.assertEqual(branch.is_tap_changer_Q(),new_branch.is_tap_changer_Q())
                    if branch.is_tap_changer_v():
                        self.assertEqual(branch.reg_bus.index,new_branch.reg_bus.index)
                    self.assertLess(norminf(branch.g-new_branch.g),eps*(1+norminf(branch.g)))
                    self.assertLess(norminf(branch.g_k-new_branch.g_k),eps*(1+norminf(branch.g_k)))
                    self.assertLess(norminf(branch.g_m-new_branch.g_m),eps*(1+norminf(branch.g_m)))
                    self.assertLess(norminf(branch.b-new_branch.b),eps*(1+norminf(branch.b)))
                    self.assertLess(norminf(branch.b_k-new_branch.b_k),eps*(1+norminf(branch.b_k)))
                    self.assertLess(norminf(branch.b_m-new_branch.b_m),eps*(1+norminf(branch.b_m)))
                    self.assertLess(norminf(branch.ratio-new_branch.ratio),eps)
                    self.assertLess(norminf(branch.ratio_max-new_branch.ratio_max),eps)
                    self.assertLess(norminf(branch.ratio_min-new_branch.ratio_min),eps)
                    self.assertLess(norminf(branch.phase-new_branch.phase),eps)
                    self.assertLess(norminf(branch.phase_max-new_branch.phase_max),eps)
                    self.assertLess(norminf(branch.phase_min-new_branch.phase_min),eps)
                    self.assertLess(norminf(branch.ratingA-new_branch.ratingA),eps)
                    self.assertLess(norminf(branch.ratingB-new_branch.ratingB),eps)
                    self.assertLess(norminf(branch.ratingC-new_branch.ratingC),eps)
                    self.assertEqual(branch.is_on_outage(),new_branch.is_on_outage())
                    self.assertEqual(branch.has_pos_ratio_v_sens(),new_branch.has_pos_ratio_v_sens())

                # Generators
                self.assertEqual(net.num_generators,new_net.num_generators)
                for i in range(net.num_generators):
                    gen = net.generators[i]
                    new_gen = new_net.generators[i]
                    self.assertTrue(gen is not new_gen)
                    self.assertEqual(gen.num_periods,T)
                    self.assertEqual(gen.num_periods,new_gen.num_periods)
                    self.assertEqual(gen.bus.index,new_gen.bus.index)
                    self.assertEqual(gen.is_on_outage(),new_gen.is_on_outage())
                    self.assertEqual(gen.is_slack(),new_gen.is_slack())
                    self.assertEqual(gen.is_regulator(),new_gen.is_regulator())
                    self.assertEqual(gen.is_P_adjustable(),new_gen.is_P_adjustable())
                    if gen.is_regulator():
                        self.assertEqual(gen.reg_bus.index,new_gen.reg_bus.index)
                    self.assertLess(norminf(gen.P-new_gen.P),eps)
                    self.assertLess(norminf(gen.P_max-new_gen.P_max),eps)
                    self.assertLess(norminf(gen.P_min-new_gen.P_min),eps)
                    self.assertLess(norminf(gen.dP_max-new_gen.dP_max),eps)
                    self.assertLess(norminf(gen.P_prev-new_gen.P_prev),eps)
                    self.assertLess(norminf(gen.Q-new_gen.Q),eps)
                    self.assertLess(norminf(gen.Q_max-new_gen.Q_max),eps)
                    self.assertLess(norminf(gen.Q_min-new_gen.Q_min),eps)
                    self.assertLess(norminf(gen.cost_coeff_Q0-new_gen.cost_coeff_Q0),eps)
                    self.assertLess(norminf(gen.cost_coeff_Q1-new_gen.cost_coeff_Q1),eps)
                    self.assertLess(norminf(gen.cost_coeff_Q2-new_gen.cost_coeff_Q2),eps)

                # Var generators
                self.assertEqual(net.num_var_generators,new_net.num_var_generators)
                for i in range(net.num_var_generators):
                    vargen = net.var_generators[i]
                    new_vargen = new_net.var_generators[i]
                    self.assertTrue(vargen is not new_vargen)
                    self.assertEqual(vargen.num_periods,T)
                    self.assertEqual(vargen.num_periods,new_vargen.num_periods)
                    self.assertEqual(vargen.bus.index,new_vargen.bus.index)
                    self.assertEqual(vargen.name,new_vargen.name)
                    self.assertLess(norminf(vargen.P-new_vargen.P),eps)
                    self.assertLess(norminf(vargen.P_ava-new_vargen.P_ava),eps)
                    self.assertLess(norminf(vargen.P_max-new_vargen.P_max),eps)
                    self.assertLess(norminf(vargen.P_min-new_vargen.P_min),eps)
                    self.assertLess(norminf(vargen.P_std-new_vargen.P_std),eps)
                    self.assertLess(norminf(vargen.Q-new_vargen.Q),eps)
                    self.assertLess(norminf(vargen.Q_max-new_vargen.Q_max),eps)
                    self.assertLess(norminf(vargen.Q_min-new_vargen.Q_min),eps)

                # Shunts
                self.assertEqual(net.num_shunts,new_net.num_shunts)
                for i in range(net.num_shunts):
                    shunt = net.shunts[i]
                    new_shunt = new_net.shunts[i]
                    self.assertTrue(shunt is not new_shunt)
                    self.assertEqual(shunt.num_periods,T)
                    self.assertEqual(shunt.num_periods,new_shunt.num_periods)
                    self.assertEqual(shunt.bus.index,new_shunt.bus.index)
                    self.assertEqual(shunt.is_fixed(),new_shunt.is_fixed())
                    self.assertEqual(shunt.is_switched_v(),new_shunt.is_switched_v())
                    if shunt.is_switched_v():
                        self.assertEqual(shunt.reg_bus.index,new_shunt.reg_bus.index)
                    self.assertLess(norminf(shunt.g-new_shunt.g),eps*(1+norminf(shunt.g)))
                    self.assertLess(norminf(shunt.b-new_shunt.b),eps*(1+norminf(shunt.b)))
                    self.assertLess(norminf(shunt.b_max-new_shunt.b_max),eps*(1+norminf(shunt.b_max)))
                    self.assertLess(norminf(shunt.b_min-new_shunt.b_min),eps*(1+norminf(shunt.b_min)))
                    
                # Loads
                self.assertEqual(net.num_loads,new_net.num_loads)
                for i in range(net.num_loads):
                    load = net.loads[i]
                    new_load = new_net.loads[i]
                    self.assertTrue(load is not new_load)
                    self.assertEqual(load.num_periods,T)
                    self.assertEqual(load.num_periods,new_load.num_periods)
                    self.assertEqual(load.bus.index,new_load.bus.index)
                    self.assertLess(norminf(load.P-new_load.P),eps)
                    self.assertLess(norminf(load.P_max-new_load.P_max),eps)
                    self.assertLess(norminf(load.P_min-new_load.P_min),eps)
                    self.assertLess(norminf(load.Q-new_load.Q),eps)
                    self.assertLess(norminf(load.target_power_factor-new_load.target_power_factor),eps)
                    self.assertLess(norminf(load.util_coeff_Q0-new_load.util_coeff_Q0),eps)
                    self.assertLess(norminf(load.util_coeff_Q1-new_load.util_coeff_Q1),eps)
                    self.assertLess(norminf(load.util_coeff_Q2-new_load.util_coeff_Q2),eps)                    

                # Batteries
                self.assertEqual(net.num_batteries,new_net.num_batteries)
                for i in range(net.num_batteries):
                    bat = net.batteries[i]
                    new_bat = new_net.batteries[i]
                    self.assertTrue(bat is not new_bat)
                    self.assertEqual(bat.num_periods,T)
                    self.assertEqual(bat.num_periods,new_bat.num_periods)
                    self.assertEqual(bat.bus.index,new_bat.bus.index)
                    self.assertLess(norminf(bat.P-new_bat.P),eps)
                    self.assertLess(norminf(bat.P_max-new_bat.P_max),eps)
                    self.assertLess(norminf(bat.P_min-new_bat.P_min),eps)
                    self.assertLess(norminf(bat.eta_c-new_bat.eta_c),eps)
                    self.assertLess(norminf(bat.eta_d-new_bat.eta_d),eps)
                    self.assertLess(norminf(bat.E-new_bat.E),eps)
                    self.assertLess(norminf(bat.E_init-new_bat.E_init),eps)
                    self.assertLess(norminf(bat.E_final-new_bat.E_final),eps)
                    self.assertLess(norminf(bat.E_max-new_bat.E_max),eps)
                    
                # Hashes                 
                for bus in net.buses:
                    self.assertEqual(bus.index,net.get_bus_by_number(bus.number).index)
                    self.assertEqual(bus.name,net.get_bus_by_name(bus.name).name)
                for vargen in net.var_generators:
                    self.assertEqual(vargen.index,net.get_var_generator_by_name(vargen.name).index)
                    
            finally:
                
                os.remove("temp_json.json")
