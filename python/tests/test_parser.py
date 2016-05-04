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

class TestParser(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_sys_problem2(self):

        for case in test_cases.CASES:
            if case == '../data/sys_problem2.mat':

                net = self.net

                net.load(case)

                self.assertEqual(net.base_power,100.)
                
                self.assertEqual(net.num_buses,3)
                self.assertEqual(net.num_gens,4)
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
                
                self.assertEqual(len(bus1.gens),2)
                self.assertEqual(len(bus2.gens),1)
                self.assertEqual(len(bus3.gens),1)
                
                self.assertEqual(len(bus1.loads),1)
                self.assertEqual(len(bus2.loads),1)
                self.assertEqual(len(bus3.loads),1)
                
                gen1a = bus1.gens[0]
                gen1b = bus1.gens[1]
                gen2 = bus2.gens[0]
                gen3 = bus3.gens[0]
                
                load1 = bus1.loads[0]
                load2 = bus2.loads[0]
                load3 = bus3.loads[0]
                
                self.assertEqual(load1.bus,bus1)
                self.assertEqual(load2.bus,bus2)
                self.assertEqual(load3.bus,bus3)
                
                self.assertEqual(gen1a.P_max,50)
                self.assertEqual(gen1b.P_max,50)
                self.assertEqual(gen2.P_max,25)
                self.assertEqual(gen3.P_max,19)
                self.assertEqual(gen1a.bus,bus1)
                self.assertEqual(gen1b.bus,bus1)
                self.assertEqual(gen2.bus,bus2)
                self.assertEqual(gen3.bus,bus3)
                for gen in net.generators:
                    self.assertEqual(gen.P_min,0.)
                    
                self.assertEqual(branch13.bus_from.number,bus1.number)
                self.assertEqual(branch13.bus_to.number,bus3.number)

                self.assertEqual(branch23.bus_from.number,bus2.number)
                self.assertEqual(branch23.bus_to.number,bus3.number)
                
                self.assertEqual(branch12.bus_from.number,bus1.number)
                self.assertEqual(branch12.bus_to.number,bus2.number)
                
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
                
                self.assertEqual(gen1a.cost_coeff_Q0,0)
                self.assertEqual(gen1a.cost_coeff_Q1,5.*net.base_power)
                self.assertEqual(gen1a.cost_coeff_Q2,0.02*(net.base_power**2.))
                
                self.assertEqual(gen1b.cost_coeff_Q0,0)
                self.assertEqual(gen1b.cost_coeff_Q1,6.*net.base_power)
                self.assertEqual(gen1b.cost_coeff_Q2,0.03*(net.base_power**2.))
                
                self.assertEqual(gen2.cost_coeff_Q0,0)
                self.assertEqual(gen2.cost_coeff_Q1,12.*net.base_power)
                self.assertEqual(gen2.cost_coeff_Q2,0.06*(net.base_power**2.))
                
                self.assertEqual(gen3.cost_coeff_Q0,0)
                self.assertEqual(gen3.cost_coeff_Q1,10.*net.base_power)
                self.assertEqual(gen3.cost_coeff_Q2,0.08*(net.base_power**2.))
                
                # Load utility
                self.assertEqual(load1.util_coeff_Q0,0)
                self.assertEqual(load1.util_coeff_Q1,400.*net.base_power)
                self.assertEqual(load1.util_coeff_Q2,-0.03*(net.base_power**2.))
                
                self.assertEqual(load2.util_coeff_Q0,0)
                self.assertEqual(load2.util_coeff_Q1,450.*net.base_power)
                self.assertEqual(load2.util_coeff_Q2,-0.02*(net.base_power**2.))
                
                self.assertEqual(load3.util_coeff_Q0,0)
                self.assertEqual(load3.util_coeff_Q1,300.*net.base_power)
                self.assertEqual(load3.util_coeff_Q2,-0.03*(net.base_power**2.))
                
    def test_sys_problem3(self):

        for case in test_cases.CASES:            
            if case == '../data/sys_problem3.mat':

                net = self.net

                net.load(case)

                self.assertEqual(net.base_power,100.)

                # numbers
                self.assertEqual(net.num_buses,10)
                self.assertEqual(net.num_gens,5)
                self.assertEqual(net.num_loads,10)
                self.assertEqual(net.num_shunts,0)
                self.assertEqual(net.num_vargens,0)
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
                self.assertEqual(len(bus1.gens),0)
                self.assertEqual(len(bus2.gens),1)
                self.assertEqual(len(bus3.gens),1)
                self.assertEqual(len(bus4.gens),0)
                self.assertEqual(len(bus5.gens),1)
                self.assertEqual(len(bus6.gens),0)
                self.assertEqual(len(bus7.gens),1)
                self.assertEqual(len(bus8.gens),1)
                self.assertEqual(len(bus9.gens),0)
                self.assertEqual(len(bus10.gens),0)
                gen1 = bus2.gens[0]
                gen2 = bus3.gens[0]
                gen3 = bus5.gens[0]
                gen4 = bus7.gens[0]
                gen5 = bus8.gens[0]
                
                self.assertEqual(gen1.bus,bus2)
                self.assertEqual(gen2.bus,bus3)
                self.assertEqual(gen3.bus,bus5)
                self.assertEqual(gen4.bus,bus7)
                self.assertEqual(gen5.bus,bus8)
                
                self.assertEqual(gen1.P_min,0)
                self.assertEqual(gen1.P_max,1200./100.)
                self.assertEqual(gen1.cost_coeff_Q0,0)
                self.assertEqual(gen1.cost_coeff_Q1,6.9*100)
                self.assertEqual(gen1.cost_coeff_Q2,0.00067*(100**2.))    
                
                self.assertEqual(gen2.P_min,0)
                self.assertEqual(gen2.P_max,8000./100.)
                self.assertEqual(gen2.cost_coeff_Q0,0)
                self.assertEqual(gen2.cost_coeff_Q1,24.3*100)
                self.assertEqual(gen2.cost_coeff_Q2,0.00040*(100**2.))    
                
                self.assertEqual(gen3.P_min,0)
                self.assertEqual(gen3.P_max,3000./100.)
                self.assertEqual(gen3.cost_coeff_Q0,0)
                self.assertEqual(gen3.cost_coeff_Q1,29.1*100)
                self.assertEqual(gen3.cost_coeff_Q2,0.00006*(100**2.))    
                
                self.assertEqual(gen4.P_min,0)
                self.assertEqual(gen4.P_max,800./100.)
                self.assertEqual(gen4.cost_coeff_Q0,0)
                self.assertEqual(gen4.cost_coeff_Q1,6.9*100)
                self.assertEqual(gen4.cost_coeff_Q2,0.00026*(100**2.))    
                
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
                
                self.assertEqual(branch1.bus_from,bus1)
                self.assertEqual(branch1.bus_to,bus3)
                self.assertLess(abs(branch1.b + 1./0.1),1e-10)
                self.assertEqual(branch1.ratingA,3000./100.)
                self.assertEqual(branch1.ratingB,3000./100.)
                self.assertEqual(branch1.ratingC,3000./100.)
            
                self.assertEqual(branch2.bus_from,bus1)
                self.assertEqual(branch2.bus_to,bus10)
                self.assertLess(abs(branch2.b + 1./0.27),1e-10)
                self.assertEqual(branch2.ratingA,2000./100.)
                self.assertEqual(branch2.ratingB,2000./100.)
                self.assertEqual(branch2.ratingC,2000./100.)
                
                self.assertEqual(branch3.bus_from,bus2)
                self.assertEqual(branch3.bus_to,bus3)
                self.assertLess(abs(branch3.b + 1./0.12),1e-10)
                self.assertEqual(branch3.ratingA,6500./100.)
                self.assertEqual(branch3.ratingB,6500./100.)
                self.assertEqual(branch3.ratingC,6500./100.)
                
                self.assertEqual(branch4.bus_from,bus2)
                self.assertEqual(branch4.bus_to,bus9)
                self.assertLess(abs(branch4.b + 1./0.07),1e-10)
                self.assertEqual(branch4.ratingA,5500./100.)
                self.assertEqual(branch4.ratingB,5500./100.)
                self.assertEqual(branch4.ratingC,5500./100.)
                
                self.assertEqual(branch5.bus_from,bus2)
                self.assertEqual(branch5.bus_to,bus10)
                self.assertLess(abs(branch5.b + 1./0.14),1e-10)
                self.assertEqual(branch5.ratingA,5500./100.)
                self.assertEqual(branch5.ratingB,5500./100.)
                self.assertEqual(branch5.ratingC,5500./100.)
                
                self.assertEqual(branch6.bus_from,bus3)
                self.assertEqual(branch6.bus_to,bus4)
                self.assertLess(abs(branch6.b + 1./0.1),1e-10)
                self.assertEqual(branch6.ratingA,3000./100.)
                self.assertEqual(branch6.ratingB,3000./100.)
                self.assertEqual(branch6.ratingC,3000./100.)
                
                self.assertEqual(branch7.bus_from,bus3)
                self.assertEqual(branch7.bus_to,bus5)
                self.assertLess(abs(branch7.b + 1./0.17),1e-10)
                self.assertEqual(branch7.ratingA,4000./100.)
                self.assertEqual(branch7.ratingB,4000./100.)
                self.assertEqual(branch7.ratingC,4000./100.)
                
                self.assertEqual(branch8.bus_from,bus4)
                self.assertEqual(branch8.bus_to,bus5)
                self.assertLess(abs(branch8.b + 1./0.17),1e-10)
                self.assertEqual(branch8.ratingA,4000./100.)
                self.assertEqual(branch8.ratingB,4000./100.)
                self.assertEqual(branch8.ratingC,4000./100.)
                
                self.assertEqual(branch9.bus_from,bus5)
                self.assertEqual(branch9.bus_to,bus6)
                self.assertLess(abs(branch9.b + 1./0.17),1e-10)
                self.assertEqual(branch9.ratingA,5000./100.)
                self.assertEqual(branch9.ratingB,5000./100.)
                self.assertEqual(branch9.ratingC,5000./100.)
                
                self.assertEqual(branch10.bus_from,bus6)
                self.assertEqual(branch10.bus_to,bus7)
                self.assertLess(abs(branch10.b + 1./0.16),1e-10)
                self.assertEqual(branch10.ratingA,2000./100.)
                self.assertEqual(branch10.ratingB,2000./100.)
                self.assertEqual(branch10.ratingC,2000./100.)
                
                self.assertEqual(branch11.bus_from,bus7)
                self.assertEqual(branch11.bus_to,bus8)
                self.assertLess(abs(branch11.b + 1./0.25),1e-10)
                self.assertEqual(branch11.ratingA,3000./100.)
                self.assertEqual(branch11.ratingB,3000./100.)
                self.assertEqual(branch11.ratingC,3000./100.)
                
                self.assertEqual(branch12.bus_from,bus8)
                self.assertEqual(branch12.bus_to,bus9)
                self.assertLess(abs(branch12.b + 1./0.25),1e-10)
                self.assertEqual(branch12.ratingA,2500./100.)
                self.assertEqual(branch12.ratingB,2500./100.)
                self.assertEqual(branch12.ratingC,2500./100.)
                
                self.assertEqual(branch13.bus_from,bus8)
                self.assertEqual(branch13.bus_to,bus10)
                self.assertLess(abs(branch13.b + 1./0.07),1e-10)
                self.assertEqual(branch13.ratingA,4000./100.)
                self.assertEqual(branch13.ratingB,4000./100.)
                self.assertEqual(branch13.ratingC,4000./100.)

    def test_cas32art(self):

        for case in test_cases.CASES:
            if case == '../data/case32.art':

                net = self.net

                net.load(case)

                self.assertEqual(net.num_buses,31)
                self.assertEqual(net.num_bats,3)
                b1 = net.get_bat(0)
                b2 = net.get_bat(1)
                b3 = net.get_bat(2)
                self.assertRaises(pf.NetworkError,net.get_bat,3)

                self.assertEqual(b1.bus.name,"N18")
                self.assertEqual(b1.index,0)
                self.assertTrue(all(map(lambda y: isinstance(y,pf.Battery),b1.bus.bats)))
                self.assertEqual(len(b1.bus.bats),2)
                self.assertEqual(map(lambda y: y.index,b1.bus.bats),[b1.index,b2.index])
                self.assertEqual(b1.P,6./100.)
                self.assertEqual(b1.P_min,-7./100.)
                self.assertEqual(b1.P_max,8./100.)
                self.assertEqual(b1.E,14./100.)
                self.assertEqual(b1.E_max,22./100.)
                self.assertEqual(b1.eta_c,0.93)
                self.assertEqual(b1.eta_d,0.97)

                self.assertEqual(b2.bus.name,"N18")
                self.assertEqual(b2.index,1)
                self.assertTrue(all(map(lambda y: isinstance(y,pf.Battery),b2.bus.bats)))
                self.assertEqual(len(b2.bus.bats),2)
                self.assertEqual(map(lambda y: y.index,b2.bus.bats),[b1.index,b2.index])
                self.assertEqual(b2.P,3./100.)
                self.assertEqual(b2.P_min,-3./100.)
                self.assertEqual(b2.P_max,9./100.)
                self.assertEqual(b2.E,12./100.)
                self.assertEqual(b2.E_max,21./100.)
                self.assertEqual(b2.eta_c,0.94)
                self.assertEqual(b2.eta_d,0.92)

                self.assertEqual(b3.bus.name,"N15")
                self.assertEqual(b3.index,2)
                self.assertTrue(all(map(lambda y: isinstance(y,pf.Battery),b3.bus.bats)))
                self.assertEqual(len(b3.bus.bats),1)
                self.assertEqual(map(lambda y: y.index,b3.bus.bats),[b3.index])
                self.assertEqual(b3.P,2./100.)
                self.assertEqual(b3.P_min,-4./100.)
                self.assertEqual(b3.P_max,5./100.)
                self.assertEqual(b3.E,10./100.)
                self.assertEqual(b3.E_max,20./100.)
                self.assertEqual(b3.eta_c,0.95)
                self.assertEqual(b3.eta_d,0.93)
