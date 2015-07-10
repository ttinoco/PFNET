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

class TestGraph(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_construction(self):

        net = self.net

        for case in test_cases.CASES:
            
            net.load(case)
            
            g = pf.Graph(net)
            
            # Missing layout
            self.assertRaises(pf.GraphError,g.write,"png","foo")

            # Bad sens type
            self.assertRaises(pf.GraphError,g.color_nodes_by_sensitivity,pf.BUS_MIS_ACTIVE)

            # Bad mis type
            self.assertRaises(pf.GraphError,g.color_nodes_by_mismatch,pf.BUS_SENS_Q_BALANCE)
            
    def tearDown(self):
        
        pass
                
                
