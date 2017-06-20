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

class TestGraph(unittest.TestCase):
    
    def setUp(self):
        
        pass

    def test_construction(self):

        for case in test_cases.CASES:
            
            net = pf.Parser(case).parse(case)
            self.assertEqual(net.num_periods,1)
            
            g = pf.Graph(net)
            
            # Missing layout
            self.assertRaises(pf.GraphError,g.write,"png","foo")
            self.assertTrue(g.has_error())
            g.clear_error()
            self.assertFalse(g.has_error())
            
            # Bad sens type
            self.assertRaises(pf.GraphError,g.color_nodes_by_sensitivity,pf.BUS_MIS_ACTIVE)
            self.assertTrue(g.has_error())
            g.clear_error()
            self.assertFalse(g.has_error())

            # Bad mis type
            self.assertRaises(pf.GraphError,g.color_nodes_by_mismatch,pf.BUS_SENS_Q_BALANCE)
            self.assertTrue(g.has_error())
            g.clear_error()
            self.assertFalse(g.has_error())

            g.set_layout()

            if g.has_viz():
                g.write("dot","foo.dot")
            else:
                self.assertRaises(pf.GraphError,g.write,"dot","foo.dot")         
                
    def tearDown(self):
        
        pass
                
                
