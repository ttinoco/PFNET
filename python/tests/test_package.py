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

class TestPackage(unittest.TestCase):
    
    def setUp(self):
        
        pass

    def test_info(self):

        self.assertTrue(isinstance(pf.info,dict))
        self.assertTrue('graphviz' in pf.info)
        self.assertTrue('raw parser' in pf.info)
        
        
