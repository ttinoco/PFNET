#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import unittest
import numpy as np
import pfnet as pf
from . import test_cases

class TestOutages(unittest.TestCase):

    def test_generators(self):

        # set/get outage

        # json

        # copy

        pass

    def test_branches(self):

        # set/get outage

        # json

        # copy

        pass

    def test_buses(self):

        # total gen P

        # total gen Q

        # total gen Qmin
        
        # total gen Qmax

        # total reg gen Qmin

        # total reg gen Qmax

        # reg by gen

        # reg by tran
        
        pass
    
    def test_network(self):

        # clear outages

        # num branches not on outage

        # num gens not on outage

        pass

    def test_network_properties(self):

        # tap ratio vio and actions

        # phase shift vio and actions

        # mismatches (branch and gen outages)

        # v reg limit violations

        # v set deviations and actions

        # gen active power cost

        # gen Q vio

        # gen P vio and actions

        pass

    def test_functions(self):

        # gen cost

        # load util

        # netcon cost

        # reg phase

        # reg pq

        # reg ratio

        # reg susc

        # reg vang

        # reg var

        # reg vmag

        # reg slim vmag

        # sp controls

        pass

    def test_constraints(self):

        pass

    
