#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

from pfnet import CustomFunction

class DummyGenCost(CustomFunction):

    def count_step(self,branch,t):

        print 'count step'
        
    def allocate(self):

        print 'allocate'

    def clear(self):

        print 'clear'
        
    def analyze_step(self,branch,t):

        print 'analzye step'

    def eval_step(self,branch,t):

        print 'eval step'
