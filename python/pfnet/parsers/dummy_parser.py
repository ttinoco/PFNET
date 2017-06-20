#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

from __future__ import print_function
from pfnet import CustomParser, Network, ParserError

class DummyParser(CustomParser):

    def init(self):

        self.some_init_data = 3

    def parse(self,filename,num_periods=1):
        
        if filename.split('.')[-1] != 'dummy':
            raise ParserError('invalid extension')
        return Network(num_periods)

    def set(self,key,value):

        pass

    def show(self):

        print('this is a dummy parser')

    def write(self,net,filename):
        
        pass
        
    
