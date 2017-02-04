#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

from __future__ import print_function
from parser import BaseParser

class DummyParser(BaseParser):

    def __init__(self):
        
        BaseParser.__init__(self)
        
    def __del__(self):

        pass

    def read(self,filename):
        
        pass
        
    def set(self,key,value):

        if key == 'output_level':
            self.output_level = value

    def show(self):

        pass

    def load(self,network):
        
        pass
        
    def has_error(self):

        return self.error_flag

    def get_error_string(self):

        return self.error_string

    
