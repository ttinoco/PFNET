#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

from parser import PythonParser

class DummyPythonParser(PythonParser):

    def __init__(self):

        print "DummyPythonParser created"
        
        self.error_flag = False
        self.error_string = ""

    def __del__(self):

        print "DummyPythonParser destroyed"

    def read(self,filename):
        
        print "DummyPythonParser asked to read file %s" %filename

    def show(self):

        print "DummyPythonParser asked to show"

    def load(self,network):
        
        print "DummyPythonParser asked to load network"
        print "Network received by DummyPythonParser is of type",type(network)
        network.show_components()

    def has_error(self):

        return self.error_flag

    def get_error_string(self):

        return self.error_string

    
