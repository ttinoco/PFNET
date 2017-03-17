#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

class BaseParser:
    """
    Parser base class.
    """

    def __init__(self):
        """
        Constructor.
        """

        #: Output level (int)
        self.output_level = 0

        #: Error flag (boolean)
        self.error_flag = False

        #: Error string (string)
        self.error_string = ""

    def read(self,filename):
        """
        Reads file and stores data.

        Parameters
        ----------
        filename : string
        """
        
        pass

    def set(self,key,value):
        """
        Sets parser option.

        Parameters
        ----------
        key : string
        value : float
        """

        pass

    def show(self):
        """
        Shows information about parsed data.
        """

        pass

    def load(self,network):
        """
        Loads network using parsed data.

        Parameters
        ----------
        network : Network
        """
        
        pass

    def has_error(self):
        """
        Checks whether parser encountered 
        any error.

        Returns
        -------
        flag : {True,False}
        """

        return True

    def get_error_string(self):
        """
        Gives error message.
        
        Returns
        -------
        error : string
        """

        return "This is just a dummy python parser"
        
    
