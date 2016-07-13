#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

class PythonParser:
    """
    Python parser base class.
    """

    def read(self,filename):
        """
        Reads file and stores data.

        Parameters
        ----------
        filename : string
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
        
    
