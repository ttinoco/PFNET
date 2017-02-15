#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cheur

# Types
HEUR_TYPE_PVPQ = cheur.HEUR_TYPE_PVPQ

class HeuristicError(Exception):
    """
    Heuristic error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Heuristic:
    """
    Heuristic class.
    """

    pass
