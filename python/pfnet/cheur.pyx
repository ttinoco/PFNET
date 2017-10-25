#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
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

    pass

cdef class Heuristic:
    """
    Heuristic class.
    """

    pass
