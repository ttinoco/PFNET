#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/load.h":

    ctypedef struct Load
    ctypedef struct Bus
    ctypedef double REAL
        
    char LOAD_get_obj_type(Load* load)
    int LOAD_get_index(Load* load)
    Bus* LOAD_get_bus(Load* load)
    REAL LOAD_get_P(Load* load)
    REAL LOAD_get_Q(Load* load)
    Load* LOAD_get_next(Load* load)
    Load* LOAD_new()
    void LOAD_set_P(Load* load, REAL P)
    void LOAD_set_Q(Load* load, REAL Q)

