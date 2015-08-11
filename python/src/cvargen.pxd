#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/vargen.h":

    ctypedef struct Vargen
    ctypedef struct Bus
    ctypedef double REAL
    
    cdef char VARGEN_VAR_P

    cdef char VARGEN_PROP_ANY

    Vargen* VARGEN_array_new(int num)
    int VARGEN_get_index(Vargen* gen)
    int VARGEN_get_index_P(Vargen* gen)
    Bus* VARGEN_get_bus(Vargen* gen)
    REAL VARGEN_get_P(Vargen* gen)
    REAL VARGEN_get_P_max(Vargen* gen)
    Vargen* VARGEN_get_next(Vargen* gen)
    bint VARGEN_has_flags(Vargen* gen, char flag_type, char mask)
    Vargen* VARGEN_new()
