#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/vargen.h":

    ctypedef struct Vargen
    ctypedef struct Bus
    ctypedef double REAL
    
    cdef char VARGEN_VAR_P
    cdef char VARGEN_VAR_Q

    cdef char VARGEN_PROP_ANY

    char* VARGEN_get_name(Vargen* gen)
    int VARGEN_get_index(Vargen* gen)
    int VARGEN_get_index_P(Vargen* gen)
    int VARGEN_get_index_Q(Vargen* gen)
    Bus* VARGEN_get_bus(Vargen* gen)
    REAL VARGEN_get_P(Vargen* gen)
    REAL VARGEN_get_P_max(Vargen* gen)
    REAL VARGEN_get_P_min(Vargen* gen)
    REAL VARGEN_get_P_std(Vargen* gen)
    REAL VARGEN_get_Q(Vargen* gen)
    REAL VARGEN_get_Q_max(Vargen* gen)
    REAL VARGEN_get_Q_min(Vargen* gen)
    Vargen* VARGEN_get_next(Vargen* gen)
    bint VARGEN_has_flags(Vargen* gen, char flag_type, char mask)
    Vargen* VARGEN_new()
    void VARGEN_set_name(Vargen* gen, char* name)
    void VARGEN_set_P(Vargen* gen, REAL P)
    void VARGEN_set_P_max(Vargen* gen, REAL P_max)
    void VARGEN_set_P_min(Vargen* gen, REAL P_min)
    void VARGEN_set_P_std(Vargen* gen, REAL P_std)
    void VARGEN_set_Q(Vargen* gen, REAL Q)
    void VARGEN_set_Q_max(Vargen* gen, REAL Q)
    void VARGEN_set_Q_min(Vargen* gen, REAL Q)
