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

    cdef char LOAD_VAR_P

    cdef double LOAD_INF_P
       
    cdef char LOAD_PROP_ANY
    cdef char LOAD_PROP_P_ADJUST
 
    REAL LOAD_get_sens_P_u_bound(Load* load)
    REAL LOAD_get_sens_P_l_bound(Load* load)
    REAL LOAD_get_P_util(Load* load)
    REAL LOAD_get_util_coeff_Q0(Load* load)
    REAL LOAD_get_util_coeff_Q1(Load* load)
    REAL LOAD_get_util_coeff_Q2(Load* load)
    char LOAD_get_obj_type(void* load)
    int LOAD_get_index(Load* load)
    int LOAD_get_index_P(Load* load)
    Bus* LOAD_get_bus(Load* load)
    REAL LOAD_get_P(Load* load)
    REAL LOAD_get_P_max(Load* load)
    REAL LOAD_get_P_min(Load* load)
    REAL LOAD_get_Q(Load* load)
    Load* LOAD_get_next(Load* load)
    bint LOAD_is_P_adjustable(Load* load)
    bint LOAD_has_flags(Load* load, char flag_type, char mask)
    Load* LOAD_new()
    void LOAD_set_P(Load* load, REAL P)
    void LOAD_set_P_min(Load* load, REAL P_min)
    void LOAD_set_P_max(Load* load, REAL P_max)
    void LOAD_set_Q(Load* load, REAL Q)
    void LOAD_set_util_coeff_Q0(Load* load, REAL c)
    void LOAD_set_util_coeff_Q1(Load* load, REAL c)
    void LOAD_set_util_coeff_Q2(Load* load, REAL c)

