#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cvec

cdef extern from "pfnet/shunt.h":

    ctypedef struct Shunt
    ctypedef struct Bus
    ctypedef double REAL
    
    cdef char SHUNT_VAR_SUSC

    cdef double SHUNT_INF_SUSC

    cdef char SHUNT_PROP_ANY
    cdef char SHUNT_PROP_SWITCHED_V
    
    int SHUNT_get_num_periods(Shunt* shunt)
    char SHUNT_get_obj_type(void* shunt)
    int SHUNT_get_index(Shunt* shunt)
    int SHUNT_get_index_b(Shunt* shunt, int t)
    Bus* SHUNT_get_bus(Shunt* shunt)
    Bus* SHUNT_get_reg_bus(Shunt* shunt)
    REAL SHUNT_get_g(Shunt* shunt)
    REAL SHUNT_get_b(Shunt* shunt, int t)
    REAL SHUNT_get_b_max(Shunt* shunt)
    REAL SHUNT_get_b_min(Shunt* shunt)
    REAL* SHUNT_get_b_values(Shunt* shunt)
    int SHUNT_get_num_b_values(Shunt* shunt)
    Shunt* SHUNT_get_next(Shunt* shunt)
    Shunt* SHUNT_get_reg_next(Shunt* shunt)
    bint SHUNT_is_fixed(Shunt* shunt)
    bint SHUNT_is_switched_v(Shunt* shunt)
    bint SHUNT_has_flags(Shunt* shunt, char flag_type, char mask)
    Shunt* SHUNT_new(int num_periods)
    Shunt* SHUNT_array_new(int size, int num_periods)
    void SHUNT_set_bus(Shunt* shunt, Bus* bus)
    void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus)
    void SHUNT_set_g(Shunt* shunt, REAL g)
    void SHUNT_set_b(Shunt* shunt, REAL b, int t)
    void SHUNT_set_b_max(Shunt* shunt, REAL b_max)
    void SHUNT_set_b_min(Shunt* shunt, REAL b_min)
    void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num, REAL norm)
