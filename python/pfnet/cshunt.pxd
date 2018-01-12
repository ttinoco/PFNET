#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
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

    char SHUNT_get_flags_vars(Shunt* shunt)
    char SHUNT_get_flags_fixed(Shunt* shunt)
    char SHUNT_get_flags_bounded(Shunt* shunt)
    char SHUNT_get_flags_sparse(Shunt* shunt)

    REAL* SHUNT_get_sens_b_u_bound_array(Shunt* shunt)
    REAL* SHUNT_get_sens_b_l_bound_array(Shunt* shunt)
  
    char* SHUNT_get_name(Shunt* shunt)    
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
    char* SHUNT_get_json_string(Shunt* shunt, char* output)
    char* SHUNT_get_var_info_string(Shunt* shunt, int index)
    bint SHUNT_is_equal(Shunt* load, Shunt* other)
    bint SHUNT_is_fixed(Shunt* shunt)
    bint SHUNT_is_switched_v(Shunt* shunt)
    bint SHUNT_has_flags(Shunt* shunt, char flag_type, char mask)
    Shunt* SHUNT_new(int num_periods)
    Shunt* SHUNT_array_new(int size, int num_periods)
    void SHUNT_array_del(Shunt* shunt_array, int size)
    void SHUNT_set_name(Shunt* shunt, char* name)
    void SHUNT_set_bus(Shunt* shunt, Bus* bus)
    void SHUNT_set_reg_bus(Shunt* shunt, Bus* reg_bus)
    void SHUNT_set_g(Shunt* shunt, REAL g)
    void SHUNT_set_b(Shunt* shunt, REAL b, int t)
    void SHUNT_set_b_max(Shunt* shunt, REAL b_max)
    void SHUNT_set_b_min(Shunt* shunt, REAL b_min)
    void SHUNT_set_b_values(Shunt* shunt, REAL* values, int num)
