#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/gen.h":

    ctypedef struct Gen
    ctypedef struct Bus
    ctypedef double REAL
    ctypedef char BOOL
    
    cdef char GEN_VAR_P
    cdef char GEN_VAR_Q

    cdef double GEN_INF_P
    cdef double GEN_INF_Q

    cdef char GEN_PROP_ANY
    cdef char GEN_PROP_SLACK
    cdef char GEN_PROP_REG
    cdef char GEN_PROP_NOT_REG
    cdef char GEN_PROP_NOT_SLACK
    cdef char GEN_PROP_NOT_OUT
    cdef char GEN_PROP_P_ADJUST

    REAL GEN_get_sens_P_u_bound(Gen* gen, int t)
    REAL GEN_get_sens_P_l_bound(Gen* gen, int t)
    REAL GEN_get_P_cost(Gen* gen, int t)
    REAL GEN_get_cost_coeff_Q0(Gen* gen)
    REAL GEN_get_cost_coeff_Q1(Gen* gen)
    REAL GEN_get_cost_coeff_Q2(Gen* gen)
    int GEN_get_num_periods(Gen* gen)
    char GEN_get_obj_type(void* gen)
    int GEN_get_index(Gen* gen)
    int GEN_get_index_P(Gen* gen, int t)
    int GEN_get_index_Q(Gen* gen, int t)
    Bus* GEN_get_bus(Gen* gen)
    Bus* GEN_get_reg_bus(Gen* gen)
    REAL GEN_get_P(Gen* gen, int t) 
    REAL GEN_get_dP_max(Gen* gen)
    REAL GEN_get_P_max(Gen* gen)
    REAL GEN_get_P_min(Gen* gen)
    REAL GEN_get_P_prev(Gen* gen)
    REAL GEN_get_Q(Gen* gen, int t)
    REAL GEN_get_Q_max(Gen* gen)
    REAL GEN_get_Q_min(Gen* gen)
    Gen* GEN_get_next(Gen* gen)
    Gen* GEN_get_reg_next(Gen* gen)
    char* GEN_get_json_string(Gen* gen, char* output)
    char* GEN_get_var_info_string(Gen* gen, int index)
    bint GEN_is_equal(Gen* gen, Gen* other)
    bint GEN_is_on_outage(Gen* gen)
    bint GEN_is_P_adjustable(Gen* gen)
    bint GEN_is_regulator(Gen* gen)
    bint GEN_is_slack(Gen* gen)
    bint GEN_has_flags(Gen* gen, char flag_type, char mask)
    Gen* GEN_new(int num_periods)
    Gen* GEN_array_new(int size, int num_periods)
    void GEN_set_bus(Gen* gen, Bus* bus)
    void GEN_set_reg_bus(Gen* gen, Bus* reg_bus)
    void GEN_set_P(Gen* gen, REAL P, int t)
    void GEN_set_P_min(Gen* gen, REAL P)
    void GEN_set_P_prev(Gen* gen, REAL P)
    void GEN_set_P_max(Gen* gen, REAL P)
    void GEN_set_dP_max(Gen* gen, REAL P)
    void GEN_set_Q(Gen* gen, REAL Q, int t)
    void GEN_set_Q_max(Gen* gen, REAL Q)
    void GEN_set_Q_min(Gen* gen, REAL Q)
    void GEN_set_cost_coeff_Q0(Gen* gen, REAL c)
    void GEN_set_cost_coeff_Q1(Gen* gen, REAL c)
    void GEN_set_cost_coeff_Q2(Gen* gen, REAL c)
