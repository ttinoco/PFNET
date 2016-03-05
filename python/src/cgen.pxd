#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/gen.h":

    ctypedef struct Gen
    ctypedef struct Bus
    ctypedef double REAL
    
    cdef char GEN_VAR_P
    cdef char GEN_VAR_Q

    cdef double GEN_INF_P
    cdef double GEN_INF_Q

    cdef char GEN_PROP_ANY
    cdef char GEN_PROP_SLACK
    cdef char GEN_PROP_REG
    cdef char GEN_PROP_NOT_REG
    cdef char GEN_PROP_NOT_SLACK
    cdef char GEN_PROP_P_ADJUST

    REAL GEN_get_sens_P_u_bound(Gen* gen)
    REAL GEN_get_sens_P_l_bound(Gen* gen)
    REAL GEN_get_P_cost(Gen* gen)
    REAL GEN_get_cost_coeff_Q0(Gen* gen)
    REAL GEN_get_cost_coeff_Q1(Gen* gen)
    REAL GEN_get_cost_coeff_Q2(Gen* gen)
    char GEN_get_obj_type(void* gen)
    int GEN_get_index(Gen* gen)
    int GEN_get_index_P(Gen* gen)
    int GEN_get_index_Q(Gen* gen)
    Bus* GEN_get_bus(Gen* gen)
    Bus* GEN_get_reg_bus(Gen* gen)
    REAL GEN_get_P(Gen* gen)
    REAL GEN_get_P_max(Gen* gen)
    REAL GEN_get_P_min(Gen* gen)
    REAL GEN_get_Q(Gen* gen)
    REAL GEN_get_Q_max(Gen* gen)
    REAL GEN_get_Q_min(Gen* gen)
    Gen* GEN_get_next(Gen* gen)
    Gen* GEN_get_reg_next(Gen* gen)
    bint GEN_is_P_adjustable(Gen* gen)
    bint GEN_is_regulator(Gen* gen)
    bint GEN_is_slack(Gen* gen)
    bint GEN_has_flags(Gen* gen, char flag_type, char mask)
    Gen* GEN_new()
    void GEN_set_P_min(Gen* gen, REAL P_min)
    void GEN_set_P_max(Gen* gen, REAL P_max)
    void GEN_set_cost_coeff_Q0(Gen* gen, REAL c)
    void GEN_set_cost_coeff_Q1(Gen* gen, REAL c)
    void GEN_set_cost_coeff_Q2(Gen* gen, REAL c)
