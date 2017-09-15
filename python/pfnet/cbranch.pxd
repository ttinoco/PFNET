#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cvec

cdef extern from "pfnet/branch.h":

    ctypedef struct Branch
    ctypedef struct Bus
    ctypedef double REAL
    ctypedef char BOOL

    cdef char BRANCH_VAR_RATIO
    cdef char BRANCH_VAR_PHASE

    cdef double BRANCH_INF_RATIO
    cdef double BRANCH_INF_PHASE
    cdef double BRANCH_INF_FLOW

    cdef char BRANCH_PROP_ANY
    cdef char BRANCH_PROP_TAP_CHANGER
    cdef char BRANCH_PROP_TAP_CHANGER_V
    cdef char BRANCH_PROP_TAP_CHANGER_Q
    cdef char BRANCH_PROP_PHASE_SHIFTER
    cdef char BRANCH_PROP_NOT_OUT
    
    cdef int BRANCH_TYPE_LINE
    cdef int BRANCH_TYPE_TRAN_FIXED
    cdef int BRANCH_TYPE_TRAN_TAP_V
    cdef int BRANCH_TYPE_TRAN_TAP_Q
    cdef int BRANCH_TYPE_TRAN_PHASE

    char BRANCH_get_flags_vars(Branch* br)
    char BRANCH_get_flags_fixed(Branch* br)
    char BRANCH_get_flags_bounded(Branch* br)
    char BRANCH_get_flags_sparse(Branch* br)

    REAL BRANCH_get_sens_P_u_bound(Branch* br, int t)
    REAL BRANCH_get_sens_P_l_bound(Branch* br, int t)
    char BRANCH_get_obj_type(void* br)
    int BRANCH_get_num_periods(Branch* br)
    int BRANCH_get_index(Branch* br)
    int BRANCH_get_index_ratio(Branch* br, int t)
    int BRANCH_get_index_phase(Branch* br, int t)
    Bus* BRANCH_get_bus_k(Branch* br)
    Bus* BRANCH_get_bus_m(Branch* br)
    Bus* BRANCH_get_reg_bus(Branch* br)
    REAL BRANCH_get_ratio(Branch* br, int t)
    REAL BRANCH_get_ratio_max(Branch* br)
    REAL BRANCH_get_ratio_min(Branch* br)
    REAL BRANCH_get_b(Branch* br)
    REAL BRANCH_get_b_k(Branch* br)
    REAL BRANCH_get_b_m(Branch* br)
    REAL BRANCH_get_g(Branch* br)
    REAL BRANCH_get_g_k(Branch* br)
    REAL BRANCH_get_g_m(Branch* br)
    REAL BRANCH_get_phase(Branch* br, int t)
    REAL BRANCH_get_phase_max(Branch* br)
    REAL BRANCH_get_phase_min(Branch* br)
    REAL BRANCH_get_P_max(Branch* br)
    REAL BRANCH_get_P_min(Branch* br)
    REAL BRANCH_get_Q_max(Branch* br)
    REAL BRANCH_get_Q_min(Branch* br)
    REAL BRANCH_get_P_km(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_km(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_P_mk(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_mk(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_P_km_series(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_km_series(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_P_mk_series(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_mk_series(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_P_k_shunt(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_k_shunt(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_P_m_shunt(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_Q_m_shunt(Branch* br, cvec.Vec* values, int t)
    REAL BRANCH_get_ratingA(Branch* br)
    REAL BRANCH_get_ratingB(Branch* br)
    REAL BRANCH_get_ratingC(Branch* br)
    REAL BRANCH_get_P_km_DC(Branch* br, int t)
    REAL BRANCH_get_P_mk_DC(Branch* br, int t)
    Branch* BRANCH_get_reg_next(Branch* br)
    Branch* BRANCH_get_next_k(Branch* br)
    Branch* BRANCH_get_next_m(Branch* br)
    char* BRANCH_get_json_string(Branch* br, char* output)
    char* BRANCH_get_var_info_string(Branch* br, int index)  
    bint BRANCH_has_pos_ratio_v_sens(Branch* br)
    bint BRANCH_is_equal(Branch* br, Branch* other)
    bint BRANCH_is_on_outage(Branch* br)
    bint BRANCH_is_fixed_tran(Branch* br)
    bint BRANCH_is_line(Branch* br)
    bint BRANCH_is_phase_shifter(Branch* br)
    bint BRANCH_is_tap_changer(Branch* br)
    bint BRANCH_is_tap_changer_v(Branch* br)
    bint BRANCH_is_tap_changer_Q(Branch* br)
    bint BRANCH_has_flags(Branch* br, char flag_type, char mask)
    Branch* BRANCH_new(int num_periods)
    Branch* BRANCH_array_new(int size, int num_periods)
    void BRANCH_array_show(Branch* br_array, int size, int t)
    void BRANCH_set_type(Branch* br, int type)
    void BRANCH_set_bus_k(Branch* br, Bus* bus_k)
    void BRANCH_set_bus_m(Branch* br, Bus* bus_m)
    void BRANCH_set_reg_bus(Branch* br, Bus* reg_bus)
    void BRANCH_set_g(Branch* br, REAL g)
    void BRANCH_set_g_k(Branch* br, REAL g_k)
    void BRANCH_set_g_m(Branch* br, REAL g_m)
    void BRANCH_set_b(Branch* br, REAL b)
    void BRANCH_set_b_k(Branch* br, REAL b_k)
    void BRANCH_set_b_m(Branch* br, REAL b_m)
    void BRANCH_set_ratio(Branch* br, REAL ratio, int t)
    void BRANCH_set_ratio_max(Branch* br, REAL ratio)
    void BRANCH_set_ratio_min(Branch* br, REAL ratio)
    void BRANCH_set_pos_ratio_v_sens(Branch* br, BOOL flag)
    void BRANCH_set_phase(Branch* br, REAL phase, int t)
    void BRANCH_set_phase_max(Branch* br, REAL phase) 
    void BRANCH_set_phase_min(Branch* br, REAL phase) 
    void BRANCH_set_P_max(Branch* br, REAL P_max)
    void BRANCH_set_P_min(Branch* br, REAL P_min)
    void BRANCH_set_Q_max(Branch* br, REAL Q_max)
    void BRANCH_set_Q_min(Branch* br, REAL Q_min)
    void BRANCH_set_ratingA(Branch* br, REAL r)
    void BRANCH_set_ratingB(Branch* br, REAL r)
    void BRANCH_set_ratingC(Branch* br, REAL r)
    void BRANCH_show(Branch* br, int t)
