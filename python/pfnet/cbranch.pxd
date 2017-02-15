#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cvec

cdef extern from "pfnet/branch.h":

    ctypedef struct Branch
    ctypedef struct Bus
    ctypedef double REAL

    cdef char BRANCH_VAR_RATIO
    cdef char BRANCH_VAR_RATIO_DEV
    cdef char BRANCH_VAR_PHASE

    cdef double BRANCH_INF_RATIO
    cdef double BRANCH_INF_FLOW

    cdef char BRANCH_PROP_ANY
    cdef char BRANCH_PROP_TAP_CHANGER
    cdef char BRANCH_PROP_TAP_CHANGER_V
    cdef char BRANCH_PROP_TAP_CHANGER_Q
    cdef char BRANCH_PROP_PHASE_SHIFTER
    cdef char BRANCH_PROP_NOT_OUT

    REAL BRANCH_get_sens_P_u_bound(Branch* br, int t)
    REAL BRANCH_get_sens_P_l_bound(Branch* br, int t)
    char BRANCH_get_obj_type(void* br)
    int BRANCH_get_num_periods(Branch* br)
    int BRANCH_get_index(Branch* br)
    int BRANCH_get_index_ratio(Branch* br, int t)
    int BRANCH_get_index_ratio_y(Branch* br, int t)
    int BRANCH_get_index_ratio_z(Branch* br, int t)
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
    void BRANCH_set_ratio_max(Branch* br, REAL ratio)
    void BRANCH_set_ratio_min(Branch* br, REAL ratio)
    void BRANCH_set_ratingA(Branch* br, REAL r)
    void BRANCH_set_ratingB(Branch* br, REAL r)
    void BRANCH_set_ratingC(Branch* br, REAL r)
