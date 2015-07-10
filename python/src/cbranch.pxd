#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/branch.h":

    ctypedef struct Branch 
    ctypedef struct Bus
    ctypedef double REAL

    cdef char BRANCH_VAR_RATIO
    cdef char BRANCH_VAR_RATIO_DEV
    cdef char BRANCH_VAR_PHASE

    cdef char BRANCH_PROP_ANY
    cdef char BRANCH_PROP_TAP_CHANGER
    cdef char BRANCH_PROP_TAP_CHANGER_V
    cdef char BRANCH_PROP_TAP_CHANGER_Q
    cdef char BRANCH_PROP_PHASE_SHIFTER
    
    int BRANCH_get_index(Branch* br)
    int BRANCH_get_index_ratio(Branch* br)
    int BRANCH_get_index_ratio_y(Branch* br)
    int BRANCH_get_index_ratio_z(Branch* br)
    int BRANCH_get_index_phase(Branch* br)
    Bus* BRANCH_get_bus_from(Branch* br)
    Bus* BRANCH_get_bus_to(Branch* br)
    Bus* BRANCH_get_reg_bus(Branch* br)
    REAL BRANCH_get_ratio(Branch* br)
    REAL BRANCH_get_ratio_max(Branch* br)
    REAL BRANCH_get_ratio_min(Branch* br)
    REAL BRANCH_get_b(Branch* br)
    REAL BRANCH_get_b_from(Branch* br)
    REAL BRANCH_get_b_to(Branch* br)
    REAL BRANCH_get_g(Branch* br)
    REAL BRANCH_get_g_from(Branch* br)
    REAL BRANCH_get_g_to(Branch* br)
    REAL BRANCH_get_phase(Branch* br)
    REAL BRANCH_get_phase_max(Branch* br)
    REAL BRANCH_get_phase_min(Branch* br)
    Branch* BRANCH_get_reg_next(Branch* br)
    Branch* BRANCH_get_from_next(Branch* br)
    Branch* BRANCH_get_to_next(Branch* br)
    bint BRANCH_has_pos_ratio_v_sens(Branch* br)
    bint BRANCH_is_fixed_tran(Branch* br)
    bint BRANCH_is_line(Branch* br)
    bint BRANCH_is_phase_shifter(Branch* br)
    bint BRANCH_is_tap_changer(Branch* br)
    bint BRANCH_is_tap_changer_v(Branch* br)
    bint BRANCH_is_tap_changer_Q(Branch* br)
    bint BRANCH_has_flags(Branch* br, char flag_type, char mask)
    Branch* BRANCH_new()
    void BRANCH_set_ratio_max(Branch* br, REAL ratio)
    void BRANCH_set_ratio_min(Branch* br, REAL ratio)
    
