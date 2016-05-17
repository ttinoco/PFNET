#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/bat.h":

    ctypedef struct Bat
    ctypedef struct Bus
    ctypedef double REAL

    cdef char BAT_VAR_P
    cdef char BAT_VAR_E

    cdef double BAT_INF_P
    cdef double BAT_INF_E
       
    cdef char BAT_PROP_ANY
    
    char BAT_get_obj_type(void* bat)
    int BAT_get_index(Bat* bat)
    int BAT_get_index_P(Bat* bat)
    int BAT_get_index_E(Bat* bat)
    Bus* BAT_get_bus(Bat* bat)
    REAL BAT_get_P(Bat* bat)
    REAL BAT_get_P_max(Bat* bat)
    REAL BAT_get_P_min(Bat* bat)
    REAL BAT_get_E(Bat* bat)
    REAL BAT_get_E_max(Bat* bat)
    REAL BAT_get_eta_c(Bat* bat)
    REAL BAT_get_eta_d(Bat* bat)
    Bat* BAT_get_next(Bat* bat)
    bint BAT_has_flags(Bat* bat, char flag_type, char mask)
    Bat* BAT_new()
    void BAT_set_P(Bat* bat, REAL P)
    void BAT_set_P_min(Bat* bat, REAL P_min)
    void BAT_set_P_max(Bat* bat, REAL P_max)
    void BAT_set_E(Bat* bat, REAL E)
    void BAT_set_E_max(Bat* bat, REAL E_max)
    void BAT_set_eta_c(Bat* bat, REAL eta_c)
    void BAT_set_eta_d(Bat* bat, REAL eta_d)
