#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/func.h":

    ctypedef struct Func
    ctypedef struct Net
    ctypedef struct Vec
    ctypedef struct Mat
    ctypedef double REAL
    
    cdef char FUNC_TYPE_REG_VMAG
    cdef char FUNC_TYPE_REG_VANG
    cdef char FUNC_TYPE_REG_PQ
    cdef char FUNC_TYPE_REG_RATIO
    cdef char FUNC_TYPE_REG_PHASE
    cdef char FUNC_TYPE_REG_SUSC
    cdef char FUNC_TYPE_GEN_COST
    cdef char FUNC_TYPE_SP_CONTROLS
    cdef char FUNC_TYPE_SLIM_VMAG
    cdef char FUNC_TYPE_LOAD_UTIL
    
    void FUNC_del(Func* f)
    void FUNC_del_matvec(Func* f)
    int FUNC_get_type(Func* f)
    REAL FUNC_get_weight(Func* f)
    REAL FUNC_get_phi(Func* f)
    Vec* FUNC_get_gphi(Func* f)
    Mat* FUNC_get_Hphi(Func* f)
    int FUNC_get_Hcounter(Func* f)
    Func* FUNC_get_next(Func* f)
    Func* FUNC_new(int type, REAL weight, Net* net)
    void FUNC_count(Func* f)
    void FUNC_allocate(Func* f)
    void FUNC_analyze(Func* f)
    void FUNC_eval(Func* f, Vec* var_values)
    bint FUNC_has_error(Func* f)
    void FUNC_clear_error(Func * f)
    char* FUNC_get_error_string(Func* f)
    void FUNC_update_network(Func* f)
