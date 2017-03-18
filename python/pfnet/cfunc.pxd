#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/pfnet.h":

    ctypedef struct Func
    ctypedef struct Net
    ctypedef struct Vec
    ctypedef struct Mat
    ctypedef double REAL
    ctypedef struct Branch
        
    void FUNC_del(Func* f)
    void FUNC_del_matvec(Func* f)
    REAL FUNC_get_weight(Func* f)
    REAL FUNC_get_phi(Func* f)
    Vec* FUNC_get_gphi(Func* f)
    Mat* FUNC_get_Hphi(Func* f)
    int FUNC_get_Hphi_nnz(Func* f)
    Net* FUNC_get_network(Func* f)
    Func* FUNC_get_next(Func* f)
    Func* FUNC_new(REAL weight, Net* net)
    void FUNC_count(Func* f)
    void FUNC_allocate(Func* f)
    void FUNC_analyze(Func* f)
    void FUNC_eval(Func* f, Vec* var_values)
    bint FUNC_has_error(Func* f)
    void FUNC_clear_error(Func * f)
    char* FUNC_get_name(Func* f)
    char* FUNC_get_error_string(Func* f)
    void FUNC_update_network(Func* f)
    void FUNC_set_name(Func* f, char*)
    void FUNC_set_phi(Func* f, REAL phi)
    void FUNC_set_gphi(Func* f, Vec* gphi)
    void FUNC_set_Hphi(Func* f, Mat* Hphi)
    void FUNC_set_Hphi_nnz(Func* f, int nnz)

    void FUNC_set_func_init(Func* f, void (*func)(Func* f))
    void FUNC_set_func_count_step(Func* f, void (*func)(Func* f, Branch* br, int t))
    void FUNC_set_func_allocate(Func* f, void (*func)(Func* f))
    void FUNC_set_func_clear(Func* f, void (*func)(Func* f))
    void FUNC_set_func_analyze_step(Func* f, void (*func)(Func* f, Branch* br, int t))
    void FUNC_set_func_eval_step(Func* f, void (*func)(Func* f, Branch* br, int t, Vec* v))

    Func* FUNC_GEN_COST_new(REAL w, Net* net)
    Func* FUNC_LOAD_UTIL_new(REAL w, Net* net)
    Func* FUNC_NETCON_COST_new(REAL w, Net* net)
    Func* FUNC_REG_PHASE_new(REAL w, Net* net)
    Func* FUNC_REG_PQ_new(REAL w, Net* net)
    Func* FUNC_REG_RATIO_new(REAL w, Net* net)
    Func* FUNC_REG_SUSC_new(REAL w, Net* net)
    Func* FUNC_REG_VANG_new(REAL w, Net* net)
    Func* FUNC_REG_VMAG_new(REAL w, Net* net)
    Func* FUNC_SLIM_VMAG_new(REAL w, Net* net)
    Func* FUNC_SP_CONTROLS_new(REAL w, Net* net)

    void* FUNC_get_data(Func* f)
    void FUNC_set_data(Func* f, void* data)
    
    char* FUNC_get_bus_counted(Func* f)
    int FUNC_get_bus_counted_size(Func* f)
