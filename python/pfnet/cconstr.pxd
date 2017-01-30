#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/constr.h":

    ctypedef struct Constr
    ctypedef struct Net
    ctypedef struct Vec
    ctypedef struct Mat

    ctypedef double REAL
    
    cdef char CONSTR_TYPE_PF
    cdef char CONSTR_TYPE_DCPF
    cdef char CONSTR_TYPE_LINPF
    cdef char CONSTR_TYPE_FIX
    cdef char CONSTR_TYPE_BOUND
    cdef char CONSTR_TYPE_PAR_GEN_P
    cdef char CONSTR_TYPE_PAR_GEN_Q
    cdef char CONSTR_TYPE_REG_GEN
    cdef char CONSTR_TYPE_REG_TRAN
    cdef char CONSTR_TYPE_REG_SHUNT
    cdef char CONSTR_TYPE_DC_FLOW_LIM
    cdef char CONSTR_TYPE_LBOUND
    cdef char CONSTR_TYPE_GEN_RAMP
    
    void CONSTR_combine_H(Constr* c, Vec* coeff, bint ensure_psd)
    void CONSTR_del(Constr* c)
    void CONSTR_del_matvec(Constr* constr)
    Constr* CONSTR_new(int type, Net* net)
    int CONSTR_get_Acounter(Constr* c)
    int CONSTR_get_Gcounter(Constr* c)
    int CONSTR_get_Jcounter(Constr* c)
    int CONSTR_get_Aconstr_index(Constr* c)
    int CONSTR_get_Gconstr_index(Constr* c)
    int CONSTR_get_Jconstr_index(Constr* c)
    Vec* CONSTR_get_f(Constr* c)
    Mat* CONSTR_get_J(Constr* c)
    Vec* CONSTR_get_b(Constr* c)
    Mat* CONSTR_get_A(Constr* c)
    Vec* CONSTR_get_l(Constr* c)
    Vec* CONSTR_get_u(Constr* c)
    Mat* CONSTR_get_G(Constr* c)
    Mat* CONSTR_get_H_single(Constr* c, int i)
    Mat* CONSTR_get_H_combined(Constr* c)
    int CONSTR_get_type(Constr* c)
    Constr* CONSTR_get_next(Constr* c)
    void CONSTR_count(Constr* c)
    void CONSTR_allocate(Constr* c)
    void CONSTR_analyze(Constr* c)
    void CONSTR_eval(Constr* c, Vec* values)
    void CONSTR_store_sens(Constr* c, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl)
    bint CONSTR_has_error(Constr* c)
    void CONSTR_clear_error(Constr * c)
    char* CONSTR_get_error_string(Constr* c)
    void CONSTR_update_network(Constr* c)
    
