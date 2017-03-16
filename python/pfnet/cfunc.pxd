#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/func.h":

    ctypedef struct Func
    ctypedef struct Net
    ctypedef struct Vec
    ctypedef struct Mat
    ctypedef double REAL
        
    void FUNC_del(Func* f)
    void FUNC_del_matvec(Func* f)
    REAL FUNC_get_weight(Func* f)
    REAL FUNC_get_phi(Func* f)
    Vec* FUNC_get_gphi(Func* f)
    Mat* FUNC_get_Hphi(Func* f)
    int FUNC_get_Hcounter(Func* f)
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
