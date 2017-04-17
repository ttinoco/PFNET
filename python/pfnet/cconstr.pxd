#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/constr.h":

    ctypedef struct Constr
    ctypedef struct Net
    ctypedef struct Vec
    ctypedef struct Mat
    ctypedef struct Branch
    ctypedef double REAL
        
    void CONSTR_combine_H(Constr* c, Vec* coeff, bint ensure_psd)
    void CONSTR_del(Constr* c)
    void CONSTR_del_matvec(Constr* constr)
    Constr* CONSTR_new(Net* net)
    int CONSTR_get_A_nnz(Constr* c)
    int CONSTR_get_G_nnz(Constr* c)
    int CONSTR_get_J_nnz(Constr* c)
    int CONSTR_get_A_row(Constr* c)
    int CONSTR_get_G_row(Constr* c)
    int CONSTR_get_J_row(Constr* c)
    char* CONSTR_get_name(Constr* c)
    Vec* CONSTR_get_f(Constr* c)
    Mat* CONSTR_get_J(Constr* c)
    Mat* CONSTR_get_Jbar(Constr* c)
    Vec* CONSTR_get_b(Constr* c)
    Mat* CONSTR_get_A(Constr* c)
    Vec* CONSTR_get_l(Constr* c)
    Vec* CONSTR_get_u(Constr* c)
    Mat* CONSTR_get_G(Constr* c)
    Mat* CONSTR_get_Gbar(Constr* c)
    Mat* CONSTR_get_H_single(Constr* c, int i)
    Mat* CONSTR_get_H_combined(Constr* c)
    int CONSTR_get_type(Constr* c)
    Constr* CONSTR_get_next(Constr* c)
    void CONSTR_init(Constr* c)
    void CONSTR_count(Constr* c)
    void CONSTR_allocate(Constr* c)
    void CONSTR_analyze(Constr* c)
    void CONSTR_eval(Constr* c, Vec* values)
    void CONSTR_store_sens(Constr* c, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl)
    bint CONSTR_has_error(Constr* c)
    void CONSTR_clear_error(Constr * c)
    char* CONSTR_get_error_string(Constr* c)
    void CONSTR_update_network(Constr* c)
    int CONSTR_get_num_extra_vars(Constr* c)
    Net* CONSTR_get_network(Constr* c)

    void CONSTR_set_name(Constr* f, char*)
    void CONSTR_set_A_nnz(Constr* c, int nnz)
    void CONSTR_set_G_nnz(Constr* c, int nnz)
    void CONSTR_set_Gbar_nnz(Constr* c, int nnz)
    void CONSTR_set_J_nnz(Constr* c, int nnz)
    void CONSTR_set_Jbar_nnz(Constr* c, int nnz)
    void CONSTR_set_H_nnz(Constr* c, int* nnz, int size)
    void CONSTR_set_A_row(Constr* c, int index)
    void CONSTR_set_G_row(Constr* c, int index)
    void CONSTR_set_J_row(Constr* c, int index)
    void CONSTR_set_b(Constr* c, Vec* b)
    void CONSTR_set_A(Constr* c, Mat* A)
    void CONSTR_set_l(Constr* c, Vec* l)
    void CONSTR_set_u(Constr* c, Vec* u)
    void CONSTR_set_G(Constr* c, Mat* G)
    void CONSTR_set_f(Constr* c, Vec* f)
    void CONSTR_set_J(Constr* c, Mat* J)

    void CONSTR_set_func_init(Constr* c, void (*func)(Constr* c))
    void CONSTR_set_func_count_step(Constr* c, void (*func)(Constr* c, Branch* br, int t))
    void CONSTR_set_func_allocate(Constr* c, void (*func)(Constr* c))
    void CONSTR_set_func_clear(Constr* c, void (*func)(Constr* c))
    void CONSTR_set_func_analyze_step(Constr* c, void (*func)(Constr* c, Branch* br, int t))
    void CONSTR_set_func_eval_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* v))
    void CONSTR_set_func_store_sens_step(Constr* c, void (*func)(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl))

    Constr* CONSTR_ACPF_new(Net* net)
    Constr* CONSTR_DCPF_new(Net* net)
    Constr* CONSTR_LINPF_new(Net* net)
    Constr* CONSTR_FIX_new(Net* net)
    Constr* CONSTR_LBOUND_new(Net* net)
    Constr* CONSTR_NBOUND_new(Net* net)
    Constr* CONSTR_PAR_GEN_P_new(Net* net)
    Constr* CONSTR_PAR_GEN_Q_new(Net* net)
    Constr* CONSTR_GEN_RAMP_new(Net* net)
    Constr* CONSTR_REG_GEN_new(Net* net)
    Constr* CONSTR_REG_TRAN_new(Net* net)
    Constr* CONSTR_REG_SHUNT_new(Net* net)
    Constr* CONSTR_DC_FLOW_LIM_new(Net* net)
    Constr* CONSTR_AC_FLOW_LIM_new(Net* net)
    Constr* CONSTR_BAT_DYN_new(Net* net)

    void* CONSTR_get_data(Constr* c)
    void CONSTR_set_data(Constr* c, void* data)
    
    char* CONSTR_get_bus_counted(Constr* c)
    int CONSTR_get_bus_counted_size(Constr* c)
