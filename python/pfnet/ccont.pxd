#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/contingency.h":

    ctypedef struct Cont
    ctypedef struct Net

    void CONT_apply(Cont* cont, Net* net)
    void CONT_clear(Cont* cont, Net* net)
    void CONT_add_branch_outage(Cont* cont, int br_index)
    void CONT_add_gen_outage(Cont* cont, int gen_index)
    void CONT_del(Cont* cont)
    void CONT_init(Cont* cont)
    int CONT_get_num_gen_outages(Cont* cont)
    int CONT_get_num_branch_outages(Cont* cont)
    int* CONT_get_branch_outages(Cont* cont)
    int* CONT_get_gen_outages(Cont* cont)
    bint CONT_has_gen_outage(Cont* cont, int gen_index)
    bint CONT_has_branch_outage(Cont* cont, int br_index)
    Cont* CONT_new()
    void CONT_show(Cont* cont)
    char* CONT_get_show_str(Cont* cont)
    char* CONT_get_json_string(Cont* cont)
           


