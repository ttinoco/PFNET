#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/contingency.h":

    ctypedef struct Cont
    ctypedef struct Gen
    ctypedef struct Branch

    void CONT_apply(Cont* cont)
    void CONT_clear(Cont* cont)
    void CONT_add_branch_outage(Cont* cont, Branch* br)
    void CONT_add_gen_outage(Cont* cont, Gen* gen)
    void CONT_del(Cont* cont)
    void CONT_init(Cont* cont)
    int CONT_get_num_gen_outages(Cont* cont)
    int CONT_get_num_branch_outages(Cont* cont)
    bint CONT_has_gen_outage(Cont* cont, Gen* gen)
    bint CONT_has_branch_outage(Cont* cont, Branch* br)
    Cont* CONT_new()
    void CONT_show(Cont* cont)


        


