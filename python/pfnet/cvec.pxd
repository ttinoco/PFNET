#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/vector.h":

    ctypedef struct Vec:
        pass

    ctypedef double REAL
    
    Vec* VEC_new_from_array(REAL* data, int size)
    int VEC_get_size(Vec* v)
    REAL* VEC_get_data(Vec* v)


    
