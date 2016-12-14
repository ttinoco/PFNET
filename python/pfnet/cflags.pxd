#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/flag_types.h":

    cdef char ALL_VARS
    
    cdef char FLAG_NONE
    cdef char FLAG_VARS
    cdef char FLAG_FIXED
    cdef char FLAG_BOUNDED 
    cdef char FLAG_SPARSE



    
