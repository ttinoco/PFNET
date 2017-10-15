#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/matrix.h":

    ctypedef struct Mat:
        pass

    ctypedef double REAL
    
    int MAT_get_size1(Mat* m)
    int MAT_get_size2(Mat* m)
    int MAT_get_nnz(Mat* m)
    int* MAT_get_row_array(Mat* m)
    int* MAT_get_col_array(Mat* m)
    REAL* MAT_get_data_array(Mat* m)
    Mat* MAT_new_from_arrays(int size1, int size2, int nnz, int* row, int* col, REAL* data)

    
