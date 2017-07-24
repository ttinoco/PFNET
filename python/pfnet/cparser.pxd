#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cdef extern from "pfnet/pfnet.h":

    ctypedef struct Net
    ctypedef struct Parser
    ctypedef double REAL

    Parser* PARSER_new()
    Parser* PARSER_new_for_file(char* f)
    void PARSER_init(Parser* p)
    Net* PARSER_parse(Parser* p, char* f, int num_periods)
    void PARSER_set(Parser* p, char* key, REAL value)
    void PARSER_show(Parser* p)
    void PARSER_write(Parser* p, Net* net, char* f)
    void PARSER_del(Parser* p)
    
    bint PARSER_has_error(Parser* p)
    void PARSER_clear_error(Parser* p)
    char* PARSER_get_error_string(Parser* p)
    void* PARSER_get_data(Parser* p)
    void PARSER_set_error(Parser* p, char* string)
    void PARSER_set_data(Parser* p, void* data)

    void PARSER_set_func_init(Parser* p, void (*func)(Parser* p))
    void PARSER_set_func_parse(Parser* p, Net* (*func)(Parser* p, char* f, int n))
    void PARSER_set_func_set(Parser* p, void (*func)(Parser* p, char* key, REAL v))
    void PARSER_set_func_show(Parser* p, void (*func)(Parser* p))
    void PARSER_set_func_write(Parser* p, void (*func)(Parser* p, Net* net, char* f))
    void PARSER_set_func_free(Parser* p, void (*func)(Parser* p))

    Parser* MAT_PARSER_new()
    Parser* ART_PARSER_new()
    Parser* RAW_PARSER_new()
    Parser* JSON_PARSER_new()
    
    

