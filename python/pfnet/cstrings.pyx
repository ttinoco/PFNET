#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cobjs
cimport cflags

# Objects
str2obj = {'all' : cobjs.OBJ_ALL,
           'bus' : cobjs.OBJ_BUS,
           'generator' : cobjs.OBJ_GEN,
           'branch' : cobjs.OBJ_BRANCH,
           'shunt' : cobjs.OBJ_SHUNT,
           'load' : cobjs.OBJ_LOAD,
           'variable generator' : cobjs.OBJ_VARGEN,
           'battery' : cobjs.OBJ_BAT,
           'unknown' : cobjs.OBJ_UNKNOWN}

obj2str = dict([(v,k) for k,v in str2obj.items()])

# Flags
str2flag = {'none' : cflags.FLAG_NONE,
            'variable' : cflags.FLAG_VARS,
            'fixed' : cflags.FLAG_FIXED,
            'bounded' : cflags.FLAG_BOUNDED,
            'sparse' : cflags.FLAG_SPARSE}

flag2str = dict([(v,k) for k,v in str2flag.items()])

# Variables
ALL_VARS = cflags.ALL_VARS
