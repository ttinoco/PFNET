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
cimport cconstants

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
str2flag = {'variable' : cflags.FLAG_VARS,
            'fixed' : cflags.FLAG_FIXED,
            'bounded' : cflags.FLAG_BOUNDED,
            'sparse' : cflags.FLAG_SPARSE}

# Variable values
str2const = {'current' : cconstants.CURRENT,
             'upper limits' : cconstants.UPPER_LIMITS,
             'lower limits' : cconstants.LOWER_LIMITS}

# Variables
str2q_bus = {'all' : cflags.ALL_VARS, 
             'voltage magnitude' : cbus.BUS_VAR_VMAG, 
             'voltage angle' : cbus.BUS_VAR_VANG,
             'voltage magnitude deviation' : cbus.BUS_VAR_VDEV,
             'voltage magnitude violation' : cbus.BUS_VAR_VVIO}

str2q_branch = {'all' : cflags.ALL_VARS,
                'tap ratio' : cbranch.BRANCH_VAR_RATIO,
                'tap ratio deviation' : cbranch.BRANCH_VAR_RATIO_DEV,
                'phase shift' : cbranch.BRANCH_VAR_PHASE}

str2q_gen = {'all' : cflags.ALL_VARS,
             'active power' : cgen.GEN_VAR_P,
             'reactive power' : cgen.GEN_VAR_Q}

str2q_shunt = {'all' : cflags.ALL_VARS,
               'susceptance' : cshunt.SHUNT_VAR_SUSC,
               'susceptance deviation' : cshunt.SHUNT_VAR_SUSC_DEV}

str2q_load = {'all' : cflags.ALL_VARS,
              'active power' : cload.LOAD_VAR_P}

str2q_vargen = {'all' : cflags.ALL_VARS,
                'active power' : cvargen.VARGEN_VAR_P,
                'reactive power' : cvargen.VARGEN_VAR_Q}

str2q_bat = {'all' : cflags.ALL_VARS,
             'charging power' : cbat.BAT_VAR_P,
             'energy level' : cbat.BAT_VAR_E}

str2q = {'all' : {'all' : cflags.ALL_VARS},
         'bus' : str2q_bus,
         'branch' : str2q_branch,
         'generator' : str2q_gen,
         'shunt' : str2q_shunt,
         'load' : str2q_load,
         'variable generator' : str2q_vargen,
         'battery' : str2q_bat}



