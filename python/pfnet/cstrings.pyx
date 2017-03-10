#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cobjs
cimport cflags
cimport cconstants
cimport cbus
cimport cbranch
cimport cgen
cimport cshunt
cimport cload
cimport cvargen
cimport cbat

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

# Quantities
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

# Properties
str2prop_bus = {'any' :  cbus.BUS_PROP_ANY,
                'slack' : cbus.BUS_PROP_SLACK,
                'regulated by generator' : cbus.BUS_PROP_REG_BY_GEN,
                'regulated by transformer' : cbus.BUS_PROP_REG_BY_TRAN,
                'regulated by shunt' : cbus.BUS_PROP_REG_BY_SHUNT,
                'not slack' : cbus.BUS_PROP_NOT_SLACK,
                'not regulated by generator' : cbus.BUS_PROP_NOT_REG_BY_GEN}

str2prop_branch = {'any' : cbranch.BRANCH_PROP_ANY,
                   'tap changer' : cbranch.BRANCH_PROP_TAP_CHANGER,
                   'tap changer - v' : cbranch.BRANCH_PROP_TAP_CHANGER_V,
                   'tap changer - Q' : cbranch.BRANCH_PROP_TAP_CHANGER_Q,
                   'phase shifter' : cbranch.BRANCH_PROP_PHASE_SHIFTER,
                   'not on outage' : cbranch.BRANCH_PROP_NOT_OUT}

str2prop_gen = {'any' : cgen.GEN_PROP_ANY,
                'slack' : cgen.GEN_PROP_SLACK,
                'regulator' : cgen.GEN_PROP_REG,
                'not slack' : cgen.GEN_PROP_NOT_SLACK,
                'not regulator' : cgen.GEN_PROP_NOT_REG,
                'not on outage' : cgen.GEN_PROP_NOT_OUT,
                'adjustable active power' : cgen.GEN_PROP_P_ADJUST}

str2prop_shunt =  {'any' : cshunt.SHUNT_PROP_ANY,
                   'switching - v' : cshunt.SHUNT_PROP_SWITCHED_V}

str2prop_load = {'any' : cload.LOAD_PROP_ANY,
                 'adjustable active power' : cload.LOAD_PROP_P_ADJUST}

str2prop_vargen = {'any' : cvargen.VARGEN_PROP_ANY}

str2prop_bat = {'any' : cbat.BAT_PROP_ANY}

str2prop = {'bus' : str2prop_bus,
            'branch' : str2prop_branch,
            'generator' : str2prop_gen,
            'shunt' : str2prop_shunt,
            'load' : str2prop_load,
            'variable generator' : str2prop_vargen,
            'battery' : str2prop_bat}

# Functions
str2func = {'voltage magnitude regularization' : cfunc.FUNC_TYPE_REG_VMAG,
            'voltage angle regularization' : cfunc.FUNC_TYPE_REG_VANG,
            'generator powers regularization' : cfunc.FUNC_TYPE_REG_PQ,
            'tap ratio regularization' : cfunc.FUNC_TYPE_REG_RATIO,
            'phase shift regularization' : cfunc.FUNC_TYPE_REG_PHASE,
            'susceptance regularization' : cfunc.FUNC_TYPE_REG_SUSC,
            'generation cost' : cfunc.FUNC_TYPE_GEN_COST,
            'sparse controls penalty' : cfunc.FUNC_TYPE_SP_CONTROLS,
            'soft voltage magnitude limits' : cfunc.FUNC_TYPE_SLIM_VMAG,
            'consumption utility' : cfunc.FUNC_TYPE_LOAD_UTIL,
            'net consumption cost' : cfunc.FUNC_TYPE_NETCON_COST}

func2str = dict([(v,k) for k,v in str2func.items()])

# Constraints
str2constr = {'AC power balance' : cconstr.CONSTR_TYPE_PF,
              'DC power balance' : cconstr.CONSTR_TYPE_DCPF,
              'linearized AC power balance' : cconstr.CONSTR_TYPE_LINPF,
              'variable fixing' : cconstr.CONSTR_TYPE_FIX,
              'variable nonlinear bounds' : cconstr.CONSTR_TYPE_BOUND,
              'generator active power participation' : cconstr.CONSTR_TYPE_PAR_GEN_P,
              'generator reactive power participation' : cconstr.CONSTR_TYPE_PAR_GEN_Q,
              'voltage regulation by generators' : cconstr.CONSTR_TYPE_REG_GEN,
              'voltage regulation by transformers' : cconstr.CONSTR_TYPE_REG_TRAN,
              'voltage regulation by shunts' : cconstr.CONSTR_TYPE_REG_SHUNT,
              'DC branch flow limits' : cconstr.CONSTR_TYPE_DC_FLOW_LIM,
              'AC branch flow limits' : cconstr.CONSTR_TYPE_AC_FLOW_LIM,
              'variable bounds' : cconstr.CONSTR_TYPE_LBOUND,
              'generator ramp limits' : cconstr.CONSTR_TYPE_GEN_RAMP}

constr2str = dict([(v,k) for k,v in str2constr.items()])
