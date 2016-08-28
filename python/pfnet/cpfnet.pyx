#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import numpy as np
cimport numpy as np

from functools import reduce

cimport cconstants
cimport cflags
cimport cobjs
cimport cvec
cimport cmat
cimport cgen
cimport cshunt
cimport cbus
cimport cbranch
cimport cload
cimport cvargen
cimport cbat
cimport cnet
cimport ccont
cimport cgraph
cimport cconstr
cimport cfunc
cimport cheur
cimport cprob

from scipy import misc
import tempfile

from scipy.sparse import coo_matrix, bmat

np.import_array()

# Options
#########

OPTION_SP_SCALARS = True # for scalar outputs when num periods is 1 (single period)

# Constants
###########

# Pi
PI = cconstants.PI

# Var values
CURRENT= cconstants.CURRENT
UPPER_LIMITS = cconstants.UPPER_LIMITS
LOWER_LIMITS = cconstants.LOWER_LIMITS

# Objects
OBJ_BUS = cobjs.OBJ_BUS
OBJ_GEN = cobjs.OBJ_GEN
OBJ_BRANCH = cobjs.OBJ_BRANCH
OBJ_SHUNT = cobjs.OBJ_SHUNT
OBJ_LOAD = cobjs.OBJ_LOAD
OBJ_VARGEN = cobjs.OBJ_VARGEN
OBJ_BAT = cobjs.OBJ_BAT
OBJ_UNKNOWN = cobjs.OBJ_UNKNOWN

# Flags
FLAG_VARS = cflags.FLAG_VARS
FLAG_FIXED = cflags.FLAG_FIXED
FLAG_BOUNDED = cflags.FLAG_BOUNDED
FLAG_SPARSE = cflags.FLAG_SPARSE

# C pointer
###########

cdef class CPtr:
    """
    C Pointer class.
    """
    
    cdef void* _c_ptr

cdef new_CPtr(void* ptr):
    cptr = CPtr()
    cptr._c_ptr = ptr
    return cptr

# Vector
########

cdef extern from "numpy/arrayobject.h":
     void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef Vector(cvec.Vec* v, owndata=False):
     cdef np.npy_intp shape[1]
     if v is not NULL:
         shape[0] = <np.npy_intp> cvec.VEC_get_size(v)
         arr = np.PyArray_SimpleNewFromData(1,shape,np.NPY_DOUBLE,cvec.VEC_get_data(v))
         if owndata:
             PyArray_ENABLEFLAGS(arr,np.NPY_OWNDATA)
         return arr
     else:
         return np.zeros(0)

# Matrix
########

cdef Matrix(cmat.Mat* m, owndata=False):
     cdef np.npy_intp shape[1]
     if m is not NULL:
         shape[0] = <np.npy_intp> cmat.MAT_get_nnz(m)
         size1 = cmat.MAT_get_size1(m)
         size2 = cmat.MAT_get_size2(m)
         row = np.PyArray_SimpleNewFromData(1,shape,np.NPY_INT,cmat.MAT_get_row_array(m))
         col = np.PyArray_SimpleNewFromData(1,shape,np.NPY_INT,cmat.MAT_get_col_array(m))
         data = np.PyArray_SimpleNewFromData(1,shape,np.NPY_DOUBLE,cmat.MAT_get_data_array(m))
         if owndata:
             PyArray_ENABLEFLAGS(row,np.NPY_OWNDATA)
             PyArray_ENABLEFLAGS(col,np.NPY_OWNDATA)
             PyArray_ENABLEFLAGS(data,np.NPY_OWNDATA)
         return coo_matrix((data,(row,col)),shape=(size1,size2))
     else:
         return coo_matrix(([],([],[])),shape=(0,0))

# Bus
#####

# Properties
BUS_PROP_ANY = cbus.BUS_PROP_ANY
BUS_PROP_SLACK = cbus.BUS_PROP_SLACK
BUS_PROP_REG_BY_GEN = cbus.BUS_PROP_REG_BY_GEN
BUS_PROP_REG_BY_TRAN = cbus.BUS_PROP_REG_BY_TRAN
BUS_PROP_REG_BY_SHUNT = cbus.BUS_PROP_REG_BY_SHUNT
BUS_PROP_NOT_REG_BY_GEN = cbus.BUS_PROP_NOT_REG_BY_GEN
BUS_PROP_NOT_SLACK = cbus.BUS_PROP_NOT_SLACK

# Variables
BUS_VAR_VMAG = cbus.BUS_VAR_VMAG
BUS_VAR_VANG = cbus.BUS_VAR_VANG
BUS_VAR_VDEV = cbus.BUS_VAR_VDEV
BUS_VAR_VVIO = cbus.BUS_VAR_VVIO

# Infinite
BUS_INF_V_MAG = cbus.BUS_INF_V_MAG
BUS_INF_V_ANG = cbus.BUS_INF_V_ANG

# Sensitivities
BUS_SENS_LARGEST = cbus.BUS_SENS_LARGEST
BUS_SENS_P_BALANCE = cbus.BUS_SENS_P_BALANCE
BUS_SENS_Q_BALANCE = cbus.BUS_SENS_Q_BALANCE
BUS_SENS_V_MAG_U_BOUND = cbus.BUS_SENS_V_MAG_U_BOUND
BUS_SENS_V_MAG_L_BOUND = cbus.BUS_SENS_V_MAG_L_BOUND
BUS_SENS_V_ANG_U_BOUND = cbus.BUS_SENS_V_ANG_U_BOUND
BUS_SENS_V_ANG_L_BOUND = cbus.BUS_SENS_V_ANG_L_BOUND
BUS_SENS_V_REG_BY_GEN = cbus.BUS_SENS_V_REG_BY_GEN
BUS_SENS_V_REG_BY_TRAN = cbus.BUS_SENS_V_REG_BY_TRAN
BUS_SENS_V_REG_BY_SHUNT = cbus.BUS_SENS_V_REG_BY_SHUNT
    
# Mismatches
BUS_MIS_LARGEST = cbus.BUS_MIS_LARGEST
BUS_MIS_ACTIVE = cbus.BUS_MIS_ACTIVE
BUS_MIS_REACTIVE = cbus.BUS_MIS_REACTIVE

class BusError(Exception):
    """
    Bus error exception.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Bus:
    """
    Bus class.
    """
    
    cdef cbus.Bus* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Bus class.
        
        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """
        
        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cbus.BUS_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_equal(self,other):
        """
        Determines whether bus is equal to given bus.

        Parameters
        ----------
        other : :class:`Bus <pfnet.Bus>`
        """

        cdef Bus b_other

        if not isinstance(other,Bus):
            return False

        b_other = other

        return cbus.BUS_is_equal(self._c_ptr,b_other._c_ptr)

    def is_slack(self):
        """
        Determines whether the bus is a slack bus.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_slack(self._c_ptr)

    def is_regulated_by_gen(self):
        """
        Determines whether the bus is regulated by a generator.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_gen(self._c_ptr)

    def is_regulated_by_tran(self):
        """
        Determines whether the bus is regulated by a transformer.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_tran(self._c_ptr)

    def is_regulated_by_shunt(self):
        """
        Determines whether the bus is regulated by a shunt device.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_is_regulated_by_shunt(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """
        Determines whether the bus has the flags associated with 
        certain quantities set.

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_bus_var`)

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbus.BUS_has_flags(self._c_ptr,fmask,vmask)

    def get_largest_sens(self,t=0):
        """
        Gets the bus sensitivity of largest absolute value.
       
        Parameters
        ----------
        t : int (time period)
 
        Returns
        -------
        sens : float
        """

        return cbus.BUS_get_largest_sens(self._c_ptr,t)

    def get_largest_sens_type(self,t=0):
        """
        Gets the type of bus sensitivity of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        type : int
        """

        return cbus.BUS_get_largest_sens_type(self._c_ptr,t)

    def get_largest_mis(self,t=0):
        """
        Gets the bus power mismatch of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        mis : float
        """

        return cbus.BUS_get_largest_mis(self._c_ptr,t)

    def get_largest_mis_type(self,t=0):
        """
        Gets the type of bus power mismatch of largest absolute value.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        type : int
        """

        return cbus.BUS_get_largest_mis_type(self._c_ptr,t)

    def get_quantity(self,type,t=0):
        """
        Gets the bus quantity of the given type. 

        Parameters
        ----------
        type : int (:ref:`ref_bus_sens`:, :ref:`ref_bus_mis`)
        t : int (time period)

        Returns
        -------
        value : float
        """

        return cbus.BUS_get_quantity(self._c_ptr,type,t)

    def get_total_gen_P(self,t=0):
        """
        Gets the total active power injected by generators
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        P : float
        """

        return cbus.BUS_get_total_gen_P(self._c_ptr,t)

    def get_total_gen_Q(self,t=0):
        """
        Gets the total reactive power injected by generators
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        Q : float
        """

        return cbus.BUS_get_total_gen_Q(self._c_ptr,t)

    def get_total_gen_Q_max(self):
        """ 
        Gets the largest total reactive power that can be 
        injected by generators connected to this bus.

        Returns
        -------
        Q_max : float
        """

        return cbus.BUS_get_total_gen_Q_max(self._c_ptr)

    def get_total_gen_Q_min(self):
        """
        Gets the smallest total reactive power that can be 
        injected by generators connected to this bus.

        Returns
        -------
        Q_min : float
        """

        return cbus.BUS_get_total_gen_Q_min(self._c_ptr)

    def get_total_load_P(self,t=0):
        """ 
        Gets the total active power consumed by loads
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        P : float
        """

        return cbus.BUS_get_total_load_P(self._c_ptr,t)

    def get_total_load_Q(self,t=0):
        """
        Gets the total reactive power consumed by loads
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        Q : float
        """

        return cbus.BUS_get_total_load_Q(self._c_ptr,t)

    def get_total_shunt_g(self):
        """ 
        Gets the combined conductance of shunt devices 
        connected to this bus.

        Returns
        -------
        g : float
        """

        return cbus.BUS_get_total_shunt_g(self._c_ptr)

    def get_total_shunt_b(self,t=0):
        """
        Gets the combined susceptance of shunt devices 
        connected to this bus.

        Parameters
        ----------
        t : int (time period)

        Returns
        -------
        b : float
        """

        return cbus.BUS_get_total_shunt_b(self._c_ptr,t)
    
    def show(self,t=0):
        """
        Shows bus properties.

        Parameters
        ----------
        t : int (time period)
        """
        cbus.BUS_show(self._c_ptr,t)

    def __richcmp__(self,other,op):
        """
        Compares two buses.

        Parameters
        ----------
        other : Bus
        op : comparison type

        Returns
        -------
        flag : {``True``, ``False``}
        """

        if op == 2:
            return self.is_equal(other)
        elif op == 3:
            return not self.is_equal(other)
        else:
            return False

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cbus.BUS_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cbus.BUS_get_obj_type(self._c_ptr)

    property index:
        """ Bus index (int). """
        def __get__(self): return cbus.BUS_get_index(self._c_ptr)

    property index_v_mag: 
        """ Index of voltage magnitude variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_v_mag(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_v_ang:
        """ Index of voltage angle variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_v_ang(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        
    property index_y:
        """ Index of voltage magnitude positive deviation variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_y(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)        

    property index_z:
        """ Index of voltage magnitude negative deviation variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_z(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)        

    property index_vl:
        """ Index of voltage low limit violation variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_vl(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_vh:
        """ Index of voltage high limit violation variable (int or array). """
        def __get__(self): 
            r = [cbus.BUS_get_index_vh(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_P:
        """ Index of bus active power mismatch (int). """
        def __get__(self): return cbus.BUS_get_index_P(self._c_ptr)

    property index_Q:
        """ Index for bus reactive power mismatch (int). """
        def __get__(self): return cbus.BUS_get_index_Q(self._c_ptr)

    property price:
        """ Bus energy price (float or array) ($ / (hr p.u.)). """
        def __get__(self): 
            r = [cbus.BUS_get_price(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,p):
            cdef int t
            cdef np.ndarray par = np.array(p).flatten()
            for t in range(np.minimum(par.size,self.num_periods)):
                cbus.BUS_set_price(self._c_ptr,par[t],t)

    property number:
        """ Bus number (int). """
        def __get__(self): return cbus.BUS_get_number(self._c_ptr)

    property name:
        """ Bus name (sting). """
        def __get__(self): return cbus.BUS_get_name(self._c_ptr).decode('UTF-8')
        def __set__(self,name): 
            name = name.encode('UTF-8')
            cbus.BUS_set_name(self._c_ptr,name)

    property degree:
        """ Bus degree (number of incident branches) (float). """
        def __get__(self): return cbus.BUS_get_degree(self._c_ptr)    

    property v_mag:
        """ Bus volatge magnitude (p.u. bus base kv) (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_v_mag(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,v):
            cdef int t
            cdef np.ndarray var = np.array(v).flatten()
            for t in range(np.minimum(var.size,self.num_periods)):
                cbus.BUS_set_v_mag(self._c_ptr,var[t],t)

    property v_ang:
        """ Bus voltage angle (radians) (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_v_ang(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,v):
            cdef int t
            cdef np.ndarray var = np.array(v).flatten()
            for t in range(np.minimum(var.size,self.num_periods)):
                cbus.BUS_set_v_ang(self._c_ptr,var[t],t)

    property v_set:
        """ Bus voltage set point (p.u. bus base kv) (float or array). Equals one if bus is not regulated by a generator. """
        def __get__(self): 
            r = [cbus.BUS_get_v_set(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property v_max:
        """ Bus volatge upper bound (p.u. bus base kv) (float). """
        def __get__(self): return cbus.BUS_get_v_max(self._c_ptr)

    property v_min:
        """ Bus voltage lower bound (p.u. bus base kv) (float). """
        def __get__(self): return cbus.BUS_get_v_min(self._c_ptr)

    property P_mis:
        """ Bus active power mismatch (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_P_mis(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property Q_mis:
        """ Bus reactive power mismatch (p.u. system base MVA) (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_Q_mis(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_P_balance:
        """ Objective function sensitivity with respect to bus active power balance (float or array). """
        def __get__(self):
            r = [cbus.BUS_get_sens_P_balance(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_Q_balance:
        """ Objective function sensitivity with respect to bus reactive power balance (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_Q_balance(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_mag_u_bound:
        """ Objective function sensitivity with respect to voltage magnitude upper bound (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_mag_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_mag_l_bound:
        """ Objective function sensitivity with respect to voltage magnitude lower bound (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_mag_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_ang_u_bound:
        """ Objective function sensitivity with respect to voltage angle upper bound (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_ang_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_ang_l_bound:
        """ Objective function sensitivity with respect to voltage angle lower bound (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_ang_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_reg_by_gen:
        """ Objective function sensitivity with respect to bus voltage regulation by generators (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_reg_by_gen(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_reg_by_tran:
        """ Objective function sensitivity with respect to bus voltage regulation by transformers (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_reg_by_tran(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_v_reg_by_shunt:
        """ Objective function sensitivity with respect to bus voltage regulation by shunts (float or array). """
        def __get__(self): 
            r = [cbus.BUS_get_sens_v_reg_by_shunt(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property generators:
        """ List of :class:`generators <pfnet.Generator>` connected to this bus (list). """
        def __get__(self):
            gens = []
            cdef cgen.Gen* g = cbus.BUS_get_gen(self._c_ptr)
            while g is not NULL:
                gens.append(new_Generator(g))
                g = cgen.GEN_get_next(g)
            return gens

    property reg_generators: 
        """ List of :class:`generators <pfnet.Generator>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_gens = []
            cdef cgen.Gen* g = cbus.BUS_get_reg_gen(self._c_ptr)
            while g is not NULL:
                reg_gens.append(new_Generator(g))
                g = cgen.GEN_get_reg_next(g)
            return reg_gens

    property reg_trans:
        """ List of :class:`tap-changing transformers <pfnet.Branch>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_trans = []
            cdef cbranch.Branch* br = cbus.BUS_get_reg_tran(self._c_ptr)
            while br is not NULL:
                reg_trans.append(new_Branch(br))
                br = cbranch.BRANCH_get_reg_next(br)
            return reg_trans

    property reg_shunts:
        """ List of :class:`switched shunt devices <pfnet.Shunt>` regulating the voltage magnitude of this bus (list). """
        def __get__(self):
            reg_shunts = []
            cdef cshunt.Shunt* s = cbus.BUS_get_reg_shunt(self._c_ptr)
            while s is not NULL:
                reg_shunts.append(new_Shunt(s))
                s = cshunt.SHUNT_get_reg_next(s)
            return reg_shunts

    property branches_from:
        """ List of :class:`branches <pfnet.Branch>` that have this bus on the "from" side (list). """
        def __get__(self):
            branches = []
            cdef cbranch.Branch* br = cbus.BUS_get_branch_from(self._c_ptr)
            while br is not NULL:
                branches.append(new_Branch(br))
                br = cbranch.BRANCH_get_from_next(br)
            return branches

    property branches_to:
        """ List of :class:`branches <pfnet.Branch>` that have this bus on the "to" side (list). """
        def __get__(self):
            branches = []
            cdef cbranch.Branch* br = cbus.BUS_get_branch_to(self._c_ptr)
            while br is not NULL:
                branches.append(new_Branch(br))
                br = cbranch.BRANCH_get_to_next(br)
            return branches

    property branches:
        """ List of :class:`branches <pfnet.Branch>` incident on this bus (list). """
        def __get__(self):
            return self.branches_from+self.branches_to

    property loads:
        """ List of :class:`loads <pfnet.Load>` connected to this bus (list). """
        def __get__(self):
            loads = []
            cdef cload.Load* l = cbus.BUS_get_load(self._c_ptr)
            while l is not NULL:
                loads.append(new_Load(l))
                l = cload.LOAD_get_next(l)
            return loads

    property var_generators:
        """ List of :class:`variable generators <pfnet.VarGenerator>` connected to this bus (list). """
        def __get__(self):
            vargens = []
            cdef cvargen.Vargen* g = cbus.BUS_get_vargen(self._c_ptr)
            while g is not NULL:
                vargens.append(new_VarGenerator(g))
                g = cvargen.VARGEN_get_next(g)
            return vargens

    property batteries:
        """ List of :class:`batteries <pfnet.Battery>` connected to this bus (list). """
        def __get__(self):
            bats = []
            cdef cbat.Bat* b = cbus.BUS_get_bat(self._c_ptr)
            while b is not NULL:
                bats.append(new_Battery(b))
                b = cbat.BAT_get_next(b)
            return bats

cdef new_Bus(cbus.Bus* b):
    if b is not NULL:
        bus = Bus(alloc=False)
        bus._c_ptr = b
        return bus
    else:
        raise BusError('no bus data')

# Branch
########

# Properties
BRANCH_PROP_ANY = cbranch.BRANCH_PROP_ANY
BRANCH_PROP_TAP_CHANGER = cbranch.BRANCH_PROP_TAP_CHANGER
BRANCH_PROP_TAP_CHANGER_V = cbranch.BRANCH_PROP_TAP_CHANGER_V
BRANCH_PROP_TAP_CHANGER_Q = cbranch.BRANCH_PROP_TAP_CHANGER_Q
BRANCH_PROP_PHASE_SHIFTER = cbranch.BRANCH_PROP_PHASE_SHIFTER
BRANCH_PROP_NOT_OUT = cbranch.BRANCH_PROP_NOT_OUT

# Variables
BRANCH_VAR_RATIO = cbranch.BRANCH_VAR_RATIO
BRANCH_VAR_RATIO_DEV = cbranch.BRANCH_VAR_RATIO_DEV
BRANCH_VAR_PHASE = cbranch.BRANCH_VAR_PHASE

# Infinite
BRANCH_INF_RATIO = cbranch.BRANCH_INF_RATIO
BRANCH_INF_FLOW = cbranch.BRANCH_INF_FLOW

class BranchError(Exception):
    """
    Branch error exception.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Branch:
    """
    Branch class.
    """

    cdef cbranch.Branch* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Branch class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cbranch.BRANCH_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def has_pos_ratio_v_sens(self):
        """
        Determines whether tap-changing transformer has positive
        sensitivity between tap ratio and controlled bus voltage magnitude.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_has_pos_ratio_v_sens(self._c_ptr)

    def is_equal(self,other):
        """
        Determines whether branch is equal to given branch.

        Parameters
        ----------
        other : :class:`Branch <pfnet.Branch>`
        """

        cdef Branch b_other

        if not isinstance(other,Branch):
            return False

        b_other = other

        return cbranch.BRANCH_is_equal(self._c_ptr,b_other._c_ptr)

    def __richcmp__(self,other,op):
        """
        Compares two branches.

        Parameters
        ----------
        other : Branch
        op : comparison type

        Returns
        -------
        flag : {``True``, ``False``}
        """

        if op == 2:
            return self.is_equal(other)
        elif op == 3:
            return not self.is_equal(other)
        else:
            return False

    def is_on_outage(self):
        """
        Determines whether branch in on outage.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_on_outage(self._c_ptr)

    def is_fixed_tran(self):
        """
        Determines whether branch is fixed transformer.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_fixed_tran(self._c_ptr)

    def is_line(self):
        """
        Determines whether branch is transmission line.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_line(self._c_ptr)

    def is_phase_shifter(self):
        """
        Determines whether branch is phase shifter.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_phase_shifter(self._c_ptr)

    def is_tap_changer(self):
        """
        Determines whether branch is tap-changing transformer.

        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cbranch.BRANCH_is_tap_changer(self._c_ptr)

    def is_tap_changer_v(self):
        """
        Determines whether branch is tap-changing transformer
        that regulates bus voltage magnitude.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_tap_changer_v(self._c_ptr)

    def is_tap_changer_Q(self):
        """
        Determines whether branch is tap-changing transformer
        that regulates reactive power flow.

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_is_tap_changer_Q(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """
        Determines whether the branch has the flags associated with
        specific quantities set.

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_branch_var`)

        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbranch.BRANCH_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cbranch.BRANCH_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cbranch.BRANCH_get_obj_type(self._c_ptr)

    property index:
        """ Branch index (int). """
        def __get__(self): return cbranch.BRANCH_get_index(self._c_ptr)

    property index_ratio:
        """ Index of transformer tap ratio variable (int or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_index_ratio(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_ratio_y:
        """ Index of transformer tap ratio positive deviation variable (int or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_index_ratio_y(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_ratio_z:
        """ Index of transformer tap ratio negative deviation variable (int or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_index_ratio_z(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_phase:
        """ Index of transformer phase shift variable (int or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_index_phase(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property ratio:
        """ Transformer tap ratio (float or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_ratio(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property ratio_max:
        """ Transformer tap ratio upper limit (float). """
        def __get__(self): return cbranch.BRANCH_get_ratio_max(self._c_ptr)
        def __set__(self,value): cbranch.BRANCH_set_ratio_max(self._c_ptr,value)

    property ratio_min:
        """ Transformer tap ratio lower limit (float). """
        def __get__(self): return cbranch.BRANCH_get_ratio_min(self._c_ptr)
        def __set__(self,value): cbranch.BRANCH_set_ratio_min(self._c_ptr,value)

    property bus_from:
        """ :class:`Bus <pfnet.Bus>` connected to the "from" side. """
        def __get__(self): return new_Bus(cbranch.BRANCH_get_bus_from(self._c_ptr))
   
    property bus_to:
        """ :class:`Bus <pfnet.Bus>` connected to the "to" side. """
        def __get__(self): return new_Bus(cbranch.BRANCH_get_bus_to(self._c_ptr))

    property reg_bus:
        """ :class:`Bus <pfnet.Bus>` whose voltage is regulated by this tap-changing transformer. """
        def __get__(self): return new_Bus(cbranch.BRANCH_get_reg_bus(self._c_ptr))

    property b:
        """ Branch series susceptance (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_b(self._c_ptr)

    property b_from:
        """ Branch shunt susceptance at the "from" side (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_b_from(self._c_ptr)

    property b_to:
        """ Branch shunt susceptance at the "to" side (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_b_to(self._c_ptr)

    property g:
        """ Branch series conductance (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_g(self._c_ptr)

    property g_from:
        """ Branch shunt conductance at the "from" side (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_g_from(self._c_ptr)

    property g_to:
        """ Branch shunt conductance at the "to" side (p.u.) (float). """
        def __get__(self): return cbranch.BRANCH_get_g_to(self._c_ptr)

    property phase:
        """ Transformer phase shift (radians) (float or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_phase(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property phase_max:
        """ Transformer phase shift upper limit (radians) (float). """
        def __get__(self): return cbranch.BRANCH_get_phase_max(self._c_ptr)

    property phase_min:
        """ Transformer phase shift lower limit (radians) (float). """
        def __get__(self): return cbranch.BRANCH_get_phase_min(self._c_ptr)

    property ratingA:
        """ Branch thermal rating A (p.u. system base power) (float). """
        def __get__(self): return cbranch.BRANCH_get_ratingA(self._c_ptr)
        def __set__(self,r): cbranch.BRANCH_set_ratingA(self._c_ptr,r)

    property ratingB:
        """ Branch thermal rating B (p.u. system base power) (float). """
        def __get__(self): return cbranch.BRANCH_get_ratingB(self._c_ptr)
        def __set__(self,r): cbranch.BRANCH_set_ratingB(self._c_ptr,r)

    property ratingC:
        """ Branch thermal rating C (p.u. system base power) (float). """
        def __get__(self): return cbranch.BRANCH_get_ratingC(self._c_ptr)
        def __set__(self,r): cbranch.BRANCH_set_ratingC(self._c_ptr,r)

    property P_flow_DC:
        """ Active power flow (DC approx.) from bus "from" to bus "to" (float or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_P_flow_DC(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power flow upper bound (float or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_sens_P_u_bound(self._c_ptr,t) for t in range(self.num_periods)] 
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
                
    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power flow lower bound (float or array). """
        def __get__(self): 
            r = [cbranch.BRANCH_get_sens_P_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property outage:
        """ Flag that indicates whehter branch is on outage. """
        def __get__(self): return cbranch.BRANCH_is_on_outage(self._c_ptr)

cdef new_Branch(cbranch.Branch* b):
    if b is not NULL:
        branch = Branch(alloc=False)
        branch._c_ptr = b
        return branch
    else:
        raise BranchError('no branch data')

# Generator
###########

# Properties
GEN_PROP_ANY = cgen.GEN_PROP_ANY
GEN_PROP_SLACK = cgen.GEN_PROP_SLACK
GEN_PROP_REG = cgen.GEN_PROP_REG
GEN_PROP_NOT_REG = cgen.GEN_PROP_NOT_REG
GEN_PROP_NOT_SLACK = cgen.GEN_PROP_NOT_SLACK
GEN_PROP_NOT_OUT = cgen.GEN_PROP_NOT_OUT
GEN_PROP_P_ADJUST = cgen.GEN_PROP_P_ADJUST

# Variables
GEN_VAR_P = cgen.GEN_VAR_P
GEN_VAR_Q = cgen.GEN_VAR_Q

# Infinity
GEN_INF_P = cgen.GEN_INF_P
GEN_INF_Q = cgen.GEN_INF_Q

class GeneratorError(Exception):
    """ 
    Generator error exception.
    """
    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Generator:
    """
    Generator class.
    """

    cdef cgen.Gen* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Generator class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cgen.GEN_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_equal(self,other):
        """
        Determines whether generator is equal to given generator.

        Parameters
        ----------
        other : :class:`Generator <pfnet.Generator>`
        """

        cdef Generator g_other

        if not isinstance(other,Generator):
            return False

        g_other = other

        return cgen.GEN_is_equal(self._c_ptr,g_other._c_ptr)

    def __richcmp__(self,other,op):
        """
        Compares two generators.

        Parameters
        ----------
        other : Generator
        op : comparison type

        Returns
        -------
        flag : {``True``, ``False``}
        """

        if op == 2:
            return self.is_equal(other)
        elif op == 3:
            return not self.is_equal(other)
        else:
            return False

    def is_on_outage(self):
        """
        Determines whether generator in on outage.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cgen.GEN_is_on_outage(self._c_ptr)

    def is_slack(self):
        """ 
        Determines whether generator is slack.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cgen.GEN_is_slack(self._c_ptr)

    def is_regulator(self):
        """ 
        Determines whether generator provides voltage regulation.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cgen.GEN_is_regulator(self._c_ptr)

    def is_P_adjustable(self):
        """ 
        Determines whether generator has adjustable active power.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cgen.GEN_is_P_adjustable(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """ 
        Determines whether the generator has the flags associated with
        certain quantities set. 

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_gen_var`)
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cgen.GEN_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cgen.GEN_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cgen.GEN_get_obj_type(self._c_ptr)

    property index:
        """ Generator index (int). """
        def __get__(self): return cgen.GEN_get_index(self._c_ptr)
        
    property index_P:
        """ Index of generator active power variable (int or array). """
        def __get__(self): 
            r = [cgen.GEN_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_Q:
        """ Index of generator reactive power variable (int or array). """
        def __get__(self): 
            r = [cgen.GEN_get_index_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which generator is connected. """
        def __get__(self): return new_Bus(cgen.GEN_get_bus(self._c_ptr))

    property reg_bus:
        """ :class:`Bus <pfnet.Bus>` whose voltage is regulated by this generator. """
        def __get__(self): return new_Bus(cgen.GEN_get_reg_bus(self._c_ptr))

    property P:
        """ Generator active power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cgen.GEN_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,P): 
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cgen.GEN_set_P(self._c_ptr,Par[t],t)

    property P_max:
        """ Generator active power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_P_max(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_P_max(self._c_ptr,P)

    property P_min:
        """ Generator active power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_P_min(self._c_ptr)
        def __set__(self,P): cgen.GEN_set_P_min(self._c_ptr,P)
            
    property Q:
        """ Generator reactive power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cgen.GEN_get_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,Q): 
            cdef int t
            cdef np.ndarray Qar = np.array(Q).flatten()
            for t in range(np.minimum(Qar.size,self.num_periods)):
                cgen.GEN_set_Q(self._c_ptr,Qar[t],t)

    property Q_max:
        """ Generator reactive power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_Q_max(self._c_ptr)

    property Q_min:
        """ Generator reactive power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cgen.GEN_get_Q_min(self._c_ptr)

    property P_cost:
        """ Active power generation cost ($/hr) (float or array). """
        def __get__(self): 
            r = [cgen.GEN_get_P_cost(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
            
    property cost_coeff_Q0:
        """ Coefficient for genertion cost function (constant term, units of $/hr). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q0(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q0(self._c_ptr,c)

    property cost_coeff_Q1:
        """ Coefficient for genertion cost function (linear term, units of $/(hr p.u.)). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q1(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q1(self._c_ptr,c)

    property cost_coeff_Q2:
        """ Coefficient for genertion cost function (quadratic term, units of $/(hr p.u.^2)). """
        def __get__(self): return cgen.GEN_get_cost_coeff_Q2(self._c_ptr)
        def __set__(self,c): cgen.GEN_set_cost_coeff_Q2(self._c_ptr,c)

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power upper bound (float or array). """
        def __get__(self): 
            r = [cgen.GEN_get_sens_P_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power lower bound (float or array). """
        def __get__(self): 
            r = [cgen.GEN_get_sens_P_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property outage:
        """ Flag that indicates whehter generator is on outage. """
        def __get__(self): return cgen.GEN_is_on_outage(self._c_ptr)

cdef new_Generator(cgen.Gen* g):
    if g is not NULL:
        gen = Generator(alloc=False)
        gen._c_ptr = g
        return gen
    else:
        raise GeneratorError('no gen data')

# Shunt
#######

# Properties
SHUNT_PROP_ANY = cshunt.SHUNT_PROP_ANY
SHUNT_PROP_SWITCHED_V = cshunt.SHUNT_PROP_SWITCHED_V

# Variables
SHUNT_VAR_SUSC = cshunt.SHUNT_VAR_SUSC
SHUNT_VAR_SUSC_DEV = cshunt.SHUNT_VAR_SUSC_DEV

# Infinite
SHUNT_INF_SUSC = cshunt.SHUNT_INF_SUSC

class ShuntError(Exception):
    """
    Shunt error exception.
    """
    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Shunt:
    """
    Shunt class.
    """

    cdef cshunt.Shunt* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Shunt class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cshunt.SHUNT_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_fixed(self):
        """
        Determines whether the shunt device is fixed (as opposed to switched).

        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cshunt.SHUNT_is_fixed(self._c_ptr)

    def is_switched_v(self):
        """
        Determines whether the shunt is switchable and regulates 
        bus voltage magnitude.

        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cshunt.SHUNT_is_switched_v(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """
        Determines whether the shunt devices has flags associated with 
        certain quantities set.
        
        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_bus_var`)

        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cshunt.SHUNT_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cshunt.SHUNT_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cshunt.SHUNT_get_obj_type(self._c_ptr)

    property index:
        """ Shunt index (int). """
        def __get__(self): return cshunt.SHUNT_get_index(self._c_ptr)    

    property index_b:
        """ Index of shunt susceptance variable (int or array). """
        def __get__(self): 
            r = [cshunt.SHUNT_get_index_b(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_y:
        """ Index of shunt susceptance positive deviation variable (int or array). """
        def __get__(self): 
            r = [cshunt.SHUNT_get_index_y(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_z:
        """ Index of shunt susceptance negative deviation variable (int or array). """
        def __get__(self): 
            r = [cshunt.SHUNT_get_index_z(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which the shunt devices is connected. """
        def __get__(self): return new_Bus(cshunt.SHUNT_get_bus(self._c_ptr))

    property reg_bus:
        """ :class:`Bus <pfnet.Bus>` whose voltage magnitude is regulated by this shunt device. """
        def __get__(self): return new_Bus(cshunt.SHUNT_get_reg_bus(self._c_ptr))

    property g:
        """ Shunt conductance (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_g(self._c_ptr)
            
    property b:
        """ Shunt susceptance (p.u.) (float or array). """
        def __get__(self): 
            r = [cshunt.SHUNT_get_b(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property b_max:
        """ Shunt susceptance upper limit (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_b_max(self._c_ptr)
        def __set__(self,value): cshunt.SHUNT_set_b_max(self._c_ptr,value)

    property b_min:
        """ Shunt susceptance lower limit (p.u.) (float). """
        def __get__(self): return cshunt.SHUNT_get_b_min(self._c_ptr)
        def __set__(self,value): cshunt.SHUNT_set_b_min(self._c_ptr,value)                

cdef new_Shunt(cshunt.Shunt* s):
    if s is not NULL:
        shunt = Shunt(alloc=False)
        shunt._c_ptr = s
        return shunt
    else:
        raise ShuntError('no shunt data')

# Load
######

# Properties
LOAD_PROP_ANY = cload.LOAD_PROP_ANY
LOAD_PROP_P_ADJUST = cload.LOAD_PROP_P_ADJUST

# Variables
LOAD_VAR_P = cload.LOAD_VAR_P

# Infinity
LOAD_INF_P = cload.LOAD_INF_P

class LoadError(Exception):
    """ 
    Load error exception.
    """
    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Load:
    """
    Load class.
    """

    cdef cload.Load* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Load class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_ptr = cload.LOAD_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def is_P_adjustable(self):
        """ 
        Determines whether the load has adjustable active power.
        
        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cload.LOAD_is_P_adjustable(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """ 
        Determines whether the load has the flags associated with
        certain quantities set. 

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_load_var`)
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cload.LOAD_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cload.LOAD_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cload.LOAD_get_obj_type(self._c_ptr)

    property index:
        """ Load index (int). """
        def __get__(self): return cload.LOAD_get_index(self._c_ptr)
       
    property index_P:
        """ Index of load active power variable (int or array). """
        def __get__(self): 
            r = [cload.LOAD_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
 
    property bus:
        """ :class:`Bus <pfnet.Bus>` to which load is connected. """
        def __get__(self): return new_Bus(cload.LOAD_get_bus(self._c_ptr))

    property P:
        """ Load active power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cload.LOAD_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,P): 
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cload.LOAD_set_P(self._c_ptr,Par[t],t)
       
    property P_max:
        """ Load active power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cload.LOAD_get_P_max(self._c_ptr)
        def __set__(self,P): cload.LOAD_set_P_max(self._c_ptr,P)

    property P_min:
        """ Load active power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cload.LOAD_get_P_min(self._c_ptr)
        def __set__(self,P): cload.LOAD_set_P_min(self._c_ptr,P)
     
    property Q:
        """ Load reactive power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cload.LOAD_get_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,Q): 
            cdef int t
            cdef np.ndarray Qar = np.array(Q).flatten()
            for t in range(np.minimum(Qar.size,self.num_periods)):
                cload.LOAD_set_Q(self._c_ptr,Qar[t],t)

    property P_util:
        """ Active power load utility ($/hr) (float or array). """
        def __get__(self): 
            r = [cload.LOAD_get_P_util(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property util_coeff_Q0:
        """ Coefficient for consumption utility function (constant term, units of $/hr). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q0(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q0(self._c_ptr,c)

    property util_coeff_Q1:
        """ Coefficient for consumption utility function (linear term, units of $/(hr p.u.)). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q1(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q1(self._c_ptr,c)

    property util_coeff_Q2:
        """ Coefficient for consumption utility function (quadratic term, units of $/(hr p.u.^2)). """
        def __get__(self): return cload.LOAD_get_util_coeff_Q2(self._c_ptr)
        def __set__(self,c): cload.LOAD_set_util_coeff_Q2(self._c_ptr,c)

    property sens_P_u_bound:
        """ Objective function sensitivity with respect to active power upper bound (float or array). """
        def __get__(self): 
            r = [cload.LOAD_get_sens_P_u_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property sens_P_l_bound:
        """ Objective function sensitivity with respect to active power lower bound (float or array). """
        def __get__(self): 
            r = [cload.LOAD_get_sens_P_l_bound(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

cdef new_Load(cload.Load* l):
    if l is not NULL:
        load = Load(alloc=False)
        load._c_ptr = l
        return load
    else:
        raise LoadError('no load data')

# Variable Generator
####################

# Properties
VARGEN_PROP_ANY = cvargen.VARGEN_PROP_ANY

# Variables
VARGEN_VAR_P = cvargen.VARGEN_VAR_P
VARGEN_VAR_Q = cvargen.VARGEN_VAR_Q

# Infinity
VARGEN_INF_P = cvargen.VARGEN_INF_P
VARGEN_INF_Q = cvargen.VARGEN_INF_Q

class VarGeneratorError(Exception):
    """ 
    Variable generator error exception.
    """
    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class VarGenerator:
    """
    Variable generator class.
    """

    cdef cvargen.Vargen* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Variable generator class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):
        
        if alloc:
            self._c_ptr = cvargen.VARGEN_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """ 
        Determines whether the variable generator has the flags associated with
        certain quantities set. 

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_vargen_var`)
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cvargen.VARGEN_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cvargen.VARGEN_get_num_periods(self._c_ptr)

    property name:
        """ Variable generator name (string). """
        def __get__(self): return cvargen.VARGEN_get_name(self._c_ptr).decode('UTF-8')
        def __set__(self,name): 
            name = name.encode('UTF-8')
            cvargen.VARGEN_set_name(self._c_ptr,name)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cvargen.VARGEN_get_obj_type(self._c_ptr)

    property index:
        """ Variable generator index (int). """
        def __get__(self): return cvargen.VARGEN_get_index(self._c_ptr)
        
    property index_P:
        """ Index of variable generator active power variable (int or array). """
        def __get__(self): 
            r = [cvargen.VARGEN_get_index_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_Q:
        """ Index of variable generator reactive power variable (int or array). """
        def __get__(self): 
            r = [cvargen.VARGEN_get_index_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which variable generator is connected. """
        def __get__(self): return new_Bus(cvargen.VARGEN_get_bus(self._c_ptr))

    property P:
        """ Variable generator active power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cvargen.VARGEN_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,P): 
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cvargen.VARGEN_set_P(self._c_ptr,Par[t],t)

    property P_max:
        """ Variable generator active power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cvargen.VARGEN_get_P_max(self._c_ptr)
        def __set__(self,P): cvargen.VARGEN_set_P_max(self._c_ptr,P)

    property P_min:
        """ Variable generator active power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cvargen.VARGEN_get_P_min(self._c_ptr)
        def __set__(self,P): cvargen.VARGEN_set_P_min(self._c_ptr,P)

    property P_std:
        """ Variable generator active power standard deviation (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cvargen.VARGEN_get_P_std(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,P): 
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cvargen.VARGEN_set_P_std(self._c_ptr,Par[t],t)

    property Q:
        """ Variable generator reactive power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cvargen.VARGEN_get_Q(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,Q): 
            cdef int t
            cdef np.ndarray Qar = np.array(Q).flatten()
            for t in range(np.minimum(Qar.size,self.num_periods)):
                cvargen.VARGEN_set_Q(self._c_ptr,Qar[t],t)

    property Q_max:
        """ Variable generator maximum reactive power (p.u. system base MVA) (float). """
        def __get__(self): return cvargen.VARGEN_get_Q_max(self._c_ptr)
        def __set__(self,Q): cvargen.VARGEN_set_Q_max(self._c_ptr,Q)

    property Q_min:
        """ Variable generator minimum reactive power (p.u. system base MVA) (float). """
        def __get__(self): return cvargen.VARGEN_get_Q_min(self._c_ptr)
        def __set__(self,Q): cvargen.VARGEN_set_Q_min(self._c_ptr,Q)

cdef new_VarGenerator(cvargen.Vargen* g):
    if g is not NULL:
        gen = VarGenerator(alloc=False)
        gen._c_ptr = g
        return gen
    else:
        raise VarGeneratorError('no vargen data')

# Battery
#########

# Properties
BAT_PROP_ANY = cbat.BAT_PROP_ANY

# Variables
BAT_VAR_P = cbat.BAT_VAR_P
BAT_VAR_E = cbat.BAT_VAR_E

# Infinity
BAT_INF_P = cbat.BAT_INF_P
BAT_INF_E = cbat.BAT_INF_E

class BatteryError(Exception):
    """ 
    Battery error exception.
    """
    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Battery:
    """
    Battery class.
    """

    cdef cbat.Bat* _c_ptr

    def __init__(self,num_periods=1,alloc=True):
        """
        Battery class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass

    def __cinit__(self,num_periods=1,alloc=True):
        
        if alloc:
            self._c_ptr = cbat.BAT_new(num_periods)
        else:
            self._c_ptr = NULL

    def _get_c_ptr(self):

        return new_CPtr(self._c_ptr)

    def has_flags(self,fmask,vmask):
        """ 
        Determines whether the battery has the flags associated with
        certain quantities set. 

        Parameters
        ----------
        fmask : int (:ref:`ref_net_flag`)
        vmask : int (:ref:`ref_bat_var`)
        
        Returns
        -------
        flag : {``True``, ``False``}
        """

        return cbat.BAT_has_flags(self._c_ptr,fmask,vmask)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cbat.BAT_get_num_periods(self._c_ptr)

    property obj_type:
        """ Object type (int). """
        def __get__(self): return cbat.BAT_get_obj_type(self._c_ptr)

    property index:
        """ Battery index (int). """
        def __get__(self): return cbat.BAT_get_index(self._c_ptr)
        
    property index_Pc:
        """ Index of battery charging power variable (int or array). """
        def __get__(self): 
            r = [cbat.BAT_get_index_Pc(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_Pd:
        """ Index of battery discharging power variable (int or array). """
        def __get__(self): 
            r = [cbat.BAT_get_index_Pd(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property index_E:
        """ Index of battery energy level variable (int or array). """
        def __get__(self): 
            r = [cbat.BAT_get_index_E(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus:
        """ :class:`Bus <pfnet.Bus>` to which battery is connected. """
        def __get__(self): return new_Bus(cbat.BAT_get_bus(self._c_ptr))

    property P:
        """ Battery charging power (p.u. system base MVA) (float or array). """
        def __get__(self): 
            r = [cbat.BAT_get_P(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,P): 
            cdef int t
            cdef np.ndarray Par = np.array(P).flatten()
            for t in range(np.minimum(Par.size,self.num_periods)):
                cbat.BAT_set_P(self._c_ptr,Par[t],t)

    property P_max:
        """ Battery charging power upper limit (p.u. system base MVA) (float). """
        def __get__(self): return cbat.BAT_get_P_max(self._c_ptr)
        def __set__(self,P): cbat.BAT_set_P_max(self._c_ptr,P)

    property P_min:
        """ Battery charging power lower limit (p.u. system base MVA) (float). """
        def __get__(self): return cbat.BAT_get_P_min(self._c_ptr)
        def __set__(self,P): cbat.BAT_set_P_min(self._c_ptr,P)

    property E:
        """ Battery energy level (p.u. system base MVA times time unit) (float or array). """
        def __get__(self): 
            r = [cbat.BAT_get_E(self._c_ptr,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
        def __set__(self,E): 
            cdef int t
            cdef np.ndarray Ear = np.array(E).flatten()
            for t in range(np.minimum(Ear.size,self.num_periods)):
                cbat.BAT_set_E(self._c_ptr,Ear[t],t)

    property E_max:
        """ Battery energy level upper limit (p.u. system base MVA times time unit) (float). """
        def __get__(self): return cbat.BAT_get_E_max(self._c_ptr)
        def __set__(self,E): cbat.BAT_set_E_max(self._c_ptr,E)

    property eta_c:
        """ Battery charging efficiency (unitless) (float). """
        def __get__(self): return cbat.BAT_get_eta_c(self._c_ptr)
        def __set__(self,eta_c): cbat.BAT_set_eta_c(self._c_ptr,eta_c)

    property eta_d:
        """ Battery discharging efficiency (unitless) (float). """
        def __get__(self): return cbat.BAT_get_eta_d(self._c_ptr)
        def __set__(self,eta_d): cbat.BAT_set_eta_d(self._c_ptr,eta_d)

cdef new_Battery(cbat.Bat* b):
    if b is not NULL:
        bat = Battery(alloc=False)
        bat._c_ptr = b
        return bat
    else:
        raise BatteryError('no battery data')

# Network
#########

class NetworkError(Exception):
    """
    Network error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Network:
    """
    Network class.
    """

    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self,num_periods=1,alloc=True):
        """
        Network class.

        Parameters
        ----------
        alloc : {``True``, ``False``}
        num_periods : int
        """

        pass
     
    def __cinit__(self,num_periods=1,alloc=True):

        if alloc:
            self._c_net = cnet.NET_new(num_periods)
        else:
            self._c_net = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees network C data structure. 
        """
        
        if self.alloc:
            cnet.NET_del(self._c_net)
            self._c_net = NULL

    def add_vargens(self,buses,penetration,uncertainty,corr_radius,corr_value):
        """
        Adds variable generators to the network.

        Parameters
        ----------

        buses : list of :class:`Buses <pfnet.Bus>`
        penetration : float
                      percentage
        uncertainty : float
                      percentage
        corr_radius : int
                      number of branches
        corr_value : float
                     correlation coefficient
        """
        
        cdef Bus b = buses[0] if buses else None
        if b:
            cnet.NET_add_vargens(self._c_net,b._c_ptr,penetration,uncertainty,corr_radius,corr_value)
        else:
            cnet.NET_add_vargens(self._c_net,NULL,penetration,uncertainty,corr_radius,corr_value)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def adjust_generators(self):
        """
        Adjusts powers of slack and regulator generators connected to or regulating the
        same bus to correct generator participations without modifying the total power injected.
        """
        
        cnet.NET_adjust_generators(self._c_net);

    def clear_error(self):
        """ 
        Clear error flag and message string.
        """

        cnet.NET_clear_error(self._c_net);

    def clear_flags(self):
        """
        Clears all the flags of all the network components.
        """

        cnet.NET_clear_flags(self._c_net)

    def clear_properties(self):
        """
        Clears all the network properties.
        """

        cnet.NET_clear_properties(self._c_net)

    def clear_sensitivities(self):
        """
        Clears all sensitivity information.
        """
        
        cnet.NET_clear_sensitivities(self._c_net)

    def create_sorted_bus_list(self,sort_by,t=0):
        """
        Creates list of buses sorted in descending order according to a specific quantity.

        Parameters
        ----------
        sort_by : int (:ref:`ref_bus_sens`, :ref:`ref_bus_mis`).
        t : int

        Returns
        -------
        buses : list of :class:`Buses <pfnet.Bus>`
        """

        buses = []
        cdef cbus.Bus* b = cnet.NET_create_sorted_bus_list(self._c_net,sort_by,t)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def create_vargen_P_sigma(self,spread,corr):
        """
        Creates covariance matrix (lower triangular part) for
        variable vargen active powers. 

        Parameters
        ----------
        spead : int
                Determines correlation neighborhood in terms of number of edges.
        corr : float 
               Desired correlation coefficient for neighboring vargens.
        
        Returns
        -------
        sigma : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """

        sigma = Matrix(cnet.NET_create_vargen_P_sigma(self._c_net,spread,corr),
                       owndata=True)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))
        else:
            return sigma
        
    def get_bus_by_number(self,number):
        """
        Gets bus with the given number.

        Parameters
        ----------
        number : int
        
        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        ptr = cnet.NET_bus_hash_number_find(self._c_net,number)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_bus_by_name(self,name):
        """
        Gets bus with the given name.

        Parameters
        ----------
        name : string
        
        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_bus_hash_name_find(self._c_net,name)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('bus not found')

    def get_vargen_by_name(self,name):
        """
        Gets vargen with the given name.

        Parameters
        ----------
        name : string
        
        Returns
        -------
        vargen : :class:`VarGenerator <pfnet.VarGenerator>`
        """

        name = name.encode('UTF-8')
        ptr = cnet.NET_vargen_hash_name_find(self._c_net,name)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('vargen not found')
            
    def get_bus(self,index):
        """
        Gets bus with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bus : :class:`Bus <pfnet.Bus>`
        """

        ptr = cnet.NET_get_bus(self._c_net,index)
        if ptr is not NULL:
            return new_Bus(ptr)
        else:
            raise NetworkError('invalid bus index')

    def get_branch(self,index):
        """
        Gets branch with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        branch : :class:`Branch <pfnet.Branch>`
        """
        
        ptr = cnet.NET_get_branch(self._c_net,index)
        if ptr is not NULL:
            return new_Branch(ptr)
        else:
            raise NetworkError('invalid branch index')

    def get_gen(self,index):
        """
        Gets generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Generator <pfnet.Generator>`
        """

        ptr = cnet.NET_get_gen(self._c_net,index)
        if ptr is not NULL:
            return new_Generator(ptr)
        else:
            raise NetworkError('invalid gen index')

    def get_shunt(self,index):
        """
        Gets shunt with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Shunt <pfnet.Shunt>`
        """

        ptr = cnet.NET_get_shunt(self._c_net,index)
        if ptr is not NULL:
            return new_Shunt(ptr)
        else:
            raise NetworkError('invalid shunt index')

    def get_load(self,index):
        """
        Gets load with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        gen : :class:`Load <pfnet.Load>`
        """

        ptr = cnet.NET_get_load(self._c_net,index)
        if ptr is not NULL:
            return new_Load(ptr)
        else:
            raise NetworkError('invalid load index')

    def get_vargen(self,index):
        """
        Gets variable generator with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        vargen : :class:`VarGenerator <pfnet.VarGenerator>`
        """

        ptr = cnet.NET_get_vargen(self._c_net,index)
        if ptr is not NULL:
            return new_VarGenerator(ptr)
        else:
            raise NetworkError('invalid vargen index')

    def get_bat(self,index):
        """
        Gets battery with the given index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bat : :class:`Battery <pfnet.Battery>`
        """

        ptr = cnet.NET_get_bat(self._c_net,index)
        if ptr is not NULL:
            return new_Battery(ptr)
        else:
            raise NetworkError('invalid battery index')

    def get_gen_buses(self):
        """
        Gets list of buses where generators are connected.

        Returns
        -------
        buses : list
        """
        
        buses = []
        cdef cbus.Bus* b = cnet.NET_get_gen_buses(self._c_net)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def get_load_buses(self):
        """
        Gets list of buses where loads are connected.

        Returns
        -------
        buses : list
        """
        
        buses = []
        cdef cbus.Bus* b = cnet.NET_get_load_buses(self._c_net)
        while b is not NULL:
            buses.append(new_Bus(b))
            b = cbus.BUS_get_next(b)
        return buses

    def get_var_values(self,code=CURRENT):
        """
        Gets network variable values.

        Parameters
        ----------
        code : int (See var values)

        Returns
        -------
        values : :class:`ndarray <numpy.ndarray>`
        """
        return Vector(cnet.NET_get_var_values(self._c_net,code),owndata=True)

    def get_var_projection(self,obj_type,var):
        """
        Gets projection matrix for specific object variables.

        Parameters
        ----------
        obj_type : int (:ref:`ref_net_obj`)
        var : int (:ref:`ref_bus_var`, :ref:`ref_branch_var`, :ref:`ref_gen_var`, :ref:`ref_shunt_var`, :ref:`ref_load_var`, :ref:`ref_vargen_var`, :ref:`ref_bat_var`) 
        """
        
        return Matrix(cnet.NET_get_var_projection(self._c_net,obj_type,var),owndata=True)

    def get_num_buses(self):
        """
        Gets number of buses in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_buses(self._c_net)

    def get_num_slack_buses(self):
        """
        Gets number of slack buses in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_slack_buses(self._c_net)

    def get_num_buses_reg_by_gen(self):
        """
        Gets number of buses whose voltage magnitudes are regulated by generators.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_buses_reg_by_gen(self._c_net)

    def get_num_buses_reg_by_tran(self,only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by tap-changing transformers.

        Returns
        -------
        num : int
        """

        if not only:
            return cnet.NET_get_num_buses_reg_by_tran(self._c_net)
        else:
            return cnet.NET_get_num_buses_reg_by_tran_only(self._c_net)

    def get_num_buses_reg_by_shunt(self,only=False):
        """
        Gets number of buses whose voltage magnitudes are regulated by switched shunt devices.

        Returns
        -------
        num : int
        """
        
        if not only:
            return cnet.NET_get_num_buses_reg_by_shunt(self._c_net)
        else:
            return cnet.NET_get_num_buses_reg_by_shunt_only(self._c_net)

    def get_num_branches(self):
        """
        Gets number of branches in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_branches(self._c_net)

    def get_num_branches_not_on_outage(self):
        """
        Gets number of branches in the network that are not on outage.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_branches_not_on_outage(self._c_net)

    def get_num_fixed_trans(self):
        """
        Gets number of fixed transformers in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_fixed_trans(self._c_net)

    def get_num_lines(self):
        """
        Gets number of transmission lines in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_lines(self._c_net)

    def get_num_phase_shifters(self):
        """
        Gets number of phase-shifting transformers in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_phase_shifters(self._c_net)

    def get_num_tap_changers(self):
        """
        Gets number of tap-changing transformers in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_tap_changers(self._c_net)

    def get_num_tap_changers_v(self):
        """
        Gets number of tap-changing transformers in the network that regulate voltage magnitudes.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_tap_changers_v(self._c_net)

    def get_num_tap_changers_Q(self):
        """
        Gets number of tap-changing transformers in the network that regulate reactive flows.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_tap_changers_Q(self._c_net)

    def get_num_generators(self):
        """
        Gets number of generators in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_gens(self._c_net)

    def get_num_gens_not_on_outage(self):
        """
        Gets number of generators in the network that are not on outage.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_gens_not_on_outage(self._c_net)

    def get_num_reg_gens(self):
        """
        Gets number generators in the network that provide voltage regulation.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_reg_gens(self._c_net)

    def get_num_slack_gens(self):
        """
        Gets number of slack generators in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_slack_gens(self._c_net)

    def get_num_P_adjust_gens(self):
        """
        Gets number of generators in the network that have adjustable active powers.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_P_adjust_gens(self._c_net)

    def get_num_loads(self):
        """
        Gets number of loads in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_loads(self._c_net)

    def get_num_P_adjust_loads(self):
        """
        Gets number of loads in the network that have adjustable active powers.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_P_adjust_loads(self._c_net)

    def get_num_shunts(self):
        """
        Gets number of shunts in the network.
        
        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_shunts(self._c_net)

    def get_num_fixed_shunts(self):
        """
        Gets number of fixed shunts in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_fixed_shunts(self._c_net)

    def get_num_switched_shunts(self):
        """
        Gets number of switched shunts in the network.

        Returns
        -------
        num : int
        """

        return cnet.NET_get_num_switched_shunts(self._c_net)

    def get_num_var_generators(self):
        """
        Gets number of variable generators in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_vargens(self._c_net)

    def get_num_batteries(self):
        """
        Gets number of batteries in the network.

        Returns
        -------
        num : int
        """
        
        return cnet.NET_get_num_bats(self._c_net)

    def get_properties(self):
        """
        Gets network properties.

        Returns
        -------
        properties : dict
        """

        return {'bus_v_max': self.bus_v_max,
                'bus_v_min': self.bus_v_min,
                'bus_v_vio': self.bus_v_vio,
                'bus_P_mis': self.bus_P_mis,
                'bus_Q_mis': self.bus_Q_mis,
                'gen_P_cost': self.gen_P_cost,
                'gen_v_dev': self.gen_v_dev,
                'gen_Q_vio': self.gen_Q_vio,
                'gen_P_vio': self.gen_P_vio,
                'tran_v_vio': self.tran_v_vio,
                'tran_r_vio': self.tran_r_vio,
                'tran_p_vio': self.tran_p_vio,
                'shunt_v_vio': self.shunt_v_vio,
                'shunt_b_vio': self.shunt_b_vio,
                'load_P_util': self.load_P_util,
                'load_P_vio': self.load_P_vio, 
                'num_actions': self.num_actions}

    def has_error(self):
        """
        Indicates whether the network has the error flag set due to an
        invalid operation.
        """

        return cnet.NET_has_error(self._c_net)

    def load(self,filename,output_level=0):
        """
        Loads a network data contained in a specific file.
        
        Parameters
        ----------
        filename : string
        output_level : int
        """

        filename = filename.encode('UTF-8')
        cnet.NET_load(self._c_net,filename,output_level)
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def set_flags(self,obj_type,flags,props,vals):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj_type : int (:ref:`ref_net_obj`)
        flags : int or list (:ref:`ref_net_flag`)
        props : int or list (:ref:`ref_bus_prop`, :ref:`ref_branch_prop`, :ref:`ref_gen_prop`, :ref:`ref_shunt_prop`, :ref:`ref_load_prop`, :ref:`ref_vargen_prop`, :ref:`ref_bat_prop`)
        vals : int or list (:ref:`ref_bus_var`, :ref:`ref_branch_var`, :ref:`ref_gen_var`, :ref:`ref_shunt_var`, :ref:`ref_load_var`, :ref:`ref_vargen_var`, :ref:`ref_bat_var`)
        """

        props = props if isinstance(props,list) else [props]
        vals = vals if isinstance(vals,list) else [vals]
        flags = flags if isinstance(flags,list) else [flags]
        cnet.NET_set_flags(self._c_net,
                           obj_type,
                           reduce(lambda x,y: x|y,flags,0),
                           reduce(lambda x,y: x|y,props,0),
                           reduce(lambda x,y: x|y,vals,0))
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))

    def set_flags_of_component(self,obj,flags,vals):
        """
        Sets flags of network components with specific properties.

        Parameters
        ----------
        obj : :class:`Bus <pfnet.Bus>`, :class:`Branch <pfnet.Branch>`, :class:`Generator <pfnet.Generator>`, :class:`Load <pfnet.Load>`, :class:`Shunt <pfnet.Shunt>`, :class:`VarGenerator <pfnet.VarGenerator>`, :class:`Battery <pfnet.Battery>` 
        flags : int or list (:ref:`ref_net_flag`)
        vals : int or list (:ref:`ref_bus_var`, :ref:`ref_branch_var`, :ref:`ref_gen_var`, :ref:`ref_shunt_var`, :ref:`ref_load_var`, :ref:`ref_vargen_var`, :ref:`ref_bat_var`)
        """
        
        cdef CPtr ptr = obj._get_c_ptr()
        vals = vals if isinstance(vals,list) else [vals]
        flags = flags if isinstance(flags,list) else [flags]
        cnet.NET_set_flags_of_component(self._c_net,
                                        ptr._c_ptr,
                                        obj.obj_type,
                                        reduce(lambda x,y: x|y,flags,0),
                                        reduce(lambda x,y: x|y,vals,0))
        if cnet.NET_has_error(self._c_net):
            raise NetworkError(cnet.NET_get_error_string(self._c_net))


    def set_var_values(self,values):
        """
        Sets network variable values.

        Parameters
        ----------
        values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if values.size else NULL
        cnet.NET_set_var_values(self._c_net,v)
    
    def show_components(self):
        """
        Shows information about the number of network components of each type.
        """
        
        print(cnet.NET_get_show_components_str(self._c_net).decode('UTF-8'))
	 	 
    def show_properties(self,t=0):
        """
        Shows information about the state of the network component quantities.
        
        Parameters
        ----------
        t : int (time period)
        """

        print(cnet.NET_get_show_properties_str(self._c_net,t).decode('UTF-8'))

    def show_buses(self,number,sort_by,t=0):
        """
        Shows information about the most relevant network buses sorted by a specific quantity.

        Parameters
        ----------
        number : int
        sort_by : int (:ref:`ref_bus_sens`, :ref:`ref_bus_mis`)
        t : int (time period)
        """
        
        cnet.NET_show_buses(self._c_net,number,sort_by,t)

    def update_properties(self,values=None):
        """
        Re-computes the network properties using the given values
        of the network variables. If no values are given, then the
        current values of the network variables are used.

        Parameters
        ----------
        values : :class:`ndarray <numpy.ndarray>`
        """
            
        cdef np.ndarray[double,mode='c'] x = values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if (values is not None and values.size) else NULL
        cnet.NET_update_properties(self._c_net,v)

    def update_set_points(self):
        """
        Updates voltage magnitude set points of gen-regulated buses 
        to be equal to the bus voltage magnitudes.
        """

        cnet.NET_update_set_points(self._c_net)

    property num_periods:
        """ Number of time periods (int). """
        def __get__(self): return cnet.NET_get_num_periods(self._c_net)

    property base_power:
        """ System base power (MVA) (float). """
        def __get__(self): return cnet.NET_get_base_power(self._c_net)

    property buses:
        """ List of network :class:`buses <pfnet.Bus>` (list). """
        def __get__(self):
            return [self.get_bus(i) for i in range(self.num_buses)]

    property branches:
        """ List of network :class:`branches <pfnet.Branch>` (list). """
        def __get__(self):
            return [self.get_branch(i) for i in range(self.num_branches)]

    property generators:
        """ List of network :class:`generators <pfnet.Generator>` (list). """
        def __get__(self):
            return [self.get_gen(i) for i in range(self.num_generators)]

    property shunts:
        """ List of network :class:`shunts <pfnet.Shunt>` (list). """
        def __get__(self):
            return [self.get_shunt(i) for i in range(self.num_shunts)]

    property loads:
        """ List of network :class:`loads <pfnet.Load>` (list). """
        def __get__(self):
            return [self.get_load(i) for i in range(self.num_loads)]

    property var_generators:
        """ List of network :class:`variable generators <pfnet.VarGenerator>` (list). """
        def __get__(self):
            return [self.get_vargen(i) for i in range(self.num_var_generators)]

    property batteries:
        """ List of network :class:`batteries <pfnet.Battery>` (list). """
        def __get__(self):
            return [self.get_bat(i) for i in range(self.num_batteries)]
    
    property num_buses:
        """ Number of buses in the network (int). """
        def __get__(self): return cnet.NET_get_num_buses(self._c_net)

    property num_branches:
        """ Number of branches in the network (int). """
        def __get__(self): return cnet.NET_get_num_branches(self._c_net)

    property num_generators:
        """ Number of generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_gens(self._c_net)

    property num_loads:
        """ Number of loads in the network (int). """
        def __get__(self): return cnet.NET_get_num_loads(self._c_net)

    property num_shunts:
        """ Number of shunt devices in the network (int). """
        def __get__(self): return cnet.NET_get_num_shunts(self._c_net)

    property num_var_generators:
        """ Number of variable generators in the network (int). """
        def __get__(self): return cnet.NET_get_num_vargens(self._c_net)

    property num_batteries:
        """ Number of batteries in the network (int). """
        def __get__(self): return cnet.NET_get_num_bats(self._c_net)

    property num_vars:
        """ Number of network quantities that have been set to variable (int). """
        def __get__(self): return cnet.NET_get_num_vars(self._c_net)

    property num_fixed:
        """ Number of network quantities that have been set to fixed (int). """
        def __get__(self): return cnet.NET_get_num_fixed(self._c_net)

    property num_bounded:
        """ Number of network quantities that have been set to bounded (int). """
        def __get__(self): return cnet.NET_get_num_bounded(self._c_net)

    property num_sparse:
        """ Number of network control quantities that have been set to sparse (int). """
        def __get__(self): return cnet.NET_get_num_sparse(self._c_net)

    property bus_v_max:
        """ Maximum bus voltage magnitude (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_bus_v_max(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus_v_min:
        """ Minimum bus voltage magnitude (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_bus_v_min(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus_v_vio:
        """ Maximum bus voltage magnitude limit violation (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_bus_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus_P_mis:
        """ Largest bus active power mismatch in the network (MW) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_bus_P_mis(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property bus_Q_mis:
        """ Largest bus reactive power mismatch in the network (MVAr) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_bus_Q_mis(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property gen_P_cost:
        """ Total active power generation cost ($/hr) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_gen_P_cost(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property gen_v_dev:
        """ Largest voltage magnitude deviation from set point of bus regulated by generator (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_gen_v_dev(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property gen_Q_vio:
        """ Largest generator reactive power limit violation (MVAr) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_gen_Q_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property gen_P_vio:
        """ Largest generator active power limit violation (MW) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_gen_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property tran_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by transformer (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_tran_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property tran_r_vio:
        """ Largest transformer tap ratio limit violation (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_tran_r_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property tran_p_vio:
        """ Largest transformer phase shift limit violation (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_tran_p_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property shunt_v_vio:
        """ Largest voltage magnitude band violation of voltage regulated by switched shunt device (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_shunt_v_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property shunt_b_vio:
        """ Largest switched shunt susceptance limit violation (p.u.) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_shunt_b_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property load_P_util:
        """ Total active power consumption utility ($/hr) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_load_P_util(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property load_P_vio:
        """ Largest load active power limit violation (MW) (float or array). """
        def __get__(self): 
            r = [cnet.NET_get_load_P_vio(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)
            
    property num_actions:
        """ Number of control adjustments (int or array). """
        def __get__(self): 
            r = [cnet.NET_get_num_actions(self._c_net,t) for t in range(self.num_periods)]
            if OPTION_SP_SCALARS and self.num_periods == 1:
                return r[0]
            else:
                return np.array(r)

    property vargen_corr_radius:
        """ Correlation radius of variable generators (number of edges). """
        def __get__(self): return cnet.NET_get_vargen_corr_radius(self._c_net)

    property vargen_corr_value:
        """ Correlation value (coefficient) of variable generators. """
        def __get__(self): return cnet.NET_get_vargen_corr_value(self._c_net)

cdef new_Network(cnet.Net* n):
    if n is not NULL:
        net = Network(alloc=False)
        net._c_net = n
        return net
    else:
        raise NetworkError('no network data')

# Contingency
#############

class ContingencyError(Exception):
    """
    Contingency error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Contingency:
    """
    Contingency class.
    """
    
    cdef ccont.Cont* _c_cont
    cdef bint alloc
    
    def __init__(self,gens=None,branches=None,alloc=True):
        """
        Contingency class.
        
        Parameters
        ----------
        gens : list or :class:`Generators <pfnet.Generator>`
        branches : list :class:`Branchs <pfnet.Branch>`
        alloc : {``True``, ``False``}
        """

        pass
     
    def __cinit__(self,gens=None,branches=None,alloc=True):
        
        cdef Generator g
        cdef Branch br

        if alloc:
            self._c_cont = ccont.CONT_new()
        else:
            self._c_cont = NULL
        self.alloc = alloc
        
        if gens:
            for gen in gens:
                g = gen
                ccont.CONT_add_gen_outage(self._c_cont,g._c_ptr)
        if branches:
            for branch in branches:
                br = branch
                ccont.CONT_add_branch_outage(self._c_cont,br._c_ptr)

    def __dealloc__(self):
        """
        Frees contingency C data structure. 
        """
        
        if self.alloc:
            ccont.CONT_del(self._c_cont)
            self._c_cont = NULL

    def apply(self):
        """
        Applies outages that characterize contingency.
        """
        
        ccont.CONT_apply(self._c_cont)

    def clear(self):
        """
        Clears outages that characterize contingency.
        """
        
        ccont.CONT_clear(self._c_cont)

    def show(self):
        """
        Shows contingency information.
        """
        
        print(ccont.CONT_get_show_str(self._c_cont).decode('UTF-8'))

    def add_gen_outage(self,gen):
        """
        Adds generator outage to contingency.
    
        Parameters
        ----------
        gen : :class:`Generator <pfnet.Generator>`
        """

        cdef Generator g = gen
        ccont.CONT_add_gen_outage(self._c_cont,g._c_ptr)

    def add_branch_outage(self,br):
        """
        Adds branch outage to contingency.
    
        Parameters
        ----------
        br : :class:`Branch <pfnet.Branch>`
        """

        cdef Branch b = br
        ccont.CONT_add_branch_outage(self._c_cont,b._c_ptr)

    def has_gen_outage(self,gen):
        """
        Determines whether contingency specifies the given generator as being on outage.

        Parameters
        ----------
        gen : :class:`Generator <pfnet.Generator>`

        Returns
        -------
        result : {``True``, ``False``}
        """
        
        cdef Generator g = gen
        return ccont.CONT_has_gen_outage(self._c_cont,g._c_ptr)

    def has_branch_outage(self,br):
        """
        Determines whether contingency specifies the given branch as being on outage.

        Parameters
        ----------
        branch : :class:`Branch <pfnet.Branch>`

        Returns
        -------
        result : {``True``, ``False``}
        """
        
        cdef Branch b = br
        return ccont.CONT_has_branch_outage(self._c_cont,b._c_ptr)

    property num_gen_outages:
        """ Number of generator outages. """
        def __get__(self): return ccont.CONT_get_num_gen_outages(self._c_cont)

    property num_branch_outages:
        """ Number of branch outages. """
        def __get__(self): return ccont.CONT_get_num_branch_outages(self._c_cont)

# Graph
#######

class GraphError(Exception):
    """
    Graph error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Graph:
    """
    Graph class.
    """

    cdef cgraph.Graph* _c_graph
    cdef cnet.Net* _c_net
    cdef bint alloc
    
    def __init__(self,net,alloc=True):
        """
        Graph class.
        
        Parameters
        ----------
        net : :class:`Network <pfnet.Network>`
        alloc : {``True``, ``False``}
        """

        pass
     
    def __cinit__(self,Network net,alloc=True):
        
        self._c_net = net._c_net
        if alloc:
            self._c_graph = cgraph.GRAPH_new(net._c_net)
        else:
            self._c_graph = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees graph C data structure. 
        """
        
        if self.alloc:
            cgraph.GRAPH_del(self._c_graph)
            self._c_graph = NULL
    
    def has_viz(self):
        """
        Determines whether graph has visualization 
        capabilities.

        Returns
        -------
        flag : {``True``, ``False``}
        """
        
        return cgraph.GRAPH_can_viz(self._c_graph)

    def has_error(self):
        """
        Indicates whether the graph has the error flag set due to an
        invalid operation.
        """

        return cgraph.GRAPH_has_error(self._c_graph)

    def clear_error(self):
        """ 
        Clear error flag and message string.
        """

        cgraph.GRAPH_clear_error(self._c_graph);
        
    def set_layout(self):
        """
        Determines and saves a layout for the graph nodes. 
        """

        cgraph.GRAPH_set_layout(self._c_graph)

    def set_node_property(self,bus,prop,value):
        """
        Sets property of node. See `Graphviz documentation <http://www.graphviz.org/Documentation.php>`_.

        Parameters
        ----------
        bus : :class:`Bus <pfnet.Bus>`
        prop : string
        value : string
        """
        
        cdef Bus b = bus
        cgraph.GRAPH_set_node_property(self._c_graph,b._c_ptr,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def set_nodes_property(self,prop,value):
        """
        Sets property of nodes. See `Graphviz documentation <http://www.graphviz.org/Documentation.php>`_.

        Parameters
        ----------
        prop : string
        value : string
        """
        
        cgraph.GRAPH_set_nodes_property(self._c_graph,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def set_edges_property(self,prop,value):
        """
        Sets property of edges. See `Graphviz documentation <http://www.graphviz.org/Documentation.php>`_.

        Parameters
        ----------
        prop : string
        value : string
        """
        
        cgraph.GRAPH_set_edges_property(self._c_graph,prop,value)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def color_nodes_by_mismatch(self,mis_type,t=0):
        """
        Colors the graphs nodes according to their power mismatch.
        
        Parameters
        ----------
        mis_type : int (:ref:`ref_bus_mis`)
        t : int
        """
        
        cgraph.GRAPH_color_nodes_by_mismatch(self._c_graph,mis_type,t)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def color_nodes_by_sensitivity(self,sens_type,t=0):
        """
        Colors the graphs nodes according to their sensitivity.

        Parameters
        ----------
        sens_type : int (:ref:`ref_bus_sens`)
        t : int
        """

        cgraph.GRAPH_color_nodes_by_sensitivity(self._c_graph,sens_type,t)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

    def view(self, inline=False):
        """
        Displays the graph.
        """
        
        temp = tempfile.NamedTemporaryFile(delete=True, suffix='.png')
        try:
            self.write("png",temp.name)

            if inline is True:
                from IPython.display import Image

                self.write("png",temp.name)
                return Image(filename=temp.name)
            else:                
                im = misc.imread(temp.name.encode('UTF-8'))
                misc.imshow(im)
        finally:
            temp.close()
            
    def write(self,format,filename):
        """
        Writes the graph to a file.

        Parameters
        ----------
        format : string (`Graphviz output formats <http://www.graphviz.org/content/output-formats>`_)
        filename : string
        """

        format = format.encode('UTF-8')
        filename = filename.encode('UTF-8')
        cgraph.GRAPH_write(self._c_graph,format,filename)
        if cgraph.GRAPH_has_error(self._c_graph):
            raise GraphError(cgraph.GRAPH_get_error_string(self._c_graph))

# Function
##########

# Types
FUNC_TYPE_REG_VMAG = cfunc.FUNC_TYPE_REG_VMAG
FUNC_TYPE_REG_VANG = cfunc.FUNC_TYPE_REG_VANG
FUNC_TYPE_REG_PQ = cfunc.FUNC_TYPE_REG_PQ
FUNC_TYPE_REG_RATIO = cfunc.FUNC_TYPE_REG_RATIO
FUNC_TYPE_REG_PHASE = cfunc.FUNC_TYPE_REG_PHASE
FUNC_TYPE_REG_SUSC = cfunc.FUNC_TYPE_REG_SUSC
FUNC_TYPE_GEN_COST = cfunc.FUNC_TYPE_GEN_COST
FUNC_TYPE_SP_CONTROLS = cfunc.FUNC_TYPE_SP_CONTROLS
FUNC_TYPE_SLIM_VMAG = cfunc.FUNC_TYPE_SLIM_VMAG
FUNC_TYPE_LOAD_UTIL = cfunc.FUNC_TYPE_LOAD_UTIL
FUNC_TYPE_NETCON_COST = cfunc.FUNC_TYPE_NETCON_COST

class FunctionError(Exception):
    """
    Function error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Function:
    """
    Function class.
    """

    cdef cfunc.Func* _c_func
    cdef bint alloc

    def __init__(self,int type,float weight,Network net,alloc=True):
        """
        Function class.
        
        Parameters
        ----------
        type : int (:ref:`ref_func_type`)
        weight : float
        net : :class:`Network <pfnet.Network>`
        alloc : {``True``, ``False``}
        """

        pass
     
    def __cinit__(self,int type,float weight,Network net,alloc=True):
        
        if alloc:
            self._c_func = cfunc.FUNC_new(type,weight,net._c_net)
        else:
            self._c_func = NULL
        self.alloc = alloc

    def __dealloc__(self):
        """
        Frees function C data structure. 
        """
        
        if self.alloc:
            cfunc.FUNC_del(self._c_func)
            self._c_func = NULL

    def del_matvec(self):
        """
        Deletes matrices and vectors associated with
        this function.
        """

        cfunc.FUNC_del_matvec(self._c_func)
            
    def update_network(self):
        """
        Updates internal arrays to be compatible
        with any network changes.
        """

        cfunc.FUNC_update_network(self._c_func)

    def clear_error(self):
        """
        Clears internal error flag.
        """

        cfunc.FUNC_clear_error(self._c_func)
        
    def analyze(self):
        """
        Analyzes function and allocates required vectors and matrices.
        """

        cfunc.FUNC_del_matvec(self._c_func)
        cfunc.FUNC_count(self._c_func)
        cfunc.FUNC_allocate(self._c_func)
        cfunc.FUNC_analyze(self._c_func)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))
        
    def eval(self,var_values):
        """
        Evaluates function value, gradient, and Hessian using 
        the given variable values. 
        
        Parameters
        ----------
        var_values : :class:`ndarray <numpy.ndarray>`
        """

        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if var_values.size else NULL
        cfunc.FUNC_eval(self._c_func,v)
        if cfunc.FUNC_has_error(self._c_func):
            raise FunctionError(cfunc.FUNC_get_error_string(self._c_func))

    property type:
        """ Function type (int). """
        def __get__(self): return cfunc.FUNC_get_type(self._c_func)

    property Hcounter:
        """ Number of nonzero entries in Hessian matrix (int). """
        def __get__(self): return cfunc.FUNC_get_Hcounter(self._c_func)
        
    property phi:
        """ Function value (float). """
        def __get__(self): return cfunc.FUNC_get_phi(self._c_func)

    property gphi:
        """ Function gradient vector (:class:`ndarray <numpy.ndarray>`). """   
        def __get__(self): return Vector(cfunc.FUNC_get_gphi(self._c_func))

    property Hphi:
        """ Function Hessian matrix (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cfunc.FUNC_get_Hphi(self._c_func))

    property weight:
        """ Function weight (float). """
        def __get__(self): return cfunc.FUNC_get_weight(self._c_func)
        
cdef new_Function(cfunc.Func* f, cnet.Net* n):
    if f is not NULL and n is not NULL:
        func = Function(0,0,new_Network(n),alloc=False)
        func._c_func = f
        return func
    else:
        raise FunctionError('invalid function data')
    
# Constraint
############

# Types
CONSTR_TYPE_PF = cconstr.CONSTR_TYPE_PF
CONSTR_TYPE_DCPF = cconstr.CONSTR_TYPE_DCPF
CONSTR_TYPE_LINPF = cconstr.CONSTR_TYPE_LINPF
CONSTR_TYPE_FIX = cconstr.CONSTR_TYPE_FIX
CONSTR_TYPE_BOUND = cconstr.CONSTR_TYPE_BOUND
CONSTR_TYPE_PAR_GEN_P = cconstr.CONSTR_TYPE_PAR_GEN_P
CONSTR_TYPE_PAR_GEN_Q = cconstr.CONSTR_TYPE_PAR_GEN_Q
CONSTR_TYPE_REG_GEN = cconstr.CONSTR_TYPE_REG_GEN
CONSTR_TYPE_REG_TRAN = cconstr.CONSTR_TYPE_REG_TRAN
CONSTR_TYPE_REG_SHUNT = cconstr.CONSTR_TYPE_REG_SHUNT
CONSTR_TYPE_DC_FLOW_LIM = cconstr.CONSTR_TYPE_DC_FLOW_LIM
CONSTR_TYPE_LBOUND = cconstr.CONSTR_TYPE_LBOUND

class ConstraintError(Exception):
    """
    Constraint error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Constraint:
    """
    Constraint class.
    """
    
    cdef cconstr.Constr* _c_constr
    cdef cnet.Net* _c_net
    cdef bint alloc

    def __init__(self,int type,Network net,alloc=True):
        """
        Contraint class.

        Parameters
        ----------
        type : int (:ref:`ref_constr_type`)
        net : :class:`Network <pfnet.Network>`
        alloc : {``True``, ``False``}
        """

        pass
     
    def __cinit__(self,int type,Network net,alloc=True):

        self._c_net = net._c_net
        if alloc:
            self._c_constr = cconstr.CONSTR_new(type,net._c_net)
        else:
            self._c_constr = NULL
        self.alloc = alloc
            
    def __dealloc__(self):
        """
        Frees constraint C data structure. 
        """
        
        if self.alloc:
            cconstr.CONSTR_del(self._c_constr)
            self._c_constr = NULL

    def del_matvec(self):
        """
        Deletes matrices and vectors associated with
        this constraint.
        """

        cconstr.CONSTR_del_matvec(self._c_constr)

    def update_network(self):
        """
        Updates internal arrays to be compatible
        with any network changes.
        """
        
        cconstr.CONSTR_update_network(self._c_constr)

    def clear_error(self):
        """
        Clears internal error flag.
        """

        cconstr.CONSTR_clear_error(self._c_constr)
        
    def analyze(self):
        """
        Analyzes constraint and allocates required vectors and matrices.
        """

        cconstr.CONSTR_del_matvec(self._c_constr)        
        cconstr.CONSTR_count(self._c_constr)
        cconstr.CONSTR_allocate(self._c_constr)     
        cconstr.CONSTR_analyze(self._c_constr)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def combine_H(self,coeff,ensure_psd=False):
        """
        Forms and saves a linear combination of the individual constraint Hessians.

        Parameters
        ----------
        coeff : :class:`ndarray <numpy.ndarray>`
        ensure_psd : {``True``, ``False``}
        """
        
        cdef np.ndarray[double,mode='c'] x = coeff
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if coeff.size else NULL
        cconstr.CONSTR_combine_H(self._c_constr,v,ensure_psd)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))
        
    def eval(self,var_values):
        """
        Evaluates constraint violations, Jacobian, and individual Hessian matrices.
        
        Parameters
        ----------
        var_values : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if var_values.size else NULL
        cconstr.CONSTR_eval(self._c_constr,v)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def store_sensitivities(self,sA,sf,sGu,sGl):
        """
        Stores Lagrange multiplier estimates of the constraints in 
        the power network components.

        Parameters
        ----------
        sA : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear equality constraints (:math:`Ax = b`)
        sf : :class:`ndarray <numpy.ndarray>`
             sensitivities for nonlinear equality constraints (:math:`f(x) = 0`)
        sGu : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`Gx \le u`)
        sGl : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`l \le Gx`)
        """

        cdef np.ndarray[double,mode='c'] xA = sA
        cdef np.ndarray[double,mode='c'] xf = sf
        cdef np.ndarray[double,mode='c'] xGu = sGu
        cdef np.ndarray[double,mode='c'] xGl = sGl
        cdef cvec.Vec* vA = cvec.VEC_new_from_array(&(xA[0]),len(xA)) if (sA is not None and sA.size) else NULL
        cdef cvec.Vec* vf = cvec.VEC_new_from_array(&(xf[0]),len(xf)) if (sf is not None and sf.size) else NULL
        cdef cvec.Vec* vGu = cvec.VEC_new_from_array(&(xGu[0]),len(xGu)) if (sGu is not None and sGu.size) else NULL
        cdef cvec.Vec* vGl = cvec.VEC_new_from_array(&(xGl[0]),len(xGl)) if (sGl is not None and sGl.size) else NULL
        cconstr.CONSTR_store_sens(self._c_constr,vA,vf,vGu,vGl)
        if cconstr.CONSTR_has_error(self._c_constr):
            raise ConstraintError(cconstr.CONSTR_get_error_string(self._c_constr))

    def get_H_single(self,i):
        """
        Gets the Hessian matrix (only lower triangular part) of an individual constraint.
        
        Parameters
        ----------
        i : int 
        
        Returns
        -------
        H : :class:`coo_matrix <scipy.sparse.coo_matrix>`
        """
        return Matrix(cconstr.CONSTR_get_H_single(self._c_constr,i))

    property type:
        """ Constraint type (:ref:`ref_constr_type`) (int). """
        def __get__(self): return cconstr.CONSTR_get_type(self._c_constr)

    property Acounter:
        """ Number of nonzero entries in the matrix of linear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_Acounter(self._c_constr)

    property Gcounter:
        """ Number of nonzero entries in the matrix of linear inequality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_Gcounter(self._c_constr)

    property Jcounter:
        """ Number of nonzero entries in the Jacobian matrix of the nonlinear equality constraints (int). """
        def __get__(self): return cconstr.CONSTR_get_Jcounter(self._c_constr)

    property Aconstr_index:
        """ Index of linear equality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_Aconstr_index(self._c_constr)

    property Gconstr_index:
        """ Index of linear ineqquality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_Gconstr_index(self._c_constr)

    property Jconstr_index:
        """ Index of nonlinear equality constraint (int). """
        def __get__(self): return cconstr.CONSTR_get_Jconstr_index(self._c_constr)
        
    property f:
        """ Vector of nonlinear equality constraint violations (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_f(self._c_constr))

    property J:
        """ Jacobian matrix of nonlinear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_J(self._c_constr))

    property b:
        """ Right-hand side vector of linear equality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_b(self._c_constr))
        
    property A:
        """ Matrix for linear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_A(self._c_constr))

    property l:
        """ Lower bound vector of linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_l(self._c_constr))

    property u:
        """ Upper bound vector of linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cconstr.CONSTR_get_u(self._c_constr))

    property G:
        """ Matrix for linear inequality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_G(self._c_constr))

    property H_combined:
        """ Linear combination of Hessian matrices of individual nonlinear equality constraints (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cconstr.CONSTR_get_H_combined(self._c_constr))

cdef new_Constraint(cconstr.Constr* c, cnet.Net* n):
    if c is not NULL and n is not NULL:
        constr = Constraint(0,new_Network(n),alloc=False)
        constr._c_constr = c
        return constr
    else:
        raise ConstraintError('invalid constraint data')

# Heuristic
###########

# Types
HEUR_TYPE_PVPQ = cheur.HEUR_TYPE_PVPQ

class HeuristicError(Exception):
    """
    Heuristic error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Heuristic:
    """
    Heuristic class.
    """

    pass

# Problem
#########        
    
class ProblemError(Exception):
    """
    Problem error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class Problem:
    """
    Optimization problem class.
    """
    
    cdef cprob.Prob* _c_prob
    cdef bint alloc

    def __init__(self):
        """
        Optimization problem class.
        """

        pass        

    def __cinit__(self):

        self._c_prob = cprob.PROB_new()
        self.alloc = True

    def __dealloc__(self):
        """
        Frees problem C data structure. 
        """
        
        if self.alloc:
            cprob.PROB_del(self._c_prob)
            self._c_prob = NULL
        
    def add_constraint(self,ctype):
        """
        Adds constraint to optimization problem.
        
        Parameters
        ----------
        ctype : int (:ref:`ref_constr_type`)
        """

        cprob.PROB_add_constr(self._c_prob,ctype)

    def add_function(self,ftype,weight):
        """
        Adds function to optimization problem objective.

        Parameters
        ----------
        ftype : int (:ref:`ref_func_type`)
        weight : float
        """

        cprob.PROB_add_func(self._c_prob,ftype,weight)

    def add_heuristic(self,htype):

        cprob.PROB_add_heur(self._c_prob,htype)

    def analyze(self):
        """
        Analyzes function and constraint structures and allocates
        required vectors and matrices.
        """

        cprob.PROB_analyze(self._c_prob)

    def apply_heuristics(self,var_values):
        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if var_values.size else NULL
        cprob.PROB_apply_heuristics(self._c_prob,v)

    def clear(self):
        """
        Resets optimization problem data.
        """

        cprob.PROB_clear(self._c_prob)

    def combine_H(self,coeff,ensure_psd):
        """
        Forms and saves a linear combination of the individual constraint Hessians.

        Parameters
        ----------
        coeff : :class:`ndarray <numpy.ndarray>`
        ensure_psd : {``True``, ``False``}
        """
        
        cdef np.ndarray[double,mode='c'] x = coeff
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if coeff.size else NULL
        cprob.PROB_combine_H(self._c_prob,v,ensure_psd)

    def eval(self,var_values):
        """
        Evaluates objective function and constraints as well as their first and
        second derivatives using the given variable values. 
        
        Parameters
        ----------
        var_values : :class:`ndarray <numpy.ndarray>`
        """
        
        cdef np.ndarray[double,mode='c'] x = var_values
        cdef cvec.Vec* v = cvec.VEC_new_from_array(&(x[0]),len(x)) if var_values.size else NULL
        cprob.PROB_eval(self._c_prob,v)

    def store_sensitivities(self,sA,sf,sGu,sGl):
        """
        Stores Lagrange multiplier estimates of the constraints in 
        the power network components.

        Parameters
        ----------
        sA : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear equality constraints (:math:`Ax = b`)
        sf : :class:`ndarray <numpy.ndarray>`
             sensitivities for nonlinear equality constraints (:math:`f(x) = 0`)
        sGu : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`Gx \le u`)
        sGl : :class:`ndarray <numpy.ndarray>`
             sensitivities for linear inequality constraints (:math:`l \le Gx`)
        """
        
        cdef np.ndarray[double,mode='c'] xA = sA
        cdef np.ndarray[double,mode='c'] xf = sf
        cdef np.ndarray[double,mode='c'] xGu = sGu
        cdef np.ndarray[double,mode='c'] xGl = sGl
        cdef cvec.Vec* vA = cvec.VEC_new_from_array(&(xA[0]),len(xA)) if (sA is not None and sA.size) else NULL
        cdef cvec.Vec* vf = cvec.VEC_new_from_array(&(xf[0]),len(xf)) if (sf is not None and sf.size) else NULL
        cdef cvec.Vec* vGu = cvec.VEC_new_from_array(&(xGu[0]),len(xGu)) if (sGu is not None and sGu.size) else NULL
        cdef cvec.Vec* vGl = cvec.VEC_new_from_array(&(xGl[0]),len(xGl)) if (sGl is not None and sGl.size) else NULL
        cprob.PROB_store_sens(self._c_prob,vA,vf,vGu,vGl)
        if cprob.PROB_has_error(self._c_prob):
            raise ProblemError(cprob.PROB_get_error_string(self._c_prob))

    def find_constraint(self,type):
        """
        Finds constraint of give type among the constraints of this optimization problem.

        Parameters
        ----------
        type : int (:ref:`ref_constr_type`)
        """
        
        cdef cnet.Net* n = cprob.PROB_get_network(self._c_prob)
        c = cprob.PROB_find_constr(self._c_prob,type)
        if c is not NULL:
            return new_Constraint(c,n)
        else:
            raise ProblemError('constraint not found')

    def get_init_point(self):
        """
        Gets initial solution estimate from the current value of the network variables.

        Returns
        -------
        point : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_init_point(self._c_prob),owndata=True)

    def get_upper_limits(self):
        """
        Gets vector of upper limits for the network variables.

        Returns
        -------
        limits : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_upper_limits(self._c_prob),owndata=True)

    def get_lower_limits(self):
        """
        Gets vector of lower limits for the network variables.
        
        Returns
        -------
        limits : :class:`ndarray <numpy.ndarray>`
        """

        return Vector(cprob.PROB_get_lower_limits(self._c_prob),owndata=True)

    def get_network(self):
        """
        Gets the power network associated with this optimization problem.
        """

        return new_Network(cprob.PROB_get_network(self._c_prob))

    def set_network(self,net):
        """
        Sets the power network associated with this optimization problem.
        """
        
        cdef Network n = net
        cprob.PROB_set_network(self._c_prob,n._c_net)

    def show(self):
        """
        Shows information about this optimization problem.
        """
        
        print(cprob.PROB_get_show_str(self._c_prob).decode('UTF-8'))

    def update_lin(self):
        """
        Updates linear equality constraints.
        """

        cprob.PROB_update_lin(self._c_prob)

    property network:
        """ Power network associated with this optimization problem (:class:`Network <pfnet.Network>`). """
        def __get__(self): return new_Network(cprob.PROB_get_network(self._c_prob))
        def __set__(self,net):
            cdef Network n = net
            cprob.PROB_set_network(self._c_prob,n._c_net)
        
    property constraints:
        """ List of :class:`constraints <pfnet.Constraint>` of this optimization problem (list). """
        def __get__(self):
            clist = []
            cdef cconstr.Constr* c = cprob.PROB_get_constr(self._c_prob)
            cdef cnet.Net* n = cprob.PROB_get_network(self._c_prob)
            while c is not NULL:
                clist.append(new_Constraint(c,n))
                c = cconstr.CONSTR_get_next(c)
            return clist

    property functions:
        """ List of :class:`functions <pfnet.Function>` that form the objective function of this optimization problem (list). """
        def __get__(self):
            flist = []
            cdef cfunc.Func* f = cprob.PROB_get_func(self._c_prob)
            cdef cnet.Net* n = cprob.PROB_get_network(self._c_prob)
            while f is not NULL:
                flist.append(new_Function(f,n))
                f = cfunc.FUNC_get_next(f)
            return flist
            
    property A:
        """ Constraint matrix of linear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_A(self._c_prob))

    property b:
        """ Right hand side vectors of the linear equality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_b(self._c_prob))

    property G:
        """ Constraint matrix of linear inequality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_G(self._c_prob))

    property l:
        """ Lower bound for linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_l(self._c_prob))

    property u:
        """ Upper bound for linear inequality constraints (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_u(self._c_prob))

    property J:
        """ Jacobian matrix of the nonlinear equality constraints (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_J(self._c_prob))

    property f:
        """ Vector of nonlinear equality constraints violations (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_f(self._c_prob))

    property phi:
        """ Objective function value (float). """
        def __get__(self): return cprob.PROB_get_phi(self._c_prob)
            
    property gphi:
        """ Objective function gradient vector (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return Vector(cprob.PROB_get_gphi(self._c_prob))

    property Hphi:
        """ Objective function Hessian matrix (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_Hphi(self._c_prob))

    property H_combined:
        """ Linear combination of Hessian matrices of individual nonlinear equality constraints (only the lower triangular part) (:class:`coo_matrix <scipy.sparse.coo_matrix>`). """
        def __get__(self): return Matrix(cprob.PROB_get_H_combined(self._c_prob))

    property x:
        """ Initial primal point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return self.get_init_point()

    property lam:
        """ Initial dual point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return None

    property nu:
        """ Initial dual point (:class:`ndarray <numpy.ndarray>`). """
        def __get__(self): return None
