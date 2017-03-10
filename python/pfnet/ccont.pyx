#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport ccont

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
