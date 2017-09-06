#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport ccont

import json

class ContingencyError(Exception):
    """
    Contingency error exception.
    """

    pass

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
                ccont.CONT_add_gen_outage(self._c_cont,g.index)
        if branches:
            for branch in branches:
                br = branch
                ccont.CONT_add_branch_outage(self._c_cont,br.index)

    def __dealloc__(self):
        """
        Frees contingency C data structure.
        """

        if self.alloc:
            ccont.CONT_del(self._c_cont)
            self._c_cont = NULL

    def __getstate__(self):

        return self.json_string

    def __setstate__(self, state):

        class temp: pass
        state = json.loads(state)
        for index in state['generator_outages']:
            gen = temp()
            gen.index =index
            self.add_gen_outage(gen)
        for index in state['branch_outages']:
            br = temp()
            br.index =index
            self.add_branch_outage(br)

    def apply(self, network):
        """
        Applies outages that characterize contingency.

        Paramaters
        ----------
        network : :class:`Network <pfnet.Network>`
        """

        cdef Network n = network
        ccont.CONT_apply(self._c_cont, n._c_net)

    def clear(self, network):
        """
        Clears outages that characterize contingency.

        Paramaters
        ----------
        network : :class:`Network <pfnet.Network>`
        """

        cdef Network n = network
        ccont.CONT_clear(self._c_cont, n._c_net)

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

        ccont.CONT_add_gen_outage(self._c_cont,gen.index)

    def add_branch_outage(self,br):
        """
        Adds branch outage to contingency.

        Parameters
        ----------
        br : :class:`Branch <pfnet.Branch>`
        """

        ccont.CONT_add_branch_outage(self._c_cont,br.index)

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
        return ccont.CONT_has_gen_outage(self._c_cont,g.index)

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
        return ccont.CONT_has_branch_outage(self._c_cont,b.index)

    property num_gen_outages:
        """ Number of generator outages. """
        def __get__(self): return ccont.CONT_get_num_gen_outages(self._c_cont)

    property num_branch_outages:
        """ Number of branch outages. """
        def __get__(self): return ccont.CONT_get_num_branch_outages(self._c_cont)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = ccont.CONT_get_json_string(self._c_cont)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s
