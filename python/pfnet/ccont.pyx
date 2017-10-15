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

    def __init__(self, generators=None, branches=None, alloc=True):
        """
        Contingency class.

        Parameters
        ----------
        generators : list of |Generator| objects
        branches : list of |Branch| objects
        alloc : |TrueFalse|
        """

        pass

    def __cinit__(self, generators=None, branches=None, alloc=True):

        cdef Generator g
        cdef Branch br

        if alloc:
            self._c_cont = ccont.CONT_new()
        else:
            self._c_cont = NULL
        self.alloc = alloc

        if generators:
            for gen in generators:
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
            gen.index = index
            self.add_generator_outage(gen)
        for index in state['branch_outages']:
            br = temp()
            br.index = index
            self.add_branch_outage(br)

    def apply(self, network):
        """
        Applies outages that characterize contingency to given network.

        Parameters
        ----------
        network : |Network|
        """

        cdef Network n = network
        ccont.CONT_apply(self._c_cont, n._c_net)

    def clear(self, network):
        """
        Clears outages that characterize contingency from given network.

        Parameters
        ----------
        network : |Network|
        """

        cdef Network n = network
        ccont.CONT_clear(self._c_cont, n._c_net)

    def show(self):
        """
        Shows contingency information.
        """

        print(ccont.CONT_get_show_str(self._c_cont).decode('UTF-8'))

    def add_generator_outage(self, gen):
        """
        Adds generator outage to contingency.

        Parameters
        ----------
        gen : |Generator|
        """

        ccont.CONT_add_gen_outage(self._c_cont, gen.index)

    def add_branch_outage(self, br):
        """
        Adds branch outage to contingency.

        Parameters
        ----------
        br : |Branch|
        """

        ccont.CONT_add_branch_outage(self._c_cont,br.index)

    def has_generator_outage(self, gen):
        """
        Determines whether contingency specifies the given generator as being on outage.

        Parameters
        ----------
        gen : |Generator|

        Returns
        -------
        result : |TrueFalse|
        """

        cdef Generator g = gen
        return ccont.CONT_has_gen_outage(self._c_cont,g.index)

    def has_branch_outage(self, br):
        """
        Determines whether contingency specifies the given branch as being on outage.

        Parameters
        ----------
        branch : |Branch|

        Returns
        -------
        result : |TrueFalse|
        """

        cdef Branch b = br
        return ccont.CONT_has_branch_outage(self._c_cont,b.index)

    property branch_outages:
       """ Array of outage branch indices (|Array|). """
       def __get__(self): return IntArray(ccont.CONT_get_branch_outages(self._c_cont),
                                          self.num_branch_outages,
                                          owndata=True)

    property generator_outages:
       """ Array of outage generator indices (|Array|). """
       def __get__(self): return IntArray(ccont.CONT_get_gen_outages(self._c_cont),
                                          self.num_generator_outages,
                                          owndata=True)

    property num_generator_outages:
        """ Number of generator outages (int). """
        def __get__(self): return ccont.CONT_get_num_gen_outages(self._c_cont)

    property num_branch_outages:
        """ Number of branch outages (int). """
        def __get__(self): return ccont.CONT_get_num_branch_outages(self._c_cont)

    property json_string:
        """ JSON string (string). """
        def __get__(self): 
            cdef char* json_string = ccont.CONT_get_json_string(self._c_cont)
            s = json_string.decode('UTF-8')
            free(json_string)
            return s
