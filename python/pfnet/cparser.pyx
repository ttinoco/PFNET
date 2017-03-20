#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cparser

class ParserError(Exception):
    """
    Parser error exception.
    """

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

cdef class ParserBase:
    """
    Base parser class.
    """

    cdef cparser.Parser* _c_parser
    cdef bint _alloc

    def __init__(self):
        """
        Base parser class.
        """

        pass

    def __cinit__(self):

        self._c_parser = NULL
        self._alloc = False

    def __dealloc__(self):
        """
        Frees parser C data structure.
        """

        if self._alloc:
            cparser.PARSER_del(self._c_parser)
            self._c_parser = NULL
            self._alloc = False
    
    def parse(self,filename,num_periods=1):
        """
        Parsers data file.

        Parameters
        ----------
        filename : string
        num_periods : int

        Returns
        -------
        net : :class:`Network <pfnet.Network>`
        """

        filename = filename.encode('UTF-8')
        cdef cparser.Net* net = cparser.PARSER_parse(self._c_parser,filename,num_periods)
        if cparser.PARSER_has_error(self._c_parser):
            raise ParserError(cparser.PARSER_get_error_string(self._c_parser))
        return new_Network(net)

    def set(self,key,value):
        """
        Sets parser parameter.

        Parameters
        ----------
        key : string
        value : float
        """

        key = key.encode('UTF-8')
        cparser.PARSER_set(self._c_parser,key,value)
        if cparser.PARSER_has_error(self._c_parser):
            raise ParserError(cparser.PARSER_get_error_string(self._c_parser))

    def show(self):
        """
        Shows parser information/data.
        """

        cparser.PARSER_show(self._c_parser)

    def write(self,Network net, filename):
        """
        Writes data to file.

        Parameters
        ----------
        net : :class:`Network <pfnet.Network>`
        filename : string
        """
        
        filename = filename.encode('UTF-8')
        cparser.PARSER_write(self._c_parser,net._c_net,filename)
        if cparser.PARSER_has_error(self._c_parser):
            raise ParserError(cparser.PARSER_get_error_string(self._c_parser))
        
        
