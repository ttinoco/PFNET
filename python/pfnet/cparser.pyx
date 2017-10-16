#cython: embedsignature=True

#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

cimport cparser

class ParserError(Exception):
    """
    Parser error exception.
    """

    pass

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

    def init(self):
        """
        Initializes parser data.
        """
        
        pass

    def parse(self, filename, num_periods=None):
        """
        Parsers data file.
        
        Parameters
        ----------
        filename : string
        num_periods : int

        Returns
        -------
        net : |Network|
        """

        if num_periods is None:
            num_periods = 0 # format-specific parser will use its default
        filename = filename.encode('UTF-8')
        cdef cparser.Net* net = cparser.PARSER_parse(self._c_parser, filename, num_periods)
        if cparser.PARSER_has_error(self._c_parser):
            raise ParserError(cparser.PARSER_get_error_string(self._c_parser))
        cdef Network pnet = new_Network(net)
        pnet.alloc = True
        return pnet

    def set(self, key, value):
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

    def write(self, Network net, filename):
        """
        Writes data to file.

        Parameters
        ----------
        net : |Network|
        filename : string
        """
        
        filename = filename.encode('UTF-8')
        cparser.PARSER_write(self._c_parser,net._c_net,filename)
        if cparser.PARSER_has_error(self._c_parser):
            raise ParserError(cparser.PARSER_get_error_string(self._c_parser))
        
cdef class Parser(ParserBase):

    def __init__(self,ext):
        """
        Parser class.
        
        Parameters
        ----------
        ext : string (extension of sample filename)
        """

        pass

    def __cinit__(self,ext):

        ext = ext.split('.')[-1]
        if ext == 'mat':
            self._c_parser = cparser.MAT_PARSER_new()
        elif ext == 'art':
            self._c_parser = cparser.ART_PARSER_new()
        elif ext == 'raw':
            self._c_parser = cparser.RAW_PARSER_new()
        elif ext == 'json':
            self._c_parser = cparser.JSON_PARSER_new()
        else:
            raise ParserError('invalid extension')

        self._alloc = True

cdef class ParserRAW(ParserBase):

    def __init__(self):
        """
        RAW parser class.
        """
    
        pass
        
    def __cinit__(self):
        
        self._c_parser = cparser.RAW_PARSER_new()
        self._alloc = True

cdef class ParserMAT(ParserBase):

    def __init__(self):
        """
        MAT parser class.
        """
    
        pass

    def __cinit__(self):
        
        self._c_parser = cparser.MAT_PARSER_new()
        self._alloc = True
        
cdef class ParserART(ParserBase):

    def __init__(self):
        """
        ART parser class.
        """
    
        pass

    def __cinit__(self):
        
        self._c_parser = cparser.ART_PARSER_new()
        self._alloc = True

cdef class ParserJSON(ParserBase):

    def __init__(self):
        """
        JSON parser class.
        """
    
        pass
        
    def __cinit__(self):
        
        self._c_parser = cparser.JSON_PARSER_new()
        self._alloc = True
