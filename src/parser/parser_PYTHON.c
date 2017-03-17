/** @file parser_PYTHON.c
 *  @brief This file defines the PYTHON parser interface methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser_PYTHON.h>

#if HAVE_PYTHON_PARSER

struct PYTHON_Parser {

  // Error
  BOOL error_flag;
  char error_string[PYTHON_PARSER_BUFFER_SIZE];
  
  // Instance
  PyObject* instance;
};

PYTHON_Parser* PYTHON_PARSER_new(char* classname) {

  // Local vars
  PyObject* module = NULL;
  PyObject* class = NULL;
  
  // Allocate
  PYTHON_Parser* parser = (PYTHON_Parser*)malloc(sizeof(PYTHON_Parser));
  
  // Error
  parser->error_flag = FALSE;
  strcpy(parser->error_string,"");
  
  // Init field
  parser->instance = NULL;

  // Init python
  Py_Initialize();

  // Module
  module = PyImport_ImportModule("pfnet.parsers");
  if (!module) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to import Python parser module");
    return parser;
  }

  // Class
  class = PyObject_GetAttrString(module,classname);
  if (!class) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to get Python parser class");
    return parser;
  }

  // Instance
  parser->instance = PyInstance_New(class,NULL,NULL);
  if (!parser->instance) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to instantiate Python parser");
  }

  // Return
  return parser;
}

void PYTHON_PARSER_read(PYTHON_Parser* parser, char* filename) {
  
  if (!parser)
    return;

  PyObject* result = PyObject_CallMethod(parser->instance,"read","s",filename);
  if (!result) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to call method 'read' of Python parser");
  }
}

void PYTHON_PARSER_show(PYTHON_Parser* parser) {

  if (!parser)
    return;

  PyObject* result = PyObject_CallMethod(parser->instance,"show","");
  if (!result) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to call method 'show' of Python parser");
  }
}

void PYTHON_PARSER_load(PYTHON_Parser* parser, Net* net) {
  
  if (!parser)
    return;
  
  PyObject* network = new_Network(net);

  PyObject* result = PyObject_CallMethod(parser->instance,"load","O",network);
  if (!result) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to call method 'load' of Python parser");
  }
}

void PYTHON_PARSER_del(PYTHON_Parser* parser) {
  
   if (!parser)
     return;

   // Instance
   if (parser->instance)
     Py_DECREF(parser->instance); // I hope this does it! Check!

   // Parser
   free(parser);
}

BOOL PYTHON_PARSER_has_error(PYTHON_Parser* parser) {
  if (!parser)
    return TRUE;
  else if (!parser->instance)
    return parser->error_flag;
  else {
    PyObject* result = PyObject_CallMethod(parser->instance, "has_error","");
    if (PyObject_IsTrue(result)) {
      Py_DECREF(result);
      return TRUE;
    }
    else {
      Py_DECREF(result);
      return FALSE;
    }
  }
}

char* PYTHON_PARSER_get_error_string(PYTHON_Parser* parser) {
  if (!parser)
    return "empty parser";
  else if (!parser->instance)
    return parser->error_string;
  else {
    PyObject* result = PyObject_CallMethod(parser->instance, "get_error_string","");
    if (result) {
      strcpy(parser->error_string,PyString_AsString(result));
      Py_DECREF(result);
      return parser->error_string;
    }
    else {
      strcpy(parser->error_string,"unable to call method 'get_error_string' of Python parser");
      return parser->error_string;
    }
  }
}

void PYTHON_PARSER_set(PYTHON_Parser* parser, char* key, REAL value) {
  if (!parser)
    return;
  else if (!parser->instance)
    return;
  else {
    PyObject* result = PyObject_CallMethod(parser->instance,"set","sd",key,value);
    if (!result) {
      parser->error_flag = TRUE;
      strcpy(parser->error_string,"unable to call method 'set' of Python parser");
    }
  }
}

#endif
