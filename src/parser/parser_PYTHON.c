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
  module = PyImport_ImportModule("pfnet.parser");
  if (!module) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to import Python module");
    return parser;
  }

  // Class
  class = PyObject_GetAttrString(module,classname);
  if (!class) {
    parser->error_flag = TRUE;
    strcpy(parser->error_string,"unable to get Python module class");
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

  PyObject_CallMethod(parser->instance,"read","s",filename);
}

void PYTHON_PARSER_show(PYTHON_Parser* parser) {

  if (!parser)
    return;

  PyObject_CallMethod(parser->instance,"show","");
}

void PYTHON_PARSER_load(PYTHON_Parser* parser, Net* net) {
  
  if (!parser)
    return;
  
  PyObject* network = new_Network(net);

  PyObject_CallMethod(parser->instance,"load","O",network);
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
    strcpy(parser->error_string,PyString_AsString(result));
    Py_DECREF(result);
    return parser->error_string;
  }
}
