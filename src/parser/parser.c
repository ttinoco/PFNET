/** @file parser.c
 *  @brief This file defines the Parser structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/parser.h>
#include <pfnet/parser_MAT.h>
#include <pfnet/parser_ART.h>
#include <pfnet/parser_RAW.h>

struct Parser {

  // Error
  BOOL error_flag;                        /**< @brief Parser error flag */
  char error_string[PARSER_BUFFER_SIZE];  /**< @brief Parser error string */
  
  // Data
  void * data; /**< @brief Parser data */

  // Functions
  Net* (*func_parse)(Parser* p, char* f, int n);    /**< @brief Parsing function */
  void (*func_set)(Parser* p, char* key, REAL v);   /**< @brief Configuring function */
  void (*func_show)(Parser* p);                     /**< @brief Showing function */
  void (*func_write)(Parser* p, Net* net, char* f); /**< @brief Writing function */
  void (*func_free)(Parser* p);                     /**< @brief Cleaning function */
};

Parser* PARSER_new(void) {
  
  Parser* p = (Parser*)malloc(sizeof(Parser));

  p->error_flag = FALSE;
  strcpy(p->error_string,"");

  p->data = NULL;

  p->func_parse = NULL;
  p->func_set = NULL;
  p->func_show = NULL;
  p->func_write = NULL;
  p->func_free = NULL;
}

Parser* PARSER_new_for_file(char* f) {
  char* ext;
  ext = strrchr(f,'.');
  ext = strtolower(ext);
  if (!ext)
    return NULL;
  if (strcmp(ext+1,"raw") == 0)
    return RAW_PARSER_new();
  if (strcmp(ext+1,"mat") == 0)
    return MAT_PARSER_new();
  if (strcmp(ext+1,"art") == 0)
    return ART_PARSER_new();
  return NULL;
}

Net* PARSER_parse(Parser* p, char* f, int n) {
  if (p && p->func_parse)
    return (*(p->func_parse))(p,f,n);
  else 
    return NULL;
}

void PARSER_set(Parser* p, char* key, REAL value) {
  if (p && p->func_set)
    (*(p->func_set))(p,key,value);    
}

void PARSER_show(Parser* p) {
  if (p && p->func_show)
    (*(p->func_show))(p);    
}

void PARSER_write(Parser* p, Net* net, char* f) {
  if (p && p->func_write)
    (*(p->func_write))(p,net,f);    
}

void PARSER_del(Parser* p) {
  if (p) {
    if (p->func_free)
      (*(p->func_free))(p);
    if(p->data)
      free(p->data);
    free(p);
  }
}

BOOL PARSER_has_error(Parser* p) {
  if (p)
    return p->error_flag;
  else
    return FALSE;
}

void PARSER_clear_error(Parser* p) {
  if (p) {
    p->error_flag = FALSE;
    strcpy(p->error_string,"");
  }
}

char* PARSER_get_error_string(Parser* p) {
  if (p)
    return p->error_string;
  else
    return NULL;
}

void* PARSER_get_data(Parser* p) {
  if (p)
    return p->data;
  else
    return NULL;
}

void PARSER_set_data(Parser* p, void* data) {
  if (p)
    p->data = data;
}

void PARSER_set_error(Parser* p, char* string) {
  if (p) {
    p->error_flag = TRUE;
    strcpy(p->error_string,string);
  }
}

void PARSER_set_func_parse(Parser* p, Net* (*func)(Parser* p, char* f, int n)) {
  if (p)
    p->func_parse = func;
}


void PARSER_set_func_set(Parser* p, void (*func)(Parser* p, char* key, REAL v)) {
  if (p)
    p->func_set = func;
}

void PARSER_set_func_show(Parser* p, void (*func)(Parser* p)) {
  if (p)
    p->func_show = func;
}

void PARSER_set_func_write(Parser* p, void (*func)(Parser* p, Net* net, char* f)) {
  if (p)
    p->func_write = func;
}

void PARSER_set_func_free(Parser* p, void (*func)(Parser* p)) {
  if (p)
    p->func_free = func;
}


