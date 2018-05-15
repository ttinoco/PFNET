/** @file json_macros.h
 *  @brief This file defines json macros.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#ifndef __JSON_MACROS_HEADER__
#define __JSON_MACROS_HEADER__

#include <string.h>

// constants
#define JSON_STR_BUFFER_EXTRA 10

#define JSON_start(output) { \
  strcpy(output,"{");        \
  output += 1;		     \
}

#define JSON_end(output) { \
  strcpy(output,"}");      \
  output += 1;		   \
}

#define JSON_obj(temp,output,field_name,field,index_func,end) { \
  if (field)							\
    sprintf(temp,"\"%s\":%d",field_name,index_func(field));	\
  else								\
    sprintf(temp,"\"%s\":%s",field_name,"null");		\
  strcpy(output,temp);						\
  output += strlen(temp);					\
  if (!end) {							\
    strcpy(output,",");						\
    output += 1;						\
  }								\
}

#define JSON_int(temp,output,field_name,field,end) { \
  sprintf(temp,"\"%s\":%d",field_name,field);	     \
  strcpy(output,temp);				     \
  output += strlen(temp);			     \
  if (!end) {					     \
    strcpy(output,",");				     \
    output += 1;				     \
  }						     \
}

#define JSON_str(temp,output,field_name,field,end) { \
  sprintf(temp,"\"%s\":\"%s\"",field_name,field);   \
  strcpy(output,temp);				     \
  output += strlen(temp);			     \
  if (!end) {					     \
    strcpy(output,",");				     \
    output += 1;				     \
  }						     \
}

#define JSON_float(temp,output,field_name,field,end) {\
  sprintf(temp,"\"%s\":%.15e",field_name,field);     \
  strcpy(output,temp);				      \
  output += strlen(temp);			      \
  if (!end) {					      \
    strcpy(output,",");				      \
    output += 1;				      \
  }						      \
}

#define JSON_bool(temp,output,field_name,field,end) {            \
  sprintf(temp,"\"%s\":%s",field_name,field ? "true" : "false"); \
  strcpy(output,temp);						 \
  output += strlen(temp);					 \
  if (!end) {							 \
    strcpy(output,",");						 \
    output += 1;						 \
  }								 \
}

#define JSON_array_float(temp,output,field_name,field,num,end) { \
  int i;							 \
  sprintf(temp,"\"%s\":[",field_name);				 \
  strcpy(output,temp);						 \
  output += strlen(temp);					 \
  for (i = 0; i < num; i++) {					 \
    sprintf(temp,"%.15e",field[i]);     			 \
    strcpy(output,temp);                                         \
    output += strlen(temp);                                      \
    if (i < num-1) {	                                         \
      strcpy(output,",");					 \
      output += 1;						 \
    }								 \
  }								 \
  if (!end) {							 \
    strcpy(output,"],");					 \
    output += 2;						 \
  }								 \
  else {							 \
    strcpy(output,"]");						 \
    output += 1;						 \
  }								 \
}

#define JSON_array_int(temp,output,field_name,field,num,end) { \
  int i;						       \
  sprintf(temp,"\"%s\":[",field_name);			       \
  strcpy(output,temp);					       \
  output += strlen(temp);				       \
  for (i = 0; i < num; i++) {				       \
    sprintf(temp,"%d",field[i]);			       \
    strcpy(output,temp);                                       \
    output += strlen(temp);                                    \
    if (i < num-1) {	                                       \
      strcpy(output,",");				       \
      output += 1;					       \
    }							       \
  }							       \
  if (!end) {						       \
    strcpy(output,"],");				       \
    output += 2;					       \
  }							       \
  else {						       \
    strcpy(output,"]");					       \
    output += 1;					       \
  }							       \
}

#define JSON_list_int(temp,output,field_name,obj,iter_type,list_func,field_func,next_func,end) { \
  iter_type* t;							                                 \
  sprintf(temp,"\"%s\":[",field_name);					                         \
  strcpy(output,temp);							                         \
  output += strlen(temp);						                         \
  for (t = list_func(obj); t != NULL; t = next_func(t)) {		                         \
    sprintf(temp,"%d",field_func(t));                                                            \
    strcpy(output,temp);                                                                         \
    output += strlen(temp);                                                                      \
    if (next_func(t) != NULL) {						                         \
      strcpy(output,",");						                         \
      output += 1;							                         \
    }									                         \
  }									                         \
  if (!end) {								                         \
    strcpy(output,"],");						                         \
    output += 2;							                         \
  }									                         \
  else {								                         \
    strcpy(output,"]");							                         \
    output += 1;							                         \
  }									                         \
}

#define JSON_array_json(temp,output,field_name,array,array_get,array_size,json_func,end) { \
  int i;								                   \
  sprintf(temp,"\"%s\":[",field_name);					                   \
  strcpy(output,temp);							                   \
  output += strlen(temp);						                   \
  for (i = 0; i < array_size; i++) {					                   \
    output = json_func(array_get(array,i),output);			                   \
    if (i < array_size-1) {						                   \
      strcpy(output,",");						                   \
      output += 1;							                   \
    }									                   \
  }									                   \
  if (!end) {								                   \
    strcpy(output,"],");						                   \
    output += 2;							                   \
  }									                   \
  else {								                   \
    strcpy(output,"]");							                   \
    output += 1;							                   \
  }									                   \
}
  
#endif
