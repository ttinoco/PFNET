/** @file heur.c
 *  @brief This file defines the Heur data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/heur.h>

struct Heur {

  // Error
  BOOL error_flag;                     /**< @brief Error flag */
  char error_string[HEUR_BUFFER_SIZE]; /**< @brief Error string */
  
  // Name
  char name[HEUR_BUFFER_SIZE]; /**< @brief Name string */

  // Network
  Net* net; /**< @brief Power network */

  // Type functions
  void (*func_init)(Heur* h);
  void (*func_clear)(Heur* h);
  void (*func_apply_step)(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values);
  void (*func_free)(Heur* h);

  // Type data
  void* data;

  // List
  struct Heur* next;
};

void HEUR_del(Heur* h) {
  if (h) {

    // Data
    if (h->func_free)
      (*(h->func_free))(h);
    
    // Free
    free(h);
  }
}

Net* HEUR_get_network(Heur* h) {
  if (h)
    return h->net;
  else
    return NULL;
}

char* HEUR_get_name(Heur* h) {
  if (h)
    return h->name;
  else
    return NULL;
}

void* HEUR_get_data(Heur* h) {
  if (h)
    return h->data;
  else
    return NULL;
}

Heur* HEUR_get_next(Heur* h) {
  if (h)
    return h->next;
  else
    return NULL;
}

void HEUR_list_clear_error(Heur* hlist) {
  Heur* hh;
  for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh))
    HEUR_clear_error(hh);
}

BOOL HEUR_list_has_error(Heur* hlist) {
  Heur* hh;
  for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh)) {
    if (HEUR_has_error(hh))
      return TRUE;
  }
  return FALSE;
}

char* HEUR_list_get_error_string(Heur* hlist) {
  Heur* hh;
  for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh)) {
    if (HEUR_has_error(hh))
      return HEUR_get_error_string(hh);
  }
  return "";
}

Heur* HEUR_list_add(Heur* hlist, Heur* nh) {
  LIST_add(Heur,hlist,nh,next);
  return hlist;
}

void HEUR_list_apply_step(Heur* hlist, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values) {
  Heur* hh;
  if (hlist) {
    for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh))
      HEUR_apply_step(hh,cptrs,cnum,bus,busdc,t,var_values);
  }
}

void HEUR_list_clear(Heur* hlist) {
  Heur* hh;
  if (hlist) {
    for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh))
      HEUR_clear(hh);
  }
}

void HEUR_list_del(Heur* hlist) {
  LIST_map(Heur,hlist,h,next,{HEUR_del(h);});
}

int HEUR_list_len(Heur* hlist) {
  int len;
  LIST_len(Heur,hlist,next,len);
  return len;
}

Heur* HEUR_new(Net* net) {
  Heur* h = (Heur*)malloc(sizeof(Heur));

  // Fields
  h->error_flag = FALSE;
  strcpy(h->error_string,"");
  strcpy(h->name,"unknown");
  h->net = net;
  h->data = NULL;
  h->next = NULL;
  
  // Functions
  h->func_init = NULL;
  h->func_clear = NULL;
  h->func_apply_step = NULL;
  h->func_free = NULL;

  // Update network
  HEUR_update_network(h);
  
  return h;
}

void HEUR_clear(Heur* h) {
  if (h && h->func_clear)
    (*(h->func_clear))(h);
}

void HEUR_apply(Heur* h, Constr** cptrs, int cnum, Vec* var_values) {

  // Local variables
  int i;
  int t;
  
  // No h
  if (!h)
    return;

  // Clear
  HEUR_clear(h);

  // Apply
  for (t = 0; t < NET_get_num_periods(h->net); t++) {
    for (i = 0; i < NET_get_num_buses(h->net,FALSE); i++)
      HEUR_apply_step(h,cptrs,cnum,NET_get_bus(h->net,i),NULL,t,var_values);
    for (i = 0; i < NET_get_num_dc_buses(h->net,FALSE); i++)
      HEUR_apply_step(h,cptrs,cnum,NULL,NET_get_dc_bus(h->net,i),t,var_values);
  }
}

void HEUR_apply_step(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values) {
  if (h && h->func_apply_step)
    (*(h->func_apply_step))(h,cptrs,cnum,bus,busdc,t,var_values);
}

void HEUR_clear_error(Heur * h) {
  if (h) {
    h->error_flag = FALSE;
    strcpy(h->error_string,"");
  }
}

void HEUR_init(Heur* h) {
  if (h && h->func_free)
    (*(h->func_free))(h);
  if (h && h->func_init)
    (*(h->func_init))(h);
}

BOOL HEUR_has_error(Heur* h) {
  if (h)
    return h->error_flag;
  else
    return FALSE;
}

char* HEUR_get_error_string(Heur* h) {
  if (h)
    return h->error_string;
  else
    return NULL;
}

void HEUR_set_name(Heur* h, char* name) {
  if (h)
    strcpy(h->name,name);
}

void HEUR_set_data(Heur* h, void* data) {
  if (h)
    h->data = data;
}

void HEUR_set_error(Heur* h, char* error_string) {
  if (h) {
    h->error_flag = TRUE;
    strcpy(h->error_string,error_string);
  }
}

void HEUR_set_func_init(Heur* h, void (*func)(Heur* h)) {
  if (h)
    h->func_init = func;
}

void HEUR_set_func_clear(Heur* h, void (*func)(Heur* h)) {
  if (h)
    h->func_clear = func;
}

void HEUR_set_func_apply_step(Heur* h, void (*func)(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values)) {
  if (h)
    h->func_apply_step = func;
}

void HEUR_set_func_free(Heur* h, void (*func)(Heur* h)) {
  if (h)
    h->func_free = func;
}

void HEUR_update_network(Heur* h) {
  
  // No h
  if (!h)
    return;

  // Init
  HEUR_init(h);
}
