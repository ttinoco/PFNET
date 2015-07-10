/** @file heur.c
 *  @brief This file defines the Heur data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/heur.h>
#include <pfnet/heur_PVPQ.h>

struct Heur {

  // Type
  int type;

  // Utils
  char* bus_counted;

  // Type functions
  void (*func_init)(Heur* h, Net* net);
  void (*func_clear)(Heur* h, Net* net);
  void (*func_apply_to_branch)(Heur* h, Constr* clist, Net* net, Branch* br, Vec* var_values);
  void (*func_free)(Heur* h);

  // Type data
  void* data;

  // List
  struct Heur* next;

};

void HEUR_clear_bus_counted(Heur* h, int num) {
  int i;
  if (h) {
    if (h->bus_counted) {
      for (i = 0; i < num; i++)
	h->bus_counted[i] = 0;
    }
  } 
}

void HEUR_del(Heur* h) {
  if (h) {

    // Utils
    free(h->bus_counted);

    // Data
    if (h->func_free)
      (*(h->func_free))(h);
    
    // Free
    free(h);
  }
}

int HEUR_get_type(Heur* h) {
  if (h)
    return h->type;
  else
    return -1;
}

char* HEUR_get_bus_counted(Heur* h) {
  if (h)
    return h->bus_counted;
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

Heur* HEUR_list_add(Heur* hlist, Heur* nh) {
  LIST_add(hlist,nh,next);
  return hlist;
}

void HEUR_list_apply_to_branch(Heur* hlist, Constr* clist, Net* net, Branch* br, Vec* var_values) {
  Heur* hh;
  if (hlist && net) {
    for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh))
      HEUR_apply_to_branch(hh,clist,net,br,var_values);
  }
}

void HEUR_list_clear(Heur* hlist, Net* net) {
  Heur* hh;
  if (hlist && net) {
    for (hh = hlist; hh != NULL; hh = HEUR_get_next(hh))
      HEUR_clear(hh,net);
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

Heur* HEUR_new(int type, Net* net) {
  Heur* h = (Heur*)malloc(sizeof(Heur));

  // Fields
  h->type = type;
  h->bus_counted = NULL;
  h->data = NULL;
  h->next = NULL;
  
  // Functions
  if (type == HEUR_TYPE_PVPQ) { // PV-PQ switching heuristic
    h->func_init = HEUR_PVPQ_init;
    h->func_clear = HEUR_PVPQ_clear;
    h->func_apply_to_branch = HEUR_PVPQ_apply_to_branch;
    h->func_free = HEUR_PVPQ_free;
  }
  else {
    h->func_init = NULL;
    h->func_clear = NULL;
    h->func_apply_to_branch = NULL;
    h->func_free = NULL;
  }

  // Type-specific init
  if (h->func_init)
    (*(h->func_init))(h,net);
  
  return h;
}

void HEUR_set_bus_counted(Heur* h, char* counted) {
  if (h)
    h->bus_counted = counted;
}

void HEUR_set_data(Heur* h, void* data) {
  if (h)
    h->data = data;
}

void HEUR_clear(Heur* h, Net* net) {
  if (h && h->func_clear)
    (*(h->func_clear))(h,net);
}

void HEUR_apply_to_branch(Heur* h, Constr* clist, Net* net, Branch* br, Vec* var_values) {
  if (h && h->func_apply_to_branch)
    (*(h->func_apply_to_branch))(h,clist,net,br,var_values);
}
