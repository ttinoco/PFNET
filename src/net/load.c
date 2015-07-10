/** @file load.c
 *  @brief This file defines the Load data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/load.h>
#include <pfnet/bus.h>

struct Load {

  // Bus
  Bus *bus;           /**< @brief Bus to which the load is connected */

  // Active power
  REAL P;             /**< @brief Load active power (p.u. system base power) */

  // Reactive power
  REAL Q;             /**< @brief Load reactive power (p.u. system base power) */

  // Indices
  int index;          /**< @brief Load index */

  // List
  struct Load* next;  /**< @brief List of loads connected to a bus */

};

void* LOAD_array_get(void* load, int index) { 
  if (load) 
    return (void*)&(((Load*)load)[index]);
  else
    return NULL;
}

Load* LOAD_array_new(int num) { 
  int i;
  Load* load = (Load*)malloc(sizeof(Load)*num);
  for (i = 0; i < num; i++) {
    LOAD_init(&(load[i]));
    LOAD_set_index(&(load[i]),i);
  }
  return load;
}

void LOAD_array_show(Load* load, int num) { 
  int i;
  for (i = 0; i < num; i++) 
    LOAD_show(&(load[i]));
}

void* LOAD_get_bus(Load* load) {
  if (load)
    return (void*)load->bus;
  else
    return NULL;
}

int LOAD_get_index(Load* load) {
  if (load)
    return load->index;
  else
    return 0;
}

Load* LOAD_get_next(Load* load) {
  if (!load)
    return NULL;
  else
    return load->next;
}

REAL LOAD_get_P(Load* load) {
  if (!load)
    return 0;
  else
    return load->P;
}

REAL LOAD_get_Q(Load* load) {
  if (!load)
    return 0;
  else
    return load->Q;
}

void LOAD_init(Load* load) { 
  load->bus = NULL;
  load->P = 0;
  load->Q = 0;
  load->index = 0;
  load->next = NULL;
}

Load* LOAD_list_add(Load *load_list, Load* load) {
  LIST_add(load_list,load,next);
  return load_list;
}

int LOAD_list_len(Load* load_list) {
  int len;
  LIST_len(Load,load_list,next,len);
  return len;
}

Load* LOAD_new(void) { 
  Load* load = (Load*)malloc(sizeof(Load));
  LOAD_init(load);
  return load;
}

void LOAD_set_bus(Load* load, void* bus) { 
  load->bus = (Bus*)bus;
}

void LOAD_set_index(Load* load, int index) { 
  load->index = index;
}

void LOAD_set_P(Load* load, REAL P) { 
  load->P = P;
}

void LOAD_set_Q(Load* load, REAL Q) { 
  load->Q = Q;
}

void LOAD_show(Load* load) { 
  printf("load %d\t%d\n",BUS_get_number(load->bus),load->index);
}

