/** @file reg_obj.c
 *  @brief This file defines routines for handling regulating objects.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/reg_obj.h>

void REG_OBJ_next(char* obj_type, void** obj, Bus* bus) {

  // Gen
  if (*obj_type == OBJ_GEN) {
    if (*obj && GEN_get_reg_next((Gen*)(*obj)))
      *obj = GEN_get_reg_next((Gen*)(*obj));
    else
      *obj = NULL;
    if (!(*obj)) { // Continue with VSC
      *obj_type = OBJ_CONVVSC;
      if (BUS_get_reg_vsc_conv(bus))
        *obj = BUS_get_reg_vsc_conv(bus);
    }
  }
  
  // VSC
  else if (*obj_type == OBJ_CONVVSC) {
    if (*obj && CONVVSC_get_reg_next((ConvVSC*)(*obj)))
      *obj = CONVVSC_get_reg_next((ConvVSC*)(*obj));
    else
      *obj = NULL;
    if (!(*obj)) { // Continue with FACTS
      *obj_type = OBJ_FACTS;
      if (BUS_get_reg_facts(bus))
        *obj = BUS_get_reg_facts(bus);
    }
  }

  // FACTS
  else if (*obj_type == OBJ_FACTS) {
    if (*obj && FACTS_get_reg_next((Facts*)(*obj)))
      *obj = FACTS_get_reg_next((Facts*)(*obj));
    else
      *obj = NULL;
  }

  // Unknown
  else {
    *obj = NULL;
    *obj_type = OBJ_UNKNOWN;
  }
}

void REG_OBJ_init(char* obj_type, void** obj, Bus* bus) {

  *obj_type = OBJ_GEN;
  if (BUS_get_reg_gen(bus))
    *obj = BUS_get_reg_gen(bus);
  else
    *obj = NULL;    
  
  if (!(*obj)) {
    *obj_type = OBJ_CONVVSC;
    if (BUS_get_reg_vsc_conv(bus))
      *obj = BUS_get_reg_vsc_conv(bus);
  }

  if (!(*obj)) {
    *obj_type = OBJ_FACTS;
    if (BUS_get_reg_facts(bus))
      *obj = BUS_get_reg_facts(bus);
  }
}

void REG_OBJ_set_Q(char obj_type, void* obj, REAL Q, int t) {

  // Check
  if (!obj)
    return;
  
  // Gen
  if (obj_type == OBJ_GEN)
    GEN_set_Q((Gen*)obj,Q,t);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    CONVVSC_set_Q((ConvVSC*)obj,Q,t);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    FACTS_set_Q_sh((Facts*)obj,Q,t);
}

void REG_OBJ_show(char obj_type, void* obj) {

  // Check
  if (!obj)
    return;
  
  // Gen
  if (obj_type == OBJ_GEN)
    printf("GEN %d\n", GEN_get_index((Gen*)obj));
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    printf("CONVVSC %d\n", CONVVSC_get_index((ConvVSC*)obj));

  // FACTS
  else if (obj_type == OBJ_FACTS)
    printf("FACTS %d\n", FACTS_get_index((Facts*)obj));

  // Unknown
  else
    printf("UNKONWN\n");
}

Bus* REG_OBJ_get_bus(char obj_type, void* obj) {

  // Check
  if (!obj)
    return NULL;
  
  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_bus((Gen*)obj);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_ac_bus((ConvVSC*)obj);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_bus_k((Facts*)obj);

  else
    return NULL;
}

int REG_OBJ_get_index_Q(char obj_type, void* obj, int t) {

  // Check
  if (!obj)
    return -1;

  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_index_Q((Gen*)obj,t);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_index_Q((ConvVSC*)obj,t);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_index_Q_sh((Facts*)obj,t);

  else
    return -1;
}

REAL REG_OBJ_get_Q(char obj_type, void* obj, int t) {

  // Check
  if (!obj)
    return 0.;
  
  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_Q((Gen*)obj,t);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_Q((ConvVSC*)obj,t);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_Q_sh((Facts*)obj,t);

  else
    return 0.;
}

REAL REG_OBJ_get_Q_max(char obj_type, void* obj) {

  // Check
  if (!obj)
    return 0.;

  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_Q_max((Gen*)obj);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_Q_max((ConvVSC*)obj);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_Q_max_sh((Facts*)obj);

  else
    return 0.;
}

REAL REG_OBJ_get_Q_min(char obj_type, void* obj) {

  // Check
  if (!obj)
    return 0.;

  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_Q_min((Gen*)obj);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_Q_min((ConvVSC*)obj);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_Q_min_sh((Facts*)obj);

  else
    return 0.;
}

REAL REG_OBJ_get_Q_par(char obj_type, void* obj) {

  // Check
  if (!obj)
    return 0.;

  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_get_Q_par((Gen*)obj);
      
  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_get_Q_par((ConvVSC*)obj);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_get_Q_par((Facts*)obj);

  else
    return 0.;
}

BOOL REG_OBJ_is_candidate(char obj_type, void* obj) {

  // Check
  if (!obj)
    return FALSE;

  // Gen
  if (obj_type == OBJ_GEN)
    return GEN_has_flags((Gen*)obj,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage((Gen*)obj);

  // VSC
  else if (obj_type == OBJ_CONVVSC)
    return CONVVSC_has_flags((ConvVSC*)obj,FLAG_VARS,CONVVSC_VAR_Q);

  // FACTS
  else if (obj_type == OBJ_FACTS)
    return FACTS_has_flags((Facts*)obj,FLAG_VARS,FACTS_VAR_Q);

  else
    return FALSE;
}

int REG_OBJ_count_candidates(Bus* bus) {
  char obj_type;
  void* obj;
  int num = 0;
  for(REG_OBJ_init(&obj_type,&obj,bus); obj != NULL; REG_OBJ_next(&obj_type,&obj,bus)) {
    if (REG_OBJ_is_candidate(obj_type,obj))
      num += 1;
  }
  return num;
}
