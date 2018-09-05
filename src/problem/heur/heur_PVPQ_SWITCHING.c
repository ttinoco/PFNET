/** @file heur_PVPQ_SWITCHING.c
 *  @brief This file defines the data structure and routines associated with the heuristic of type PVPQ_SWITCHING.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/reg_obj.h>
#include <pfnet/heur_PVPQ_SWITCHING.h>
#include <pfnet/constr_PVPQ_SWITCHING.h>

Heur* HEUR_PVPQ_SWITCHING_new(Net* net) {
  Heur* h = HEUR_new(net);
  HEUR_set_func_apply_step(h,&HEUR_PVPQ_SWITCHING_apply_step);
  HEUR_set_name(h,"PVPQ switching");
  HEUR_init(h);
  return h;
}

void HEUR_PVPQ_SWITCHING_apply_step(Heur* h, Constr** cptrs, int cnum, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Net* net;
  Vec* f;
  char* fix_flag;
  char all_fixed;
  Constr* pf;
  Constr* pvpq;
  REAL v;
  REAL v_set;
  REAL Q;
  REAL Qmax;
  REAL Qmin;
  int i;

  char obj_type;
  void* obj;

  // Heur data
  net = HEUR_get_network(h);
  
  // Power flow constraint
  pf = NULL;
  for (i = 0; i < cnum; i++) {
    if (strcmp(CONSTR_get_name(cptrs[i]),"AC power balance") == 0) {
      pf = cptrs[i];
      break;
    }
  }
  if (!pf) {
    HEUR_set_error(h, "unable to find AC power balance constraint");
    return;
  }

  // PVPQ switching constraint
  pvpq = NULL;
  for (i = 0; i < cnum; i++) {
    if (strcmp(CONSTR_get_name(cptrs[i]),"PVPQ switching") == 0) {
      pvpq = cptrs[i];
      break;
    }
  }
  if (!pvpq) {
    HEUR_set_error(h, "unable to find PVPQ switching constraint");
    return;
  }  

  // Fix flags
  fix_flag = CONSTR_PVPQ_SWITCHING_get_flags(pvpq);
  if (!fix_flag) {
    HEUR_set_error(h, "unable to get PVPQ switching constraint flags");
    return;
  }

  // Constr data
  f = CONSTR_get_f(pf);

  // Candidate bus
  if (BUS_is_v_set_regulated(bus) &&                 // regulated
      BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {   // v mag is variable
      
    // Voltage magnitude and set point
    v = VEC_get(var_values,BUS_get_index_v_mag(bus,t));
    v_set = BUS_get_v_set(bus,t);
      
    // CASE: v fixed
    if (fix_flag[BUS_get_index_v_mag(bus,t)]) {
        
      all_fixed = TRUE;
        
      // Regulating objects
      for (REG_OBJ_init(&obj_type,&obj,bus); obj != NULL; REG_OBJ_next(&obj_type,&obj,bus)) {
          
        // Check candidacy
        if (!REG_OBJ_is_candidate(obj_type,obj))
          continue;
          
        Q = VEC_get(var_values,REG_OBJ_get_index_Q(obj_type,obj,t)); // per unit
        Qmax = REG_OBJ_get_Q_max(obj_type,obj);                      // per unit
        Qmin = REG_OBJ_get_Q_min(obj_type,obj);                      // per unit
          
        // Q is fixed
        if (fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)]) {
            
          // Q at Qmin and v < v_set - check if free Q might help
          if (fabs(Q-Qmin) < fabs(Q-Qmax) && v < v_set) {
            Q = Q - VEC_get(f,BUS_get_index_Q(REG_OBJ_get_bus(obj_type,obj),t)); // per unit (see constr_PF)
            if (Qmin < Q) {
              fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = FALSE; // Switch to Q free
              all_fixed = FALSE;
              if (HEUR_PVPQ_SWITCHING_DEBUG)
                printf("HEUR PVPQ: free reg obj from Qmin (reg bus %d fixed)\n", BUS_get_number(bus));
            }
          }
            
          // Q at Qmax and v > v_set - check if free Q might help
          else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
            Q = Q - VEC_get(f,BUS_get_index_Q(REG_OBJ_get_bus(obj_type,obj),t)); // per unit (see constr_PF)
            if (Q < Qmax) {
              fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = FALSE; // Switch to Q free
              all_fixed = FALSE;
              if (HEUR_PVPQ_SWITCHING_DEBUG)
                printf("HEUR PVPQ: free reg obj from Qmax (reg bus %d fixed)\n", BUS_get_number(bus));
            }
          }	    
        }	    
          
        // Q is free: Qmax violation
        else if (Q >= Qmax) {
            
          // Switch to Q fixed
          REG_OBJ_set_Q(obj_type,obj,Qmax,t);
          fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = TRUE;
          if (HEUR_PVPQ_SWITCHING_DEBUG)
            printf("HEUR PVPQ: fix reg obj at Qmax (reg bus %d fixed)\n", BUS_get_number(bus));
        }
          
        // Q is free: Qmin violation
        else if (Q <= Qmin) {
            
          // Switch to Q fixed
          REG_OBJ_set_Q(obj_type,obj,Qmin,t);
          fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = TRUE;
          if (HEUR_PVPQ_SWITCHING_DEBUG)
            printf("HEUR PVPQ: fix reg obj at Qmin (reg bus %d fixed)\n", BUS_get_number(bus));
        }
          
        // Q is free: no violation
        else
          all_fixed = FALSE;
      }
        
      // All reg objs fixed
      if (all_fixed) {
          
        // Switch to v free
        fix_flag[BUS_get_index_v_mag(bus,t)] = FALSE;
        if (HEUR_PVPQ_SWITCHING_DEBUG)
          printf("HEUR PVPQ: free reg bus %d\n", BUS_get_number(bus));
      }
    }
      
    // CASE: v free
    else {
        
      all_fixed = TRUE;
        
      // Regulating objects
      for (REG_OBJ_init(&obj_type,&obj,bus); obj != NULL; REG_OBJ_next(&obj_type,&obj,bus)) {
          
        // Check candidacy
        if (!REG_OBJ_is_candidate(obj_type,obj))
          continue;
          
        // Q is free (should never happen here)
        if (!fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)]) {
          fprintf(stderr, "WARNING: PVPQ switching has both v and Q free\n");
          all_fixed = FALSE;
          continue;
        }
          
        // Reg obj data
        Q = VEC_get(var_values,REG_OBJ_get_index_Q(obj_type,obj,t)); // per unit
        Qmax = REG_OBJ_get_Q_max(obj_type,obj);                      // per unit
        Qmin = REG_OBJ_get_Q_min(obj_type,obj);                      // per unit
          
        // Q at Qmin and v < v_set - check if free Q might help
        if (fabs(Q-Qmin) < fabs(Q-Qmax) && v < v_set) {
          Q = Q - VEC_get(f,BUS_get_index_Q(REG_OBJ_get_bus(obj_type,obj),t)); // per unit (see constr_PF)
          if (Qmin < Q) {
            fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = FALSE; // Switch to Q free
            all_fixed = FALSE;
            if (HEUR_PVPQ_SWITCHING_DEBUG)
              printf("HEUR PVPQ: free reg obj from Qmin (reg bus %d free)\n", BUS_get_number(bus));		
          }
        }
          
        // Q at Qmax and v > v_set - check if free Q might help
        else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
          Q = Q - VEC_get(f,BUS_get_index_Q(REG_OBJ_get_bus(obj_type,obj),t)); // per unit (see constr_PF)
          if (Q < Qmax) {
            fix_flag[REG_OBJ_get_index_Q(obj_type,obj,t)] = FALSE; // Switch to Q free
            all_fixed = FALSE;
            if (HEUR_PVPQ_SWITCHING_DEBUG)
              printf("HEUR PVPQ: free reg obj from Qmax (reg bus %d free)\n", BUS_get_number(bus));
          }
        }
      }
        
      // Not all reg objs are fixed
      if (!all_fixed) {
          
        // Switch to v fixed
        fix_flag[BUS_get_index_v_mag(bus,t)] = TRUE;
        if (HEUR_PVPQ_SWITCHING_DEBUG)
          printf("HEUR PVPQ: fix reg bus %d\n", BUS_get_number(bus));
      }
    }
  }
      
  // Update
  if (BUS_get_index(bus) == NET_get_num_buses(net)-1 && t == NET_get_num_periods(net)-1)
    CONSTR_analyze(pvpq);
}
