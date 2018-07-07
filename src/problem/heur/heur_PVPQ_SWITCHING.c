/** @file heur_PVPQ_SWITCHING.c
 *  @brief This file defines the data structure and routines associated with the heuristic of type PVPQ_SWITCHING.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/heur_PVPQ_SWITCHING.h>
#include <pfnet/constr_PVPQ_SWITCHING.h>

Heur* HEUR_PVPQ_SWITCHING_new(Net* net) {
  Heur* h = HEUR_new(net);
  HEUR_set_func_apply_step(h,&HEUR_PVPQ_SWITCHING_apply_step);
  HEUR_set_name(h,"PVPQ switching");
  HEUR_init(h);
  return h;
}

void HEUR_PVPQ_SWITCHING_apply_step(Heur* h, Constr** cptrs, int cnum, Bus* bus, int t, Vec* var_values) {

  // Local variables
  Net* net;
  Vec* f;
  Gen* gen;
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
    if (HEUR_PVPQ_SWITCHING_DEBUG)
      printf("HEUR PVPQ SWITCHING: no AC power balance constraint\n");
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
    if (HEUR_PVPQ_SWITCHING_DEBUG)
      printf("HEUR PVPQ SWITCHING: no PVPQ switching constraint\n");
    return;
  }

  // Fix flags
  fix_flag = CONSTR_PVPQ_SWITCHING_get_flags(pvpq);
  if (!fix_flag) {
    HEUR_set_error(h, "unable to get PVPQ switching constraint flags");
    if (HEUR_PVPQ_SWITCHING_DEBUG)
      printf("HEUR PVPQ SWITCHING: no PVPQ switching constraint flags\n");
    return;
  }

  // Constr data
  f = CONSTR_get_f(pf);

  // Candidate bus
  if (!BUS_is_slack(bus) &&                          // not slack
      BUS_is_regulated_by_gen(bus) &&                // regulated
      BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {   // v mag is variable
      
    // Voltage magnitude and set point
    v = VEC_get(var_values,BUS_get_index_v_mag(bus,t));
    v_set = BUS_get_v_set(bus,t);
    
    // CASE: v fixed
    if (fix_flag[BUS_get_index_v_mag(bus,t)]) {
      
      all_fixed = TRUE;
      
      // Generators
      for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	
	// Q not variable
	if (!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) // reg gen Q is not variable
	  continue;
	
	Q = VEC_get(var_values,GEN_get_index_Q(gen,t)); // per unit
	Qmax = GEN_get_Q_max(gen);                      // per unit
	Qmin = GEN_get_Q_min(gen);                      // per unit
	
	// Q is fixed
	if (fix_flag[GEN_get_index_Q(gen,t)]) {
	  
	  // Q at Qmin and v < v_set - check if free Q might help
	  if (fabs(Q-Qmin) < fabs(Q-Qmax) && v < v_set) {
	    Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen),t)); // per unit (see constr_PF)
	    if (Qmin < Q) {
	      fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	      all_fixed = FALSE;
	      if (HEUR_PVPQ_SWITCHING_DEBUG)
		printf("HEUR PVPQ SWITCHING: free gen %d from Qmin (reg bus %d fixed)\n", GEN_get_index(gen), BUS_get_number(bus));
	    }
	  }
	  
	  // Q at Qmax and v > v_set - check if free Q might help
	  else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
	    Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen),t)); // per unit (see constr_PF)
	    if (Q < Qmax) {
	      fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	      all_fixed = FALSE;
	      if (HEUR_PVPQ_SWITCHING_DEBUG)
		printf("HEUR PVPQ SWITCHING: free gen %d from Qmax (reg bus %d fixed)\n", GEN_get_index(gen), BUS_get_number(bus));
	    }
	  }	    
	}	    
	
	// Q is free: Qmax violation
	else if (Q >= Qmax) {
	  
	  // Switch to Q fixed
	  GEN_set_Q(gen,Qmax,t);
	  fix_flag[GEN_get_index_Q(gen,t)] = TRUE;
	  if (HEUR_PVPQ_SWITCHING_DEBUG)
	    printf("HEUR PVPQ SWITCHING: fix gen %d at Qmax (reg bus %d fixed)\n", GEN_get_index(gen), BUS_get_number(bus));
	}
	
	// Q is free: Qmin violation
	else if (Q <= Qmin) {
	  
	  // Switch to Q fixed
	  GEN_set_Q(gen,Qmin,t);
	  fix_flag[GEN_get_index_Q(gen,t)] = TRUE;
	  if (HEUR_PVPQ_SWITCHING_DEBUG)
	    printf("HEUR PVPQ SWITCHING: fix gen %d at Qmin (reg bus %d fixed)\n", GEN_get_index(gen), BUS_get_number(bus));
	}
	
	// Q is free: no violation
	else
	  all_fixed = FALSE;
      }
      
      // All gens fixed
      if (all_fixed) {
	
	// Switch to v free
	fix_flag[BUS_get_index_v_mag(bus,t)] = FALSE;
	if (HEUR_PVPQ_SWITCHING_DEBUG)
	  printf("HEUR PVPQ SWITCHING: free reg bus %d\n", BUS_get_number(bus));
      }
    }
    
    // CASE: v free
    else {
      
      all_fixed = TRUE;
      
      // Generators
      for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	
	// Q not variable
	if (!GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) // reg gen Q is not variable
	  continue;
	
	// Q is free (should never happen here)
	if (!fix_flag[GEN_get_index_Q(gen,t)]) {
	  fprintf(stderr, "WARNING: PVPQ switching has both v and Q free\n");
	  all_fixed = FALSE;
	  continue;
	}
	
	// Reg gen data
	Q = VEC_get(var_values,GEN_get_index_Q(gen,t)); // per unit
	Qmax = GEN_get_Q_max(gen);                      // per unit
	Qmin = GEN_get_Q_min(gen);                      // per unit
	
	// Q at Qmin and v < v_set - check if free Q might help
	if (fabs(Q-Qmin) < fabs(Q-Qmax) && v < v_set) {
	  Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen),t)); // per unit (see constr_PF)
	  if (Qmin < Q) {
	    fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	    all_fixed = FALSE;
	    if (HEUR_PVPQ_SWITCHING_DEBUG)
	      printf("HEUR PVPQ SWITCHING: free gen %d from Qmin (reg bus %d free)\n", GEN_get_index(gen), BUS_get_number(bus));		
	  }
	}
	
	// Q at Qmax and v > v_set - check if free Q might help
	else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
	  Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen),t)); // per unit (see constr_PF)
	  if (Q < Qmax) {
	    fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	    all_fixed = FALSE;
	    if (HEUR_PVPQ_SWITCHING_DEBUG)
	      printf("HEUR PVPQ SWITCHING: free gen %d from Qmax (reg bus %d free)\n", GEN_get_index(gen), BUS_get_number(bus));
	  }
	}
      }
      
      // Not all gens are fixed
      if (!all_fixed) {
	
	// Switch to v fixed
	fix_flag[BUS_get_index_v_mag(bus,t)] = TRUE;
	if (HEUR_PVPQ_SWITCHING_DEBUG)
	  printf("HEUR PVPQ SWITCHING: fix reg bus %d\n", BUS_get_number(bus));
      }
    }
  }
  
  // Update
  if (BUS_get_index(bus) == NET_get_num_buses(net)-1 && t == BUS_get_num_periods(bus)-1)
    CONSTR_analyze(pvpq);
}
