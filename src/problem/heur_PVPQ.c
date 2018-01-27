/** @file heur_PVPQ.c
 *  @brief This file defines the data structure and routines associated with the heuristic of type PVPQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/heur_PVPQ.h>
#include <pfnet/constr_PVPQ_SWITCHING.h>

void HEUR_PVPQ_init(Heur* h, Net* net) {

  // Local variables
  int num_buses;
  int num_periods;

  // Init
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  HEUR_set_bus_counted(h,(char*)calloc(num_buses*num_periods,sizeof(char)));
  HEUR_set_data(h,NULL);
}

void HEUR_PVPQ_clear(Heur* h, Net* net) {

  // Local variables
  int num_buses;
  int num_periods;

  // Clear bus counted flags
  num_buses = NET_get_num_buses(net);
  num_periods = NET_get_num_periods(net);
  HEUR_clear_bus_counted(h,num_buses*num_periods);
}

void HEUR_PVPQ_apply_step(Heur* h, Constr* clist, Net* net, Branch* br, int t, Vec* var_values) {

  // Local variables
  Vec* f;
  Bus* bus[2];
  Gen* gen;
  char* bus_counted;
  int bus_index_t[2];
  char* fix_flag;
  char all_fixed;
  int k;
  Constr* pf;
  Constr* pvpq;
  REAL v;
  REAL v_set;
  REAL Q;
  REAL Qmax;
  REAL Qmin;
  int T;
  int num_buses;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Num buses
  num_buses = NET_get_num_buses(net);

  // Heur data
  bus_counted = HEUR_get_bus_counted(h);

  // Bus from data
  bus[0] = BRANCH_get_bus_k(br);
  bus_index_t[0] = BUS_get_index(bus[0])*T+t;

  // Bus to data
  bus[1] = BRANCH_get_bus_m(br);
  bus_index_t[1] = BUS_get_index(bus[1])*T+t;

  // DEBUG
  if (BRANCH_get_index(br) == 0 && t == 0) {
    printf("HEUR begin\n");
  }
  
  // Power flow constraint
  for (pf = clist; pf != NULL; pf = CONSTR_get_next(pf)) {
    if (strcmp(CONSTR_get_name(pf),"AC power balance") == 0)
      break;
  }
  if (!pf)
    return;

  // PVPQ switching constraint
  for (pvpq = clist; pvpq != NULL; pvpq = CONSTR_get_next(pvpq)) {
    if (strcmp(CONSTR_get_name(pvpq),"PVPQ switching") == 0)
      break;
  }
  if (!pvpq)
    return;

  // Fix flags
  fix_flag = CONSTR_PVPQ_SWITCHING_get_flags(pvpq);
  if (!fix_flag)
    return;

  // Constr data
  f = CONSTR_get_f(pf);

  // Buses
  for (k = 0; k < 2; k++) {

    // No bus (until I fix outages)
    if (!bus[k])
      continue;

    // Candidate bus
    if (!bus_counted[bus_index_t[k]] &&                   // not counted
	!BUS_is_slack(bus[k]) &&                          // not slack
	BUS_is_regulated_by_gen(bus[k]) &&                // regulated
	BUS_has_flags(bus[k],FLAG_VARS,BUS_VAR_VMAG)) {   // v mag is variable
      
      // Voltage magnitude and set point
      v = VEC_get(var_values,BUS_get_index_v_mag(bus[k],t));
      v_set = BUS_get_v_set(bus[k],t);

      // CASE: v fixed
      if (fix_flag[BUS_get_index_v_mag(bus[k],t)]) {

	all_fixed = TRUE;

	// Generators
	for (gen = BUS_get_reg_gen(bus[k]); gen != NULL; gen = GEN_get_reg_next(gen)) {
	  
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
	      Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen))+t*2*num_buses); // per unit (see constr_PF)
	      if (Qmin < Q) {
		fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
		all_fixed = FALSE;
	      }
	    }
	    
	    // Q at Qmax and v > v_set - check if free Q might help
	    else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
	      Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen))+t*2*num_buses); // per unit (see constr_PF)
	      if (Q < Qmax) {
	      fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	      all_fixed = FALSE;
	      }
	    }	    
	  }	    
	  
	  // Q is free: Qmax violation
	  else if (Q > Qmax) {

	    // Switch to Q fixed
	    GEN_set_Q(gen,Qmax,t);
	    fix_flag[GEN_get_index_Q(gen,t)] = TRUE;
	  }

	  // Q is free: Qmin violation
	  else if (Q < Qmin) {
	    
	    // Switch to Q fixed
	    GEN_set_Q(gen,Qmin,t);
	    fix_flag[GEN_get_index_Q(gen,t)] = TRUE;
	  }

	  // Q is free: no violation
	  else
	    all_fixed = FALSE;
	}
	
	// All gens fixed
	if (all_fixed) {

	  // Switch to v free
	  fix_flag[BUS_get_index_v_mag(bus[k],t)] = FALSE;
	}
      }

      // CASE: v free
      else {

	all_fixed = TRUE;

	// Generators
	for (gen = BUS_get_reg_gen(bus[k]); gen != NULL; gen = GEN_get_reg_next(gen)) {
	  
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
	    Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen))+t*2*num_buses); // per unit (see constr_PF)
	    if (Qmin < Q) {
	      fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	      all_fixed = FALSE;
	    }
	  }

	  // Q at Qmax and v > v_set - check if free Q might help
	  else if (fabs(Q-Qmax) < fabs(Q-Qmin) && v > v_set) { 
	    Q = Q - VEC_get(f,BUS_get_index_Q(GEN_get_bus(gen))+t*2*num_buses); // per unit (see constr_PF)
	    if (Q < Qmax) {
	      fix_flag[GEN_get_index_Q(gen,t)] = FALSE; // Switch to Q free
	      all_fixed = FALSE;
	    }
	  }
	}

	// Not all gens are fixed
	if (!all_fixed) {

	  // Switch to v fixed
	  fix_flag[BUS_get_index_v_mag(bus[k],t)] = TRUE;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }

  // Update
  if (BRANCH_get_index(br) == NET_get_num_branches(net)-1 && t == NET_get_num_periods(net)-1) {
    printf("HEUR done\n");
    CONSTR_analyze(pvpq);
  }
}

void HEUR_PVPQ_free(Heur* h) {
  // Nothing
}
