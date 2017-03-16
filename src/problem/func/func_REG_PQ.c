/** @file func_REG_PQ.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_PQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_PQ.h>

Func* FUNC_REG_PQ_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_name(f,"generator powers regularization");
  FUNC_set_func_init(f, &FUNC_REG_PQ_init);
  FUNC_set_func_count_step(f, &FUNC_REG_PQ_count_step);
  FUNC_set_func_allocate(f, &FUNC_REG_PQ_allocate);
  FUNC_set_func_clear(f, &FUNC_REG_PQ_clear);
  FUNC_set_func_analyze_step(f, &FUNC_REG_PQ_analyze_step);
  FUNC_set_func_eval_setp(f, &FUNC_REG_PQ_eval_step);
  FUNC_set_func_free(f, &FUNC_REG_PQ_free);
  return f;
}

void FUNC_REG_PQ_init(Func* f) {
  // Nothing
}

void FUNC_REG_PQ_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);

  // gphi
  VEC_set_zero(FUNC_get_gphi(f));

  // Hphi
  // Constant so not clear it

  // Counter
  FUNC_set_Hcounter(f,0);

  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_REG_PQ_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index_t[2];
  int* Hcounter;
  char* bus_counted;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(bus[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) // Q var
	  (*Hcounter)++;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) // P var
	  (*Hcounter)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_PQ_allocate(Func* f) {

  // Local variables
  int num_vars;
  int Hcounter;

  num_vars = NET_get_num_vars(FUNC_get_network(f));
  Hcounter = FUNC_get_Hcounter(f);

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));

  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  Hcounter));
}

void FUNC_REG_PQ_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index_t[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  H = FUNC_get_Hphi(f);
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(bus[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var

	  dv = GEN_get_Q_max(gen)-GEN_get_Q_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hcounter,GEN_get_index_Q(gen,t));
	  MAT_set_j(H,*Hcounter,GEN_get_index_Q(gen,t));
	  MAT_set_d(H,*Hcounter,1./(dv*dv));
	  (*Hcounter)++;
	}

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var

	  dv = GEN_get_P_max(gen)-GEN_get_P_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hcounter,GEN_get_index_P(gen,t));
	  MAT_set_j(H,*Hcounter,GEN_get_index_P(gen,t));
	  MAT_set_d(H,*Hcounter,1./(dv*dv));
	  (*Hcounter)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_PQ_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL Qmid;
  REAL Pmid;
  REAL P;
  REAL Q;
  REAL dP;
  REAL dQ;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !gphi || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_k(br);
  bus[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(bus[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index_t[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {

	// Mid value
	Qmid = (GEN_get_Q_max(gen)+GEN_get_Q_min(gen))/2.; // p.u.
	Pmid = (GEN_get_P_max(gen)+GEN_get_P_min(gen))/2.; // p.u.

	// Normalization factor
	dQ = GEN_get_Q_max(gen)-GEN_get_Q_min(gen); // p.u.
	if (dQ < FUNC_REG_PQ_PARAM)
	  dQ = FUNC_REG_PQ_PARAM;
	dP = GEN_get_P_max(gen)-GEN_get_P_min(gen); // p.u.
	if (dP < FUNC_REG_PQ_PARAM)
	  dP = FUNC_REG_PQ_PARAM;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var

	  // Value
	  Q = VEC_get(var_values,GEN_get_index_Q(gen,t));

	  // phi
	  (*phi) += 0.5*pow((Q-Qmid)/dQ,2.);

	  // gphi
	  gphi[GEN_get_index_Q(gen,t)] = (Q-Qmid)/(dQ*dQ);
	}
	else {

	  // Value
	  Q = GEN_get_Q(gen,t);

	  // phi
	  (*phi) += 0.5*pow((Q-Qmid)/dQ,2.);
	}

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var

	  // Value
	  P = VEC_get(var_values,GEN_get_index_P(gen,t));

	  // phi
	  (*phi) += 0.5*pow((P-Pmid)/dP,2.);

	  // gphi
	  gphi[GEN_get_index_P(gen,t)] = (P-Pmid)/(dP*dP);
	}
	else {

	  // Value
	  P = GEN_get_P(gen,t);

	  // phi
	  (*phi) += 0.5*pow((P-Pmid)/dP,2.);
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_PQ_free(Func* f) {
  // Nothing
}
