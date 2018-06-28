/** @file func_REG_PQ.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_PQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_PQ.h>

Func* FUNC_REG_PQ_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f, &FUNC_REG_PQ_count_step);
  FUNC_set_func_analyze_step(f, &FUNC_REG_PQ_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_REG_PQ_eval_step);
  FUNC_set_name(f,"generator powers regularization");
  return f;
}

void FUNC_REG_PQ_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted)
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

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) // Q var
	  (*Hphi_nnz)++;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) // P var
	  (*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_PQ_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  H = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!Hphi_nnz || !bus_counted)
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

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var

	  dv = GEN_get_Q_max(gen)-GEN_get_Q_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hphi_nnz,GEN_get_index_Q(gen,t));
	  MAT_set_j(H,*Hphi_nnz,GEN_get_index_Q(gen,t));
	  MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	  (*Hphi_nnz)++;
	}

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var

	  dv = GEN_get_P_max(gen)-GEN_get_P_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  MAT_set_j(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	  (*Hphi_nnz)++;
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

	// Outage
	if (GEN_is_on_outage(gen))
	  continue;

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
