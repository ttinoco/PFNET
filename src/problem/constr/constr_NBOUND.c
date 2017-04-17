/** @file constr_NBOUND.c
 *  @brief This file defines the data structure and routines associated with the constraint of type NBOUND.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_NBOUND.h>

Constr* CONSTR_NBOUND_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_init(c, &CONSTR_NBOUND_init);
  CONSTR_set_func_count_step(c, &CONSTR_NBOUND_count_step);
  CONSTR_set_func_allocate(c, &CONSTR_NBOUND_allocate);
  CONSTR_set_func_clear(c, &CONSTR_NBOUND_clear);
  CONSTR_set_func_analyze_step(c, &CONSTR_NBOUND_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_NBOUND_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_NBOUND_store_sens_step);
  CONSTR_set_func_free(c, &CONSTR_NBOUND_free);
  CONSTR_init(c);
  return c;
}

void CONSTR_NBOUND_init(Constr* c) {

  // Init
  CONSTR_set_H_nnz(c,NULL,0);
  CONSTR_set_name(c,"variable nonlinear bounds");
  CONSTR_set_data(c,NULL);
}

void CONSTR_NBOUND_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));

  // Counters
  CONSTR_set_J_nnz(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_NBOUND_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int* J_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*J_nnz)++; // upper bound
    (*J_nnz)++; // lower bound
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*J_nnz)++; // upper bound
    (*J_nnz)++; // lower bound
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	(*J_nnz)++; // upper bound
	(*J_nnz)++; // lower bound
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*J_nnz)++; // upper bound
	(*J_nnz)++; // lower bound
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Shunt suscepatnace (b)
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_NBOUND_allocate(Constr* c) {

  // Local variables
  int J_nnz;
  Mat* H_array;
  Mat* H;
  int num_vars;
  int i;

  J_nnz = CONSTR_get_J_nnz(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));

  // A b
  CONSTR_set_A(c,MAT_new(0,num_vars,0));
  CONSTR_set_b(c,VEC_new(0));

  // f
  CONSTR_set_f(c,VEC_new(J_nnz));

  // J
  CONSTR_set_J(c,MAT_new(J_nnz,    // size1 (rows)
			 num_vars, // size2 (cols)
			 J_nnz));  // nnz

  // H
  H_array = MAT_array_new(J_nnz);
  CONSTR_set_H_array(c,H_array,J_nnz);
  for (i = 0; i < J_nnz; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,1);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(1,sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(1,sizeof(int)));
    MAT_set_data_array(H,(REAL*)malloc(1*sizeof(REAL)));
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,   // size1 (rows)
				  num_vars,   // size2 (cols)
				  J_nnz)); // nnz
}

void CONSTR_NBOUND_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  Mat* J;
  Mat* H_array;
  Mat* H;
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int* J_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int index_var;
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {

    index_var = BRANCH_get_index_ratio(br,t);

    // J
    MAT_set_i(J,*J_nnz,*J_nnz);
    MAT_set_j(J,*J_nnz,index_var);

    MAT_set_i(J,*J_nnz+1,*J_nnz+1);
    MAT_set_j(J,*J_nnz+1,index_var);

    // H
    H = MAT_array_get(H_array,*J_nnz);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    H = MAT_array_get(H_array,*J_nnz+1);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    (*J_nnz)++;
    (*J_nnz)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {

    index_var = BRANCH_get_index_phase(br,t);

    // J
    MAT_set_i(J,*J_nnz,*J_nnz);
    MAT_set_j(J,*J_nnz,index_var);

    MAT_set_i(J,*J_nnz+1,*J_nnz+1);
    MAT_set_j(J,*J_nnz+1,index_var);

    // H
    H = MAT_array_get(H_array,*J_nnz);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    H = MAT_array_get(H_array,*J_nnz+1);
    MAT_set_i(H,0,index_var);
    MAT_set_j(H,0,index_var);

    (*J_nnz)++;
    (*J_nnz)++;
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

	index_var = BUS_get_index_v_mag(bus,t);

	// J
	MAT_set_i(J,*J_nnz,*J_nnz);
	MAT_set_j(J,*J_nnz,index_var);

	MAT_set_i(J,*J_nnz+1,*J_nnz+1);
	MAT_set_j(J,*J_nnz+1,index_var);

	// H
	H = MAT_array_get(H_array,*J_nnz);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	H = MAT_array_get(H_array,*J_nnz+1);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	(*J_nnz)++;
	(*J_nnz)++;
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) && BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {

	index_var = BUS_get_index_v_ang(bus,t);

	// J
	MAT_set_i(J,*J_nnz,*J_nnz);
	MAT_set_j(J,*J_nnz,index_var);

	MAT_set_i(J,*J_nnz+1,*J_nnz+1);
	MAT_set_j(J,*J_nnz+1,index_var);

	// H
	H = MAT_array_get(H_array,*J_nnz);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	H = MAT_array_get(H_array,*J_nnz+1);
	MAT_set_i(H,0,index_var);
	MAT_set_j(H,0,index_var);

	(*J_nnz)++;
	(*J_nnz)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {

	  index_var = GEN_get_index_P(gen,t);

	  // J
	  MAT_set_i(J,*J_nnz,*J_nnz);
	  MAT_set_j(J,*J_nnz,index_var);

	  MAT_set_i(J,*J_nnz+1,*J_nnz+1);
	  MAT_set_j(J,*J_nnz+1,index_var);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  (*J_nnz)++;
	  (*J_nnz)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {

	  index_var = GEN_get_index_Q(gen,t);

	  // J
	  MAT_set_i(J,*J_nnz,*J_nnz);
	  MAT_set_j(J,*J_nnz,index_var);

	  MAT_set_i(J,*J_nnz+1,*J_nnz+1);
	  MAT_set_j(J,*J_nnz+1,index_var);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  (*J_nnz)++;
	  (*J_nnz)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {

	  index_var = SHUNT_get_index_b(shunt,t);

	  // J
	  MAT_set_i(J,*J_nnz,*J_nnz);
	  MAT_set_j(J,*J_nnz,index_var);

	  MAT_set_i(J,*J_nnz+1,*J_nnz+1);
	  MAT_set_j(J,*J_nnz+1,index_var);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_i(H,0,index_var);
	  MAT_set_j(H,0,index_var);

	  (*J_nnz)++;
	  (*J_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }

  // Done (last branch and period)
  if ((t == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {

    // Ensure lower triangular and save struct of H comb
    Hi_comb = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hj_comb = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (k = 0; k < CONSTR_get_H_array_size(c); k++) {
      Hi = MAT_get_row_array(MAT_array_get(H_array,k));
      Hj = MAT_get_col_array(MAT_array_get(H_array,k));
      Hi_comb[k] = Hi[0];
      Hj_comb[k] = Hj[0];
    }
  }
}

void CONSTR_NBOUND_eval_step(Constr* c, Branch* br, int t, Vec* values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  Mat* H_array;
  REAL* f;
  REAL* J;
  Mat* H;
  int* J_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL u;
  REAL umin;
  REAL umax;
  REAL du;
  REAL a1;
  REAL a2;
  REAL b;
  REAL eps;
  REAL sqrterm1;
  REAL sqrterm2;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!f || !J || !J_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Param
  eps = CONSTR_NBOUND_PARAM;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {

    u = VEC_get(values,BRANCH_get_index_ratio(br,t));
    umax = BRANCH_get_ratio_max(br);
    umin = BRANCH_get_ratio_min(br);
    du = (umax-umin > eps) ? umax-umin : eps;

    a1 = umax-u;
    a2 = u-umin;
    b = eps*eps/du;
    sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
    sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

    // f
    f[*J_nnz]   = a1 + b - sqrterm1; // upper
    f[*J_nnz+1] = a2 + b - sqrterm2; // lower

    // J
    J[*J_nnz]   = -(1-a1/sqrterm1);
    J[*J_nnz+1] = (1-a2/sqrterm2);

    // H
    H = MAT_array_get(H_array,*J_nnz);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

    H = MAT_array_get(H_array,*J_nnz+1);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

    (*J_nnz)++;
    (*J_nnz)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {

    u = VEC_get(values,BRANCH_get_index_phase(br,t));
    umax = BRANCH_get_phase_max(br);
    umin = BRANCH_get_phase_min(br);
    du = (umax-umin > eps) ? umax-umin : eps;

    a1 = umax-u;
    a2 = u-umin;
    b = eps*eps/du;
    sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
    sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

    // f
    f[*J_nnz]   = a1 + b - sqrterm1; // upper
    f[*J_nnz+1] = a2 + b - sqrterm2; // lower

    // J
    J[*J_nnz]   = -(1-a1/sqrterm1);
    J[*J_nnz+1] = (1-a2/sqrterm2);

    // H
    H = MAT_array_get(H_array,*J_nnz);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

    H = MAT_array_get(H_array,*J_nnz+1);
    MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

    (*J_nnz)++;
    (*J_nnz)++;
  }

  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

	u = VEC_get(values,BUS_get_index_v_mag(bus,t));
	umax = BUS_get_v_max_reg(bus);
	umin = BUS_get_v_min_reg(bus);
	du = (umax-umin > eps) ? umax-umin : eps;

	a1 = umax-u;
	a2 = u-umin;
	b = eps*eps/du;
	sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

	// f
	f[*J_nnz]   = a1 + b - sqrterm1; // upper
	f[*J_nnz+1] = a2 + b - sqrterm2; // lower

	// J
	J[*J_nnz]   = -(1-a1/sqrterm1);
	J[*J_nnz+1] = (1-a2/sqrterm2);

	// H
	H = MAT_array_get(H_array,*J_nnz);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

	H = MAT_array_get(H_array,*J_nnz+1);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

	(*J_nnz)++;
	(*J_nnz)++;
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {

	u = VEC_get(values,BUS_get_index_v_ang(bus,t));
	umax = 2*PI;
	umin = -2*PI;
	du = (umax-umin > eps) ? umax-umin : eps;

	a1 = umax-u;
	a2 = u-umin;
	b = eps*eps/du;
	sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

	// f
	f[*J_nnz]   = a1 + b - sqrterm1; // upper
	f[*J_nnz+1] = a2 + b - sqrterm2; // lower

	// J
	J[*J_nnz]   = -(1-a1/sqrterm1);
	J[*J_nnz+1] = (1-a2/sqrterm2);

	// H
	H = MAT_array_get(H_array,*J_nnz);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

	H = MAT_array_get(H_array,*J_nnz+1);
	MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

	(*J_nnz)++;
	(*J_nnz)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {

	  u = VEC_get(values,GEN_get_index_P(gen,t));
	  umax = GEN_get_P_max(gen);
	  umin = GEN_get_P_min(gen);
	  du = (umax-umin > eps) ? umax-umin : eps;

	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

	  // f
	  f[*J_nnz]   = a1 + b - sqrterm1; // upper
	  f[*J_nnz+1] = a2 + b - sqrterm2; // lower

	  // J
	  J[*J_nnz]   = -(1-a1/sqrterm1);
	  J[*J_nnz+1] = (1-a2/sqrterm2);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

	  (*J_nnz)++;
	  (*J_nnz)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {

	  u = VEC_get(values,GEN_get_index_Q(gen,t));
	  umax = GEN_get_Q_max(gen);
	  umin = GEN_get_Q_min(gen);
	  du = (umax-umin > eps) ? umax-umin : eps;

	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

	  // f
	  f[*J_nnz]   = a1 + b - sqrterm1; // upper
	  f[*J_nnz+1] = a2 + b - sqrterm2; // lower

	  // J
	  J[*J_nnz]   = -(1-a1/sqrterm1);
	  J[*J_nnz+1] = (1-a2/sqrterm2);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

	  (*J_nnz)++;
	  (*J_nnz)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {

	  u = VEC_get(values,SHUNT_get_index_b(shunt,t));
	  umax = SHUNT_get_b_max(shunt);
	  umin = SHUNT_get_b_min(shunt);
	  du = (umax-umin > eps) ? umax-umin : eps;

	  a1 = umax-u;
	  a2 = u-umin;
	  b = eps*eps/du;
	  sqrterm1 = sqrt(a1*a1+b*b+eps*eps);
	  sqrterm2 = sqrt(a2*a2+b*b+eps*eps);

	  // f
	  f[*J_nnz]   = a1 + b - sqrterm1; // upper
	  f[*J_nnz+1] = a2 + b - sqrterm2; // lower

	  // J
	  J[*J_nnz]   = -(1-a1/sqrterm1);
	  J[*J_nnz+1] = (1-a2/sqrterm2);

	  // H
	  H = MAT_array_get(H_array,*J_nnz);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm1*sqrterm1*sqrterm1));

	  H = MAT_array_get(H_array,*J_nnz+1);
	  MAT_set_d(H,0,-(b*b+eps*eps)/(sqrterm2*sqrterm2*sqrterm2));

	  (*J_nnz)++;
	  (*J_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_NBOUND_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int* J_nnz;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  J_nnz = CONSTR_get_J_nnz_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!J_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*J_nnz)++; // upper bound
    (*J_nnz)++; // lower bound
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_BOUNDED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*J_nnz)++; // upper bound
    (*J_nnz)++; // lower bound
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	BUS_set_sens_v_mag_u_bound(bus,VEC_get(sf,*J_nnz),t);
	(*J_nnz)++; // upper bound
	BUS_set_sens_v_mag_l_bound(bus,VEC_get(sf,*J_nnz),t);
	(*J_nnz)++; // lower bound
      }

      // Volage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_BOUNDED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*J_nnz)++; // upper bound
	(*J_nnz)++; // lower bound
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_BOUNDED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_BOUNDED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*J_nnz)++; // upper bound
	  (*J_nnz)++; // lower bound
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_NBOUND_free(Constr* c) {
  // Nothing
}
