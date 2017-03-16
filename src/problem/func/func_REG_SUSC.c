/** @file func_REG_SUSC.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_SUSC.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_SUSC.h>

Func* FUNC_REG_SUSC_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_name(f,"susceptance regularization");
  FUNC_set_func_init(f, &FUNC_REG_SUSC_init);
  FUNC_set_func_count_step(f, &FUNC_REG_SUSC_count_step);
  FUNC_set_func_allocate(f, &FUNC_REG_SUSC_allocate);
  FUNC_set_func_clear(f, &FUNC_REG_SUSC_clear);
  FUNC_set_func_analyze_step(f, &FUNC_REG_SUSC_analyze_step);
  FUNC_set_func_eval_setp(f, &FUNC_REG_SUSC_eval_step);
  FUNC_set_func_free(f, &FUNC_REG_SUSC_free);
  return f;
}

void FUNC_REG_SUSC_init(Func* f) {
  // Nothing
}

void FUNC_REG_SUSC_clear(Func* f) {

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

void FUNC_REG_SUSC_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Shunt* shunt;
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

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) // b var
	  (*Hcounter)++;

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var
	  (*Hcounter)++;
	  (*Hcounter)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_SUSC_allocate(Func* f) {

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

void FUNC_REG_SUSC_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* bus[2];
  Shunt* shunt;
  int bus_index_t[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  REAL db;
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

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt); // p.u.
	if (db < FUNC_REG_SUSC_PARAM)
	  db = FUNC_REG_SUSC_PARAM;

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  MAT_set_i(H,*Hcounter,SHUNT_get_index_b(shunt,t));
	  MAT_set_j(H,*Hcounter,SHUNT_get_index_b(shunt,t));
	  MAT_set_d(H,*Hcounter,1./(db*db));
	  (*Hcounter)++;
	}

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var

	  MAT_set_i(H,*Hcounter,SHUNT_get_index_y(shunt,t));
	  MAT_set_j(H,*Hcounter,SHUNT_get_index_y(shunt,t));
	  MAT_set_d(H,*Hcounter,1./(db*db));
	  (*Hcounter)++;

	  MAT_set_i(H,*Hcounter,SHUNT_get_index_z(shunt,t));
	  MAT_set_j(H,*Hcounter,SHUNT_get_index_z(shunt,t));
	  MAT_set_d(H,*Hcounter,1./(db*db));
	  (*Hcounter)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_SUSC_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* bus[2];
  Shunt* shunt;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL b0;
  REAL b;
  REAL db;
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

      // Shunts
      for (shunt = BUS_get_shunt(bus[k]); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Normalization factor
	db = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt); // p.u.
	if (db < FUNC_REG_SUSC_PARAM)
	  db = FUNC_REG_SUSC_PARAM;

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  b0 = SHUNT_get_b(shunt,t);
	  b = VEC_get(var_values,SHUNT_get_index_b(shunt,t));
	  (*phi) += 0.5*pow((b-b0)/db,2.);
	  gphi[SHUNT_get_index_b(shunt,t)] = (b-b0)/(db*db);
	}
	else {
	  // nothing because b0 - b0 = 0
	}

	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var

	  b = VEC_get(var_values,SHUNT_get_index_y(shunt,t));
	  (*phi) += 0.5*pow(b/db,2.);
	  gphi[SHUNT_get_index_y(shunt,t)] = b/(db*db);

	  b = VEC_get(var_values,SHUNT_get_index_z(shunt,t));
	  (*phi) += 0.5*pow(b/db,2.);
	  gphi[SHUNT_get_index_z(shunt,t)] = b/(db*db);
	}
	else {
	  // nothing because b0 - b0 = 0
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_SUSC_free(Func* f) {
  // Nothing
}
