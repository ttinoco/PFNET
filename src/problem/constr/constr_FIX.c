/** @file constr_FIX.c
 *  @brief This file defines the data structure and routines associated with the constraint of type FIX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_FIX.h>

Constr* CONSTR_FIX_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c, &CONSTR_FIX_count_step);
  CONSTR_set_func_analyze_step(c, &CONSTR_FIX_analyze_step);
  CONSTR_set_func_eval_step(c, &CONSTR_FIX_eval_step);
  CONSTR_set_func_store_sens_step(c, &CONSTR_FIX_store_sens_step);
  CONSTR_set_name(c,"variable fixing");
  return c;
}

void CONSTR_FIX_count_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Load* load;
  Bat* bat;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*A_nnz)++;
    (*A_row)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*A_nnz)++;
    (*A_row)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	(*A_nnz)++;

	// Extra nz for regulating generator (for PV-PQ switching without changing nnz pattern)
	if (BUS_is_regulated_by_gen(bus)) {
	  for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	      (*A_nnz)++;
	  }
	}

	(*A_row)++;
      }

      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*A_nnz)++;
	(*A_row)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

	// Charging/discharging power (P)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_P) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
	  (*A_nnz) += 2;
	  (*A_row) += 2;
	}

	// Energy level (E)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_E) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Active power (P)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_P) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_Q) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_FIX_analyze_step(Constr* c, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Bat* bat;
  Load* load;
  Shunt* shunt;
  int* A_nnz;
  int* A_row;
  char* bus_counted;
  Vec* b;
  Mat* A;
  int i;
  REAL Pc;
  REAL Pd;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !A_row || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    VEC_set(b,*A_row,BRANCH_get_ratio(br,t));
    MAT_set_i(A,*A_nnz,*A_row);
    MAT_set_j(A,*A_nnz,BRANCH_get_index_ratio(br,t));
    MAT_set_d(A,*A_nnz,1.);
    (*A_nnz)++;
    (*A_row)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    VEC_set(b,*A_row,BRANCH_get_phase(br,t));
    MAT_set_i(A,*A_nnz,*A_row);
    MAT_set_j(A,*A_nnz,BRANCH_get_index_phase(br,t));
    MAT_set_d(A,*A_nnz,1.);
    (*A_nnz)++;
    (*A_row)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	if (BUS_is_regulated_by_gen(bus))
	  VEC_set(b,*A_row,BUS_get_v_set(bus,t));
	else
	  VEC_set(b,*A_row,BUS_get_v_mag(bus,t));
	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++;

	// Extra nz for regulating generator (for PV-PQ switching?)
	if (BUS_is_regulated_by_gen(bus)) {
	  for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
	    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	      MAT_set_i(A,*A_nnz,*A_row);
	      MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen,t));
	      MAT_set_d(A,*A_nnz,0.); // placeholder
	      (*A_nnz)++;
	    }
	  }
	}

	(*A_row)++;
      }

      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	VEC_set(b,*A_row,BUS_get_v_ang(bus,t));
	MAT_set_i(A,*A_nnz,*A_row);
	MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus,t));
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++;
	(*A_row)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  VEC_set(b,*A_row,GEN_get_P(gen,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,GEN_get_index_P(gen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  VEC_set(b,*A_row,GEN_get_Q(gen,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  VEC_set(b,*A_row,VARGEN_get_P(vargen,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,VARGEN_get_index_P(vargen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  VEC_set(b,*A_row,VARGEN_get_Q(vargen,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,VARGEN_get_index_Q(vargen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  VEC_set(b,*A_row,SHUNT_get_b(shunt,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,SHUNT_get_index_b(shunt,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

	// Charging/discharging power (P)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_P) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {

	  if (BAT_get_P(bat,t) >= 0) {
	    Pc = BAT_get_P(bat,t);
	    Pd = 0.;
	  }
	  else {
	    Pc = 0.;
	    Pd = -BAT_get_P(bat,t);
	  }

	  // Pc
	  VEC_set(b,*A_row,Pc);
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pc(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;

	  // Pd
	  VEC_set(b,*A_row,Pd);
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pd(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Energy level (E)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_E) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E)) {
	  VEC_set(b,*A_row,BAT_get_E(bat,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,BAT_get_index_E(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Active power (P)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_P) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  VEC_set(b,*A_row,LOAD_get_P(load,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}

	// Reactive power (Q)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_Q) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {
	  VEC_set(b,*A_row,LOAD_get_Q(load,t));
	  MAT_set_i(A,*A_nnz,*A_row);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_Q(load,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*A_row)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_FIX_eval_step(Constr* c, Branch* br, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_FIX_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}
