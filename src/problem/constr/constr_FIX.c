/** @file constr_FIX.c
 *  @brief This file defines the data structure and routines associated with the constraint of type FIX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_FIX.h>

void CONSTR_FIX_init(Constr* c) {

  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_FIX_clear(Constr* c) {

  // Counters
  CONSTR_set_A_nnz(c,0);
  CONSTR_set_Aconstr_index(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
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
  int* Aconstr_index;
  char* bus_counted;
  int i;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !Aconstr_index || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*A_nnz)++;
    (*Aconstr_index)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*A_nnz)++;
    (*Aconstr_index)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	(*A_nnz)++;

	// Extra nz for regulating generator (for PV-PQ switching?)
	if (BUS_is_regulated_by_gen(bus) &&
	    GEN_has_flags(BUS_get_reg_gen(bus),FLAG_VARS,GEN_VAR_Q))
	  (*A_nnz)++;

	(*Aconstr_index)++;
      }

      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*A_nnz)++;
	(*Aconstr_index)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Batteries
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

	// Charging/discharging power (P)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_P) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P)) {
	  (*A_nnz) += 2;
	  (*Aconstr_index) += 2;
	}

	// Energy level (E)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_E) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Active power (P)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_P) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_FIX_allocate(Constr* c) {

  // Local variables
  int num_constr;
  int num_vars;
  int A_nnz;

  // Constr data
  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_constr = CONSTR_get_Aconstr_index(c);
  A_nnz = CONSTR_get_A_nnz(c);

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // G l u
  CONSTR_set_G(c,MAT_new(0,num_vars,0));
  CONSTR_set_l(c,VEC_new(0));
  CONSTR_set_u(c,VEC_new(0));

  // b
  CONSTR_set_b(c,VEC_new(num_constr));

  // A
  CONSTR_set_A(c,MAT_new(num_constr, // size1 (rows)
			 num_vars,   // size2 (cols)
			 A_nnz)); // nnz
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
  int* Aconstr_index;
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
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!A_nnz || !Aconstr_index || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    VEC_set(b,*Aconstr_index,BRANCH_get_ratio(br,t));
    MAT_set_i(A,*A_nnz,*Aconstr_index);
    MAT_set_j(A,*A_nnz,BRANCH_get_index_ratio(br,t));
    MAT_set_d(A,*A_nnz,1.);
    (*A_nnz)++;
    (*Aconstr_index)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    VEC_set(b,*Aconstr_index,BRANCH_get_phase(br,t));
    MAT_set_i(A,*A_nnz,*Aconstr_index);
    MAT_set_j(A,*A_nnz,BRANCH_get_index_phase(br,t));
    MAT_set_d(A,*A_nnz,1.);
    (*A_nnz)++;
    (*Aconstr_index)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {

    bus = buses[i];

    if (!bus_counted[BUS_get_index(bus)*T+t]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	if (BUS_is_regulated_by_gen(bus))
	  VEC_set(b,*Aconstr_index,BUS_get_v_set(bus,t));
	else
	  VEC_set(b,*Aconstr_index,BUS_get_v_mag(bus,t));
	MAT_set_i(A,*A_nnz,*Aconstr_index);
	MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++;

	// Extra nz for regulating generator (for PV-PQ switching?)
	if (BUS_is_regulated_by_gen(bus) &&
	    GEN_has_flags(BUS_get_reg_gen(bus),FLAG_VARS,GEN_VAR_Q)) {
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,GEN_get_index_Q(BUS_get_reg_gen(bus),t));
	  MAT_set_d(A,*A_nnz,0.); // placeholder
	  (*A_nnz)++;
	}

	(*Aconstr_index)++;
      }

      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	VEC_set(b,*Aconstr_index,BUS_get_v_ang(bus,t));
	MAT_set_i(A,*A_nnz,*Aconstr_index);
	MAT_set_j(A,*A_nnz,BUS_get_index_v_ang(bus,t));
	MAT_set_d(A,*A_nnz,1.);
	(*A_nnz)++;
	(*Aconstr_index)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  VEC_set(b,*Aconstr_index,GEN_get_P(gen,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,GEN_get_index_P(gen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) &&
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  VEC_set(b,*Aconstr_index,GEN_get_Q(gen,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  VEC_set(b,*Aconstr_index,VARGEN_get_P(vargen,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,VARGEN_get_index_P(vargen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) &&
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  VEC_set(b,*Aconstr_index,VARGEN_get_Q(vargen,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,VARGEN_get_index_Q(vargen,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  VEC_set(b,*Aconstr_index,SHUNT_get_b(shunt,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,SHUNT_get_index_b(shunt,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
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
	  VEC_set(b,*Aconstr_index,Pc);
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pc(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;

	  // Pd
	  VEC_set(b,*Aconstr_index,Pd);
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BAT_get_index_Pd(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}

	// Energy level (E)
	if (BAT_has_flags(bat,FLAG_FIXED,BAT_VAR_E) &&
	    BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E)) {
	  VEC_set(b,*Aconstr_index,BAT_get_E(bat,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,BAT_get_index_E(bat,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Active power (P)
	if (LOAD_has_flags(load,FLAG_FIXED,LOAD_VAR_P) &&
	    LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P)) {
	  VEC_set(b,*Aconstr_index,LOAD_get_P(load,t));
	  MAT_set_i(A,*A_nnz,*Aconstr_index);
	  MAT_set_j(A,*A_nnz,LOAD_get_index_P(load,t));
	  MAT_set_d(A,*A_nnz,1.);
	  (*A_nnz)++;
	  (*Aconstr_index)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)*T+t] = TRUE;
  }
}

void CONSTR_FIX_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {
  // Nothing to do
}

void CONSTR_FIX_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_FIX_free(Constr* c) {
  // Nothing to do
}
