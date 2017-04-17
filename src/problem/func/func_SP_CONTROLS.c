/** @file func_SP_CONTROLS.c
 *  @brief This file defines the data structure and routines associated with the function of type SP_CONTROLS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_SP_CONTROLS.h>

Func* FUNC_SP_CONTROLS_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_init(f, &FUNC_SP_CONTROLS_init);
  FUNC_set_func_count_step(f, &FUNC_SP_CONTROLS_count_step);
  FUNC_set_func_allocate(f, &FUNC_SP_CONTROLS_allocate);
  FUNC_set_func_clear(f, &FUNC_SP_CONTROLS_clear);
  FUNC_set_func_analyze_step(f, &FUNC_SP_CONTROLS_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_SP_CONTROLS_eval_step);
  FUNC_set_func_free(f, &FUNC_SP_CONTROLS_free);
  FUNC_init(f);
  return f;
}

void FUNC_SP_CONTROLS_init(Func* f) {
  
  FUNC_set_name(f,"sparse controls penalty");
}

void FUNC_SP_CONTROLS_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);

  // gphi
  VEC_set_zero(FUNC_get_gphi(f));

  // Hphi
  MAT_set_zero_d(FUNC_get_Hphi(f));

  // Counter
  FUNC_set_Hphi_nnz(f,0);

  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_SP_CONTROLS_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
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

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Tap ratio of tap-changing transformer
  if (BRANCH_is_tap_changer(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_RATIO))
    (*Hphi_nnz)++;

  // Phase shift of phase-shifting transformer
  if (BRANCH_is_phase_shifter(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_PHASE))
    (*Hphi_nnz)++;

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // Voltage mag of gen-regulated bus
      if (BUS_is_regulated_by_gen(bus) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_SPARSE,BUS_VAR_VMAG))
	(*Hphi_nnz)++;

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_SPARSE,GEN_VAR_P))
	  (*Hphi_nnz)++;
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance of switched shunt device
	if (SHUNT_is_switched(shunt) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_SPARSE,SHUNT_VAR_SUSC))
	  (*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_SP_CONTROLS_allocate(Func* f) {

  // Local variables
  int num_vars;
  int Hphi_nnz;

  num_vars = NET_get_num_vars(FUNC_get_network(f));
  Hphi_nnz = FUNC_get_Hphi_nnz(f);

  // gphi
  FUNC_set_gphi(f,VEC_new(num_vars));

  // Hphi
  FUNC_set_Hphi(f,MAT_new(num_vars,
			  num_vars,
			  Hphi_nnz));
}

void FUNC_SP_CONTROLS_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;
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

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Tap ratio of tap-changing transformer
  if (BRANCH_is_tap_changer(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_RATIO)) {
    MAT_set_i(H,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
    MAT_set_j(H,*Hphi_nnz,BRANCH_get_index_ratio(br,t));
    (*Hphi_nnz)++;
  }

  // Phase shift of phase-shifting transformer
  if (BRANCH_is_phase_shifter(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_PHASE)) {
    MAT_set_i(H,*Hphi_nnz,BRANCH_get_index_phase(br,t));
    MAT_set_j(H,*Hphi_nnz,BRANCH_get_index_phase(br,t));
    (*Hphi_nnz)++;
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];
    
    if (!bus_counted[bus_index_t[k]]) {

      // Voltage mag of gen-regulated bus
      if (BUS_is_regulated_by_gen(bus) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_SPARSE,BUS_VAR_VMAG)) {
	MAT_set_i(H,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	(*Hphi_nnz)++;
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_SPARSE,GEN_VAR_P)) {
	  MAT_set_i(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  MAT_set_j(H,*Hphi_nnz,GEN_get_index_P(gen,t));
	  (*Hphi_nnz)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance of switched shunt device
	if (SHUNT_is_switched(shunt) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_SPARSE,SHUNT_VAR_SUSC)) {
	  MAT_set_i(H,*Hphi_nnz,SHUNT_get_index_b(shunt,t));
	  MAT_set_j(H,*Hphi_nnz,SHUNT_get_index_b(shunt,t));
	  (*Hphi_nnz)++;
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_SP_CONTROLS_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Shunt* shunt;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL* Hd;
  int index_val;
  REAL val;
  REAL val0;
  REAL dval;
  REAL sqrt_term;
  int k;
  int T;

  // Num periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hd = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);

  // Check pointers
  if (!phi || !gphi || !Hd || !Hphi_nnz || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Tap ratio of tap-changing transformer
  if (BRANCH_is_tap_changer(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_RATIO)) {

    index_val = BRANCH_get_index_ratio(br,t);
    val = VEC_get(var_values,index_val);
    val0 = BRANCH_get_ratio(br,t);
    dval = BRANCH_get_ratio_max(br)-BRANCH_get_ratio_min(br);
    if (dval < FUNC_SP_CONTROLS_CEPS)
      dval = FUNC_SP_CONTROLS_CEPS;
    sqrt_term = sqrt( (val-val0)*(val-val0)/(dval*dval) + FUNC_SP_CONTROLS_EPS );

    // phi
    (*phi) += sqrt_term;

    // gphi
    gphi[index_val] = (1./sqrt_term)*((val-val0)/(dval*dval));

    // Hphi
    Hd[*Hphi_nnz] = FUNC_SP_CONTROLS_EPS/(dval*dval*sqrt_term*sqrt_term*sqrt_term);
    (*Hphi_nnz)++;
  }
  else {
    // nothing because val0-val0 = 0 for constant val
  }

  // Phase shift of phase-shifting transformer
  if (BRANCH_is_phase_shifter(br) &&
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE) &&
      BRANCH_has_flags(br,FLAG_SPARSE,BRANCH_VAR_PHASE)) {

    index_val = BRANCH_get_index_phase(br,t);
    val = VEC_get(var_values,index_val);
    val0 = BRANCH_get_phase(br,t);
    dval = BRANCH_get_phase_max(br)-BRANCH_get_phase_min(br);
    if (dval < FUNC_SP_CONTROLS_CEPS)
      dval = FUNC_SP_CONTROLS_CEPS;
    sqrt_term = sqrt( (val-val0)*(val-val0)/(dval*dval) + FUNC_SP_CONTROLS_EPS );

    // phi
    (*phi) += sqrt_term;

    // gphi
    gphi[index_val] = (1./sqrt_term)*((val-val0)/(dval*dval));

    // Hphi
    Hd[*Hphi_nnz] = FUNC_SP_CONTROLS_EPS/(dval*dval*sqrt_term*sqrt_term*sqrt_term);
    (*Hphi_nnz)++;
  }
  else {
    // nothing because val0-val0 = 0 for constant val
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // Voltage mag of gen-regulated bus
      if (BUS_is_regulated_by_gen(bus) &&
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	  BUS_has_flags(bus,FLAG_SPARSE,BUS_VAR_VMAG)) {

	index_val = BUS_get_index_v_mag(bus,t);
	val = VEC_get(var_values,index_val);
	val0 = BUS_get_v_set(bus,t);
	dval = BUS_get_v_max_reg(bus)-BUS_get_v_min_reg(bus);
	if (dval < FUNC_SP_CONTROLS_CEPS)
	  dval = FUNC_SP_CONTROLS_CEPS;
	sqrt_term = sqrt( (val-val0)*(val-val0)/(dval*dval) + FUNC_SP_CONTROLS_EPS );

	// phi
	(*phi) += sqrt_term;

	// gphi
	gphi[index_val] = (1./sqrt_term)*((val-val0)/(dval*dval));

	// Hphi
	Hd[*Hphi_nnz] = FUNC_SP_CONTROLS_EPS/(dval*dval*sqrt_term*sqrt_term*sqrt_term);
	(*Hphi_nnz)++;
      }
      else {
	// nothing because val0-val0 = 0 for constant val
      }

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) &&
	    GEN_has_flags(gen,FLAG_SPARSE,GEN_VAR_P)) {

	  index_val = GEN_get_index_P(gen,t);
	  val = VEC_get(var_values,index_val);
	  val0 = GEN_get_P(gen,t);
	  dval = GEN_get_P_max(gen)-GEN_get_P_min(gen);
	  if (dval < FUNC_SP_CONTROLS_CEPS)
	    dval = FUNC_SP_CONTROLS_CEPS;
	  sqrt_term = sqrt( (val-val0)*(val-val0)/(dval*dval) + FUNC_SP_CONTROLS_EPS );

	  // phi
	  (*phi) += sqrt_term;

	  // gphi
	  gphi[index_val] = (1./sqrt_term)*((val-val0)/(dval*dval));

	  // Hphi
	  Hd[*Hphi_nnz] = FUNC_SP_CONTROLS_EPS/(dval*dval*sqrt_term*sqrt_term*sqrt_term);
	  (*Hphi_nnz)++;
	}
	else {
	  // nothing because val0-val0 = 0 for constant val
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance of switched shunt device
	if (SHUNT_is_switched(shunt) &&
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	    SHUNT_has_flags(shunt,FLAG_SPARSE,SHUNT_VAR_SUSC)) {

	  index_val = SHUNT_get_index_b(shunt,t);
	  val = VEC_get(var_values,index_val);
	  val0 = SHUNT_get_b(shunt,t);
	  dval = SHUNT_get_b_max(shunt)-SHUNT_get_b_min(shunt);
	  if (dval < FUNC_SP_CONTROLS_CEPS)
	    dval = FUNC_SP_CONTROLS_CEPS;
	  sqrt_term = sqrt( (val-val0)*(val-val0)/(dval*dval) + FUNC_SP_CONTROLS_EPS );

	  // phi
	  (*phi) += sqrt_term;

	  // gphi
	  gphi[index_val] = (1./sqrt_term)*((val-val0)/(dval*dval));

	  // Hphi
	  Hd[*Hphi_nnz] = FUNC_SP_CONTROLS_EPS/(dval*dval*sqrt_term*sqrt_term*sqrt_term);
	  (*Hphi_nnz)++;
	}
	else {
	  // nothing because val0-val0 = 0 for constant val
	}
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_SP_CONTROLS_free(Func* f) {
  // Nothing
}
