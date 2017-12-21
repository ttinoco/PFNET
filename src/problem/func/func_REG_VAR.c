/** @file func_REG_VAR.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VAR.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/func_REG_VAR.h>

struct Func_REG_VAR_Data {
  
  REAL* x0; // center
  REAL* w;  // weights
};

Func* FUNC_REG_VAR_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_init(f,&FUNC_REG_VAR_init);
  FUNC_set_func_count_step(f,&FUNC_REG_VAR_count_step);
  FUNC_set_func_allocate(f,&FUNC_REG_VAR_allocate);
  FUNC_set_func_clear(f,&FUNC_REG_VAR_clear);
  FUNC_set_func_analyze_step(f,&FUNC_REG_VAR_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_VAR_eval_step);
  FUNC_set_func_free(f,&FUNC_REG_VAR_free);
  FUNC_set_func_set_parameter(f,&FUNC_REG_VAR_set_parameter);
  FUNC_init(f);
  return f;
}

void FUNC_REG_VAR_init(Func* f) {

  // Local variables
  Net* net;
  int num_vars;
  Func_REG_VAR_Data* data;

  // Init
  net = FUNC_get_network(f);
  num_vars = NET_get_num_vars(net);
  data = (Func_REG_VAR_Data*)malloc(sizeof(Func_REG_VAR_Data));
  ARRAY_zalloc(data->x0,REAL,num_vars);
  ARRAY_zalloc(data->w,REAL,num_vars);
  FUNC_set_name(f,"variable regularization");
  FUNC_set_data(f,(void*)data);
}

void FUNC_REG_VAR_clear(Func* f) {

  // phi
  FUNC_set_phi(f,0);

  // gphi
  VEC_set_zero(FUNC_get_gphi(f));

  // Hphi
  // Constant so not clear it

  // Counter
  FUNC_set_Hphi_nnz(f,0);

  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_REG_VAR_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  Bat* bat;
  Load* load;
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

  // Tap ratio
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO))
    (*Hphi_nnz)++;

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE))
    (*Hphi_nnz)++;
  
  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // Voltage magnitude
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	(*Hphi_nnz)++;

      // Voltage angle
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG))
	(*Hphi_nnz)++;

      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

	// Active power
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P))
	  (*Hphi_nnz)++;

	// Reactive power
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q))
	  (*Hphi_nnz)++;
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {

	// Active power
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P))
	  (*Hphi_nnz)++;

	// Reactive power
	if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q))
	  (*Hphi_nnz)++;
	
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

	// Susceptance
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC))
	  (*Hphi_nnz)++;
	  
      }

      // Batteries
      for (bat = BUS_get_bat(bus); bat != NULL; bat = BAT_get_next(bat)) {

	// Charging/discharging power
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_P))
	  (*Hphi_nnz)++;
	  

	// Energy level
	if (BAT_has_flags(bat,FLAG_VARS,BAT_VAR_E))
	  (*Hphi_nnz)++;
      }

      // Loads
      for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

	// Active power
	if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_P))
	  (*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VAR_allocate(Func* f) {

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

void FUNC_REG_VAR_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv = 1.;
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

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	MAT_set_i(H,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
	MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	(*Hphi_nnz)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VAR_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_v_mag;
  REAL v;
  REAL vt;
  REAL dv = 1.;
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
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      // Set point
      vt = BUS_get_v_set(bus,t);

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	// Index
	index_v_mag = BUS_get_index_v_mag(bus,t);

	// v
	v = VEC_get(var_values,index_v_mag);

	// phi
	(*phi) += 0.5*pow((v-vt)/dv,2.);

	// gphi
	gphi[index_v_mag] = (v-vt)/(dv*dv);
      }
      else {

	// v
	v = BUS_get_v_mag(bus,t);

	// phi
	(*phi) += 0.5*pow((v-vt)/dv,2.);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VAR_free(Func* f) {

  // Local variables
  Func_REG_VAR_Data* data;

  // Get data
  data = (Func_REG_VAR_Data*)FUNC_get_data(f);

  // Free
  if (data) {
    free(data->x0);
    free(data->w);
    free(data);
  }

  // Set data
  FUNC_set_data(f,NULL);
}

void FUNC_REG_VAR_set_parameter(Func* f, char* key, void* value) {


}
