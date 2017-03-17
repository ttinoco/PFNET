/** @file func_REG_VMAG.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2017, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_VMAG.h>

Func* FUNC_REG_VMAG_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_name(f,"voltage magnitude regularization");
  FUNC_set_func_init(f, &FUNC_REG_VMAG_init);
  FUNC_set_func_count_step(f, &FUNC_REG_VMAG_count_step);
  FUNC_set_func_allocate(f, &FUNC_REG_VMAG_allocate);
  FUNC_set_func_clear(f, &FUNC_REG_VMAG_clear);
  FUNC_set_func_analyze_step(f, &FUNC_REG_VMAG_analyze_step);
  FUNC_set_func_eval_step(f, &FUNC_REG_VMAG_eval_step);
  FUNC_set_func_free(f, &FUNC_REG_VMAG_free);
  return f;
}

void FUNC_REG_VMAG_init(Func* f) {
  // Nothing
}

void FUNC_REG_VMAG_clear(Func* f) {

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

void FUNC_REG_VMAG_count_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
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

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) // v var
	(*Hphi_nnz)++;

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yv var
	(*Hphi_nnz)++; // y var
	(*Hphi_nnz)++; // z var
      }

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var
	(*Hphi_nnz)++; // vl
	(*Hphi_nnz)++; // vh
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_allocate(Func* f) {

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

void FUNC_REG_VMAG_analyze_step(Func* f, Branch* br, int t) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
  int* Hphi_nnz;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv = FUNC_REG_VMAG_PARAM;
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

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var
	MAT_set_i(H,*Hphi_nnz,BUS_get_index_y(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_y(bus,t));
	MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	(*Hphi_nnz)++; // y var

	MAT_set_i(H,*Hphi_nnz,BUS_get_index_z(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_z(bus,t));
	MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	(*Hphi_nnz)++; // z var
      }

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var
	MAT_set_i(H,*Hphi_nnz,BUS_get_index_vl(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_vl(bus,t));
	MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	(*Hphi_nnz)++; // vl var

	MAT_set_i(H,*Hphi_nnz,BUS_get_index_vh(bus,t));
	MAT_set_j(H,*Hphi_nnz,BUS_get_index_vh(bus,t));
	MAT_set_d(H,*Hphi_nnz,1./(dv*dv));
	(*Hphi_nnz)++; // vz var
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_eval_step(Func* f, Branch* br, int t, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index_t[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_v_mag;
  int index_y;
  int index_z;
  int index_vl;
  int index_vh;
  REAL v;
  REAL vt;
  REAL y;
  REAL z;
  REAL vl;
  REAL vh;
  REAL dv = FUNC_REG_VMAG_PARAM;
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

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	// Indices
	index_y = BUS_get_index_y(bus,t);
	index_z = BUS_get_index_z(bus,t);

	// y z
	y = VEC_get(var_values,index_y);
	z = VEC_get(var_values,index_z);

	// phi
	(*phi) += 0.5*pow(y/dv,2.); // y
	(*phi) += 0.5*pow(z/dv,2.); // z

	// gphi
	gphi[index_y] = y/(dv*dv);
	gphi[index_z] = z/(dv*dv);
      }
      else {
	// nothing
      }

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	// Indices
	index_vl = BUS_get_index_vl(bus,t);
	index_vh = BUS_get_index_vh(bus,t);

	// vl and vh
	vl = VEC_get(var_values,index_vl);
	vh = VEC_get(var_values,index_vh);

	// phi
	(*phi) += 0.5*pow(vl/dv,2.); // vl
	(*phi) += 0.5*pow(vh/dv,2.); // vh

	// gphi
	gphi[index_vl] = vl/(dv*dv);
	gphi[index_vh] = vh/(dv*dv);
      }
      else {
	// nothing
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_free(Func* f) {
  // Nothing
}
