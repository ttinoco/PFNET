/** @file func_REG_VANG.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VANG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_VANG.h>

void FUNC_REG_VANG_init(Func* f) {
  // Nothing
}

void FUNC_REG_VANG_clear(Func* f) {

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

void FUNC_REG_VANG_count_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  int bus_index[2];
  Bus* bus;
  int* Hcounter;
  char* bus_counted;
  int k;

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
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Branch
  if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
    (*Hcounter)++; // wk and wk

    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG))
      (*Hcounter)++; // wk and wm
  }
  if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG)) // wm var
    (*Hcounter)++; // wm and wm

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) // w var
	(*Hcounter)++; // w and w
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_VANG_allocate(Func* f) {

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

void FUNC_REG_VANG_analyze_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  int bus_index[2];
  int index_v_ang[2];
  Bus* bus;
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dw = FUNC_REG_VANG_PARAM;

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
  buses[0] = BRANCH_get_bus_k(br);
  buses[1] = BRANCH_get_bus_m(br);
  for (k = 0; k < 2; k++) {
    bus_index[k] = BUS_get_index(buses[k]);
    index_v_ang[k] = BUS_get_index_v_ang(buses[k]);
  }

  // Branch
  if (BUS_has_flags(buses[0],FLAG_VARS,BUS_VAR_VANG)) { // wk var
    MAT_set_i(H,*Hcounter,index_v_ang[0]);
    MAT_set_j(H,*Hcounter,index_v_ang[0]);
    MAT_set_d(H,*Hcounter,1./(dw*dw));
    (*Hcounter)++; // wk and wk

    if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG)) {
      if (index_v_ang[0] >= index_v_ang[1]) {
	MAT_set_i(H,*Hcounter,index_v_ang[0]);
	MAT_set_j(H,*Hcounter,index_v_ang[1]);
      }
      else {
	MAT_set_i(H,*Hcounter,index_v_ang[1]);
	MAT_set_j(H,*Hcounter,index_v_ang[0]);
      }
      MAT_set_d(H,*Hcounter,-1./(dw*dw));
      (*Hcounter)++; // wk and wm
    }
  }
  if (BUS_has_flags(buses[1],FLAG_VARS,BUS_VAR_VANG)) { // wm var
    MAT_set_i(H,*Hcounter,index_v_ang[1]);
    MAT_set_j(H,*Hcounter,index_v_ang[1]);
    MAT_set_d(H,*Hcounter,1./(dw*dw));
    (*Hcounter)++; // wm and wm
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) { // v var
	MAT_set_i(H,*Hcounter,index_v_ang[k]);
	MAT_set_j(H,*Hcounter,index_v_ang[k]);
	MAT_set_d(H,*Hcounter,1./(dw*dw));
	(*Hcounter)++;
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_VANG_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  int bus_index[2];
  int index_v_ang[2];
  REAL w[2];
  BOOL var_w[2];
  Bus* bus;
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL shift;
  int k;
  REAL wdiff;
  REAL dw = FUNC_REG_VANG_PARAM;

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
  for (k = 0; k < 2; k++) {
    bus_index[k] = BUS_get_index(buses[k]);
    index_v_ang[k] = BUS_get_index_v_ang(buses[k]);
    var_w[k] = BUS_has_flags(buses[k],FLAG_VARS,BUS_VAR_VANG);
    if (var_w[k])
      w[k] = VEC_get(var_values,index_v_ang[k]);
    else
      w[k] = BUS_get_v_ang(buses[k]);
  }

  // Branch data
  shift = BRANCH_get_phase(br); // radians

  // Branch
  if (var_w[0] || var_w[1]) { // wk or wm var

    wdiff = (w[0]-w[1]-shift);

    // phi
    (*phi) += 0.5*pow(wdiff/dw,2.);

    // gphi
    if (var_w[0]) // wk var
      gphi[index_v_ang[0]] += wdiff/(dw*dw);
    if (var_w[1]) // wm var
      gphi[index_v_ang[1]] -= wdiff/(dw*dw);
  }

  // Buses
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) { // v var

	// phi
	(*phi) += 0.5*pow(w[k]/dw,2.);

	// gphi
	gphi[index_v_ang[k]] += w[k]/(dw*dw);
      }
    }

    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_VANG_free(Func* f) {
  // Nothing
}
