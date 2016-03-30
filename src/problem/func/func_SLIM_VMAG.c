/** @file func_SLIM_VMAG.c
 *  @brief This file defines the data structure and routines associated with the function of type SLIM_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_SLIM_VMAG.h>

void FUNC_SLIM_VMAG_init(Func* f) {
  // Nothing
}

void FUNC_SLIM_VMAG_clear(Func* f) {
  
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

void FUNC_SLIM_VMAG_count_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index[2];
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
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) // v var
	(*Hcounter)++;
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_SLIM_VMAG_allocate(Func* f) {
  
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

void FUNC_SLIM_VMAG_analyze_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv;

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
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];
    
    if (!bus_counted[bus_index[k]]) {
      
      dv = BUS_get_v_max(bus)-BUS_get_v_min(bus);
      if (dv < FUNC_SLIM_VMAG_PARAM)
	dv = FUNC_SLIM_VMAG_PARAM;
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	MAT_set_i(H,*Hcounter,BUS_get_index_v_mag(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_v_mag(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++;
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }  
}

void FUNC_SLIM_VMAG_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  int index_v_mag;
  REAL v;
  REAL vmid;
  REAL dv;
  int k;

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
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];

    if (!bus_counted[bus_index[k]]) {

      dv = BUS_get_v_max(bus)-BUS_get_v_min(bus);
      if (dv < FUNC_SLIM_VMAG_PARAM)
	dv = FUNC_SLIM_VMAG_PARAM;

      vmid = 0.5*(BUS_get_v_max(bus)+BUS_get_v_min(bus));
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	// Index
	index_v_mag = BUS_get_index_v_mag(bus);
	
	// v
	v = VEC_get(var_values,index_v_mag);

	// phi
	(*phi) += 0.5*pow((v-vmid)/dv,2.);

	// gphi
	gphi[index_v_mag] = (v-vmid)/(dv*dv);
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_SLIM_VMAG_free(Func* f) {
  // Nothing
}
