/** @file func_REG_VMAG.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VMAG.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_VMAG.h>

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
  FUNC_set_Hcounter(f,0);
  
  // Flags
  FUNC_clear_bus_counted(f);
}

void FUNC_REG_VMAG_count_branch(Func* f, Branch* br) {

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

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yv var
	(*Hcounter)++; // y var
	(*Hcounter)++; // z var
      }    

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var
	(*Hcounter)++; // vl
	(*Hcounter)++; // vh
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_allocate(Func* f) {
  
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

void FUNC_REG_VMAG_analyze_branch(Func* f, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  Mat* H;
  int k;
  REAL dv = FUNC_REG_VMAG_PARAM;

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
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
	MAT_set_i(H,*Hcounter,BUS_get_index_v_mag(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_v_mag(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++;
      }
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var
	MAT_set_i(H,*Hcounter,BUS_get_index_y(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_y(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++; // y var

	MAT_set_i(H,*Hcounter,BUS_get_index_z(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_z(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++; // z var
      }

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var
	MAT_set_i(H,*Hcounter,BUS_get_index_vl(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_vl(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++; // vl var

	MAT_set_i(H,*Hcounter,BUS_get_index_vh(bus));
	MAT_set_j(H,*Hcounter,BUS_get_index_vh(bus));
	MAT_set_d(H,*Hcounter,1./(dv*dv));
	(*Hcounter)++; // vz var
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }  
}

void FUNC_REG_VMAG_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  int bus_index[2];
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
      
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var

	// Index
	index_v_mag = BUS_get_index_v_mag(bus);

	// v
	v = VEC_get(var_values,index_v_mag);
	vt = BUS_get_v_set(bus);

	// phi
	(*phi) += 0.5*pow((v-vt)/dv,2.);

	// gphi
	gphi[index_v_mag] = (v-vt)/(dv*dv);
      }

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	// Indices
	index_y = BUS_get_index_y(bus);
	index_z = BUS_get_index_z(bus);

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

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	// Indices
	index_vl = BUS_get_index_vl(bus);
	index_vh = BUS_get_index_vh(bus);

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
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_VMAG_free(Func* f) {
  // Nothing
}
