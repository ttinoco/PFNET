/** @file func_REG_PQ.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_PQ.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_REG_PQ.h>

void FUNC_REG_PQ_init(Func* f) {
  // Nothing
}

void FUNC_REG_PQ_clear(Func* f) {

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

void FUNC_REG_PQ_count_branch(Func* f, Branch *br) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index[2];
  int* Hcounter;
  char* bus_counted;
  int k;

  // Constr data
  Hcounter = FUNC_get_Hcounter_ptr(f);
  bus_counted = FUNC_get_bus_counted(f);
  if (!Hcounter || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(bus[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]]) {

      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) // Q var
	  (*Hcounter)++;

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) // P var
	  (*Hcounter)++;
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_PQ_allocate(Func* f) {
  
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

void FUNC_REG_PQ_analyze_branch(Func* f, Branch *br) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
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
  if (!Hcounter || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(bus[k]);

  // Buses
  for (k = 0; k < 2; k++) {
    
    if (!bus_counted[bus_index[k]]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var

	  dv = GEN_get_Q_max(gen)-GEN_get_Q_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hcounter,GEN_get_index_Q(gen));
	  MAT_set_j(H,*Hcounter,GEN_get_index_Q(gen));
	  MAT_set_d(H,*Hcounter,1./(dv*dv));
	  (*Hcounter)++;
	}

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
	  
	  dv = GEN_get_P_max(gen)-GEN_get_P_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  MAT_set_i(H,*Hcounter,GEN_get_index_P(gen));
	  MAT_set_j(H,*Hcounter,GEN_get_index_P(gen));
	  MAT_set_d(H,*Hcounter,1./(dv*dv));
	  (*Hcounter)++;
	}
      }  
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }  
}

void FUNC_REG_PQ_eval_branch(Func* f, Branch* br, Vec* var_values) {

  // Local variables
  Bus* bus[2];
  Gen* gen;
  int bus_index[2];
  char* bus_counted;
  REAL* phi;
  REAL* gphi;
  REAL vmid;
  REAL v;
  REAL dv;
  int k;

  // Constr data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  bus_counted = FUNC_get_bus_counted(f);
  if (!phi || !gphi || !bus_counted)
    return;

  // Bus data
  bus[0] = BRANCH_get_bus_from(br);
  bus[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(bus[k]);

  // Buses
  for (k = 0; k < 2; k++) {

    if (!bus_counted[bus_index[k]]) {
      
      // Generators
      for (gen = BUS_get_gen(bus[k]); gen != NULL; gen = GEN_get_next(gen)) {
	
	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var
	  
	  // Normalization factor
	  dv = GEN_get_Q_max(gen)-GEN_get_Q_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  // Values
	  vmid = (GEN_get_Q_max(gen)+GEN_get_Q_min(gen))/2.; // p.u.
	  v = VEC_get(var_values,GEN_get_index_Q(gen));

	  // phi
	  (*phi) += 0.5*pow((v-vmid)/dv,2.);

	  // gphi
	  gphi[GEN_get_index_Q(gen)] = (v-vmid)/(dv*dv);
	}

	if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) { // P var
	  	
	  // Normalization factor
	  dv = GEN_get_P_max(gen)-GEN_get_P_min(gen); // p.u.
	  if (dv < FUNC_REG_PQ_PARAM)
	    dv = FUNC_REG_PQ_PARAM;

	  // Values
	  vmid = (GEN_get_P_max(gen)+GEN_get_P_min(gen))/2.; // p.u.
	  v = VEC_get(var_values,GEN_get_index_P(gen));

	  // phi
	  (*phi) += 0.5*pow((v-vmid)/dv,2.);

	  // gphi
	  gphi[GEN_get_index_P(gen)] = (v-vmid)/(dv*dv);
	}
      }
    }
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void FUNC_REG_PQ_free(Func* f) {
  // Nothing
}
