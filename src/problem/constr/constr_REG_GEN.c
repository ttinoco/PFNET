/** @file constr_REG_GEN.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_GEN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_GEN.h>

void CONSTR_REG_GEN_init(Constr* c) {
  
  // Local variables
  Net* net;
  int num_Jconstr;
  
  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 2*(NET_get_num_buses_reg_by_gen(net)-NET_get_num_slack_buses(net));
  CONSTR_set_Hcounter(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_GEN_clear(Constr* c) {

  // f
  VEC_set_zero(CONSTR_get_f(c));

  // J
  MAT_set_zero_d(CONSTR_get_J(c));

  // H
  MAT_array_set_zero_d(CONSTR_get_H_array(c),CONSTR_get_H_array_size(c));
  
  // Counters
  CONSTR_set_Acounter(c,0);
  CONSTR_set_Jcounter(c,0);
  CONSTR_set_Aconstr_index(c,0);
  CONSTR_set_Jconstr_index(c,0);
  CONSTR_clear_Hcounter(c);
  CONSTR_set_branch_counter(c,0);
   
  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_REG_GEN_count_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  int* Acounter;
  int* Jcounter;
  int* Aconstr_index;
  int* Jconstr_index;
  int* Hcounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter || !bus_counted)
    return;
 
  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      if (BUS_is_regulated_by_gen(bus) && // reg gen
	  !BUS_is_slack(bus)) {           // not slack
		
	// Linear
	//*******
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && // v var
	    BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var
	
	  // A
	  (*Acounter)++; // v
	  (*Acounter)++; // y
	  (*Acounter)++; // z

	  // Inc A cosntr index
	  (*Aconstr_index)++;
	}
	
	// Nonlinear
	//**********
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var

	  // J	  
	  (*Jcounter)++; // dCompY/dy
	  (*Jcounter)++; // dCompZ/dz

	  // H
	  Hcounter[*Jconstr_index]++;     // y and y
	  Hcounter[*Jconstr_index+1]++;   // z and z
	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {
	      Hcounter[*Jconstr_index]++;   // y and Q
	      Hcounter[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}
	
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var
	  
	    // J
	    (*Jcounter)++; // dcompY/dQ
	    (*Jcounter)++; // dcompZ/dQ

	    // H
	    Hcounter[*Jconstr_index]++;   // Q and Q
	    Hcounter[*Jconstr_index+1]++; // Q and Q
	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {
		Hcounter[*Jconstr_index]++;   // Q and Q1
		Hcounter[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
      }
    }	  
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_allocate(Constr *c) {

  // Local variables
  int Acounter;
  int Jcounter;
  int Aconstr_index;
  int Jconstr_index; 
  int* Hcounter;
  Mat* H_array;
  Mat* H;
  int H_comb_nnz;
  int num_vars;
  int i;
  
  Acounter = CONSTR_get_Acounter(c);
  Jcounter = CONSTR_get_Jcounter(c);
  Aconstr_index = CONSTR_get_Aconstr_index(c);
  Jconstr_index = CONSTR_get_Jconstr_index(c);
  Hcounter = CONSTR_get_Hcounter(c);
  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  
  // b
  CONSTR_set_b(c,VEC_new(Aconstr_index));

  // A
  CONSTR_set_A(c,MAT_new(Aconstr_index, // size1 (rows)
			 num_vars,      // size2 (cols)
			 Acounter));    // nnz
  
  // f
  CONSTR_set_f(c,VEC_new(Jconstr_index));

  // J
  CONSTR_set_J(c,MAT_new(Jconstr_index, // size1 (rows)
			 num_vars,      // size2 (cols)
			 Jcounter));    // nnz
  
  // H
  H_comb_nnz = 0;
  H_array = MAT_array_new(Jconstr_index);
  CONSTR_set_H_array(c,H_array,Jconstr_index);
  for (i = 0; i < Jconstr_index; i++) {
    H = MAT_array_get(H_array,i);
    MAT_set_nnz(H,Hcounter[i]);
    MAT_set_size1(H,num_vars);
    MAT_set_size2(H,num_vars);
    MAT_set_row_array(H,(int*)calloc(Hcounter[i],sizeof(int)));
    MAT_set_col_array(H,(int*)calloc(Hcounter[i],sizeof(int)));
    MAT_set_data_array(H,(REAL*)calloc(Hcounter[i],sizeof(REAL)));
    H_comb_nnz += Hcounter[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
}

void CONSTR_REG_GEN_analyze_branch(Constr* c, Branch* br) {

  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hy;
  Mat* Hz;
  int* Hi;
  int* Hj;
  int* Hi_comb;
  int* Hj_comb;
  int* Acounter;
  int* Jcounter;
  int* Aconstr_index;
  int* Jconstr_index;
  int* Hcounter;
  int Hcounter_comb;
  char* bus_counted;
  int bus_index[2];
  int k;
  int m;
  int temp;

  // Constr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  J = CONSTR_get_J(c);
  H_array = CONSTR_get_H_array(c);
  Acounter = CONSTR_get_Acounter_ptr(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Branch
  //*******

  CONSTR_inc_branch_counter(c);
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      if (BUS_is_regulated_by_gen(bus) &&  // reg by gen
	  !BUS_is_slack(bus)) {            // not slack
	
	// Hessians
	Hy = MAT_array_get(H_array,*Jconstr_index);
	Hz = MAT_array_get(H_array,*Jconstr_index+1);
	
	// Linear
	//*******
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && // v var
	    BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // yz var
	  
	  // b
	  VEC_set(b,*Aconstr_index,BUS_get_v_set(bus));
	
	  // A
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,BUS_get_index_v_mag(bus));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; // v

	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,BUS_get_index_y(bus));
	  MAT_set_d(A,*Acounter,-1.);
	  (*Acounter)++; // y

	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,BUS_get_index_z(bus));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; // z

	  // Inc A constr index
	  (*Aconstr_index)++;
	}

	// Nonlinear
	//**********
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	  
	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index);
	  MAT_set_j(J,*Jcounter,BUS_get_index_y(bus));
	  (*Jcounter)++; // dCompY/dy

	  MAT_set_i(J,*Jcounter,*Jconstr_index+1);
	  MAT_set_j(J,*Jcounter,BUS_get_index_z(bus));
	  (*Jcounter)++; // dCompZ/dz

	  // H
	  MAT_set_i(Hy,Hcounter[*Jconstr_index],BUS_get_index_y(bus));
	  MAT_set_j(Hy,Hcounter[*Jconstr_index],BUS_get_index_y(bus));
	  Hcounter[*Jconstr_index]++;     // y and y
	  
	  MAT_set_i(Hz,Hcounter[*Jconstr_index+1],BUS_get_index_z(bus));
	  MAT_set_j(Hz,Hcounter[*Jconstr_index+1],BUS_get_index_z(bus));
	  Hcounter[*Jconstr_index+1]++;   // z and z
	  
	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {

	      MAT_set_i(Hy,Hcounter[*Jconstr_index],BUS_get_index_y(bus));
	      MAT_set_j(Hy,Hcounter[*Jconstr_index],GEN_get_index_Q(rg));
	      Hcounter[*Jconstr_index]++;   // y and Q

	      MAT_set_i(Hz,Hcounter[*Jconstr_index+1],BUS_get_index_z(bus));
	      MAT_set_j(Hz,Hcounter[*Jconstr_index+1],GEN_get_index_Q(rg));
	      Hcounter[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}
	  
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var
	  
	    // J
	    MAT_set_i(J,*Jcounter,*Jconstr_index);
	    MAT_set_j(J,*Jcounter,GEN_get_index_Q(rg));
	    (*Jcounter)++; // dcompY/dQ
	    
	    MAT_set_i(J,*Jcounter,*Jconstr_index+1);
	    MAT_set_j(J,*Jcounter,GEN_get_index_Q(rg));
	    (*Jcounter)++; // dcompZ/dQ

	    // H
	    MAT_set_i(Hy,Hcounter[*Jconstr_index],GEN_get_index_Q(rg));
	    MAT_set_j(Hy,Hcounter[*Jconstr_index],GEN_get_index_Q(rg));
	    Hcounter[*Jconstr_index]++;   // Q and Q
	    
	    MAT_set_i(Hz,Hcounter[*Jconstr_index+1],GEN_get_index_Q(rg));
	    MAT_set_j(Hz,Hcounter[*Jconstr_index+1],GEN_get_index_Q(rg));
	    Hcounter[*Jconstr_index+1]++; // Q and Q

	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {
		
		MAT_set_i(Hy,Hcounter[*Jconstr_index],GEN_get_index_Q(rg));
		MAT_set_j(Hy,Hcounter[*Jconstr_index],GEN_get_index_Q(rg1));
		Hcounter[*Jconstr_index]++;   // Q and Q1
		
		MAT_set_i(Hz,Hcounter[*Jconstr_index+1],GEN_get_index_Q(rg));
		MAT_set_j(Hz,Hcounter[*Jconstr_index+1],GEN_get_index_Q(rg1));
		Hcounter[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
      }
    }	  
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }

  // Done
  if (CONSTR_get_branch_counter(c) == NET_get_num_branches(CONSTR_get_network(c))) {
    
    // Ensure lower triangular and save struct of H comb
    Hcounter_comb = 0;
    Hi_comb = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hj_comb = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (k = 0; k < CONSTR_get_H_array_size(c); k++) {
      Hi = MAT_get_row_array(MAT_array_get(H_array,k));
      Hj = MAT_get_col_array(MAT_array_get(H_array,k));
      for (m = 0; m < MAT_get_nnz(MAT_array_get(H_array,k)); m++) {
	if (Hi[m] < Hj[m]) {
	  temp = Hi[m];
	  Hi[m] = Hj[m];
	  Hj[m] = temp;
	}
	Hi_comb[Hcounter_comb] = Hi[m];
	Hj_comb[Hcounter_comb] = Hj[m];
	Hcounter_comb++;
      }
    }
  }
}

void CONSTR_REG_GEN_eval_branch(Constr* c, Branch *br, Vec* var_values) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* rg;
  Gen* rg1;
  Mat* H_array;
  REAL* f;
  REAL* J;
  REAL* Hy;
  REAL* Hz;
  int* Jcounter;
  int* Jconstr_index;
  int* Hcounter;
  char* bus_counted;
  int bus_index[2];
  int k;
  REAL v;
  REAL v_set;
  REAL y;
  REAL z;
  REAL Qsum;
  REAL Qmin;
  REAL Qmax;
  REAL Qy;
  REAL Qz;
  REAL sqrt_termY;
  REAL sqrt_termZ;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!f || !J || !Jcounter || !Jconstr_index || 
      !Hcounter || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);

  // Branch
  //*******

  CONSTR_inc_branch_counter(c);
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      if (BUS_is_regulated_by_gen(bus) &&  // reg by gen
	  !BUS_is_slack(bus)) {            // not slack

	// Hessians
	Hy = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index));
	Hz = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+1));
	
	// Q value
	Qsum = 0;
	Qmax = 0;
	Qmin = 0;
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q))
	    Qsum += VEC_get(var_values,GEN_get_index_Q(rg)); // p.u.
	  else
	    Qsum += GEN_get_Q(rg);      // p.u.
	  Qmax += GEN_get_Q_max(rg); // p.u.
	  Qmin += GEN_get_Q_min(rg); // p.u.
	}
	Qy = (Qsum-Qmin);
	Qz = (Qmax-Qsum);

	// yz values
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	  y = VEC_get(var_values,BUS_get_index_y(bus)); // p.u.
	  z = VEC_get(var_values,BUS_get_index_z(bus));	// p.u.		
	}
	else {
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	    v = VEC_get(var_values,BUS_get_index_v_mag(bus)); // p.u.
	  else
	    v = BUS_get_v_mag(bus);   // p.u.
	  v_set = BUS_get_v_set(bus); // p.u.
	  if (v > v_set) {
	    y = v-v_set;
	    z = 0;
	  }
	  else {
	    y = 0;
	    z = v_set-v;
	  } 
	}

	// Terms
	sqrt_termY = sqrt( Qy*Qy + y*y + 2*CONSTR_REG_GEN_PARAM );
	sqrt_termZ = sqrt( Qz*Qz + z*z + 2*CONSTR_REG_GEN_PARAM );

	// f
	f[*Jconstr_index] = Qy + y - sqrt_termY;   // CompY
	f[*Jconstr_index+1] = Qz + z - sqrt_termZ; // CompZ
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	  
	  // J
	  J[*Jcounter] = 1. - y/sqrt_termY;
	  (*Jcounter)++; // dCompY/dy
	  
	  J[*Jcounter] = 1. - z/sqrt_termZ;
	  (*Jcounter)++; // dCompZ/dz

	  // H
	  Hy[Hcounter[*Jconstr_index]] = -(Qy*Qy+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	  Hcounter[*Jconstr_index]++;     // y and y

	  Hz[Hcounter[*Jconstr_index+1]] = -(Qz*Qz+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	  Hcounter[*Jconstr_index+1]++;   // z and z
	  for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	    if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) {
	      
	      Hy[Hcounter[*Jconstr_index]] = Qy*y/pow(sqrt_termY,3.);
	      Hcounter[*Jconstr_index]++;   // y and Q

	      Hz[Hcounter[*Jconstr_index+1]] = -Qz*z/pow(sqrt_termZ,3.);
	      Hcounter[*Jconstr_index+1]++; // z and Q
	    }
	  }
	}
	
	for (rg = BUS_get_reg_gen(bus); rg != NULL; rg = GEN_get_reg_next(rg)) {
	  if (GEN_has_flags(rg,FLAG_VARS,GEN_VAR_Q)) { // Qg var
	  
	    // J
	    J[*Jcounter] = 1. - Qy/sqrt_termY;
	    (*Jcounter)++; // dcompY/dQ

	    J[*Jcounter] = -1. + Qz/sqrt_termZ;
	    (*Jcounter)++; // dcompZ/dQ

	    // H
	    Hy[Hcounter[*Jconstr_index]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
	    Hcounter[*Jconstr_index]++;   // Q and Q
	    
	    Hz[Hcounter[*Jconstr_index+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
	    Hcounter[*Jconstr_index+1]++; // Q and Q
	    
	    for (rg1 = GEN_get_reg_next(rg); rg1 != NULL; rg1 = GEN_get_reg_next(rg1)) {
	      if (GEN_has_flags(rg1,FLAG_VARS,GEN_VAR_Q)) {
	      
		Hy[Hcounter[*Jconstr_index]] = -(y*y+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termY,3.);
		Hcounter[*Jconstr_index]++;   // Q and Q1
		
		Hz[Hcounter[*Jconstr_index+1]] = -(z*z+2*CONSTR_REG_GEN_PARAM)/pow(sqrt_termZ,3.);
		Hcounter[*Jconstr_index+1]++; // Q and Q1
	      }
	    }
	  }
	}

	// Inc J constr index
	(*Jconstr_index)++; // dCompY
	(*Jconstr_index)++; // dCompZ
      }
    } 
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_store_sens_branch(Constr* c, Branch* br, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  int* Jconstr_index;
  char* bus_counted;
  int bus_index[2];
  REAL lamCompY;
  REAL lamCompZ;
  int k;
  
  // Constr data
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Jconstr_index || !bus_counted)
    return;
 
  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index[k] = BUS_get_index(buses[k]);
  
  // Buses
  //******

  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index[k]]) { // not counted yet
      
      if (BUS_is_regulated_by_gen(bus) && // reg gen
	  !BUS_is_slack(bus)) {           // not slack
	
	lamCompY = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // dCompY
	lamCompZ = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // dCompZ
	
	if (fabs(lamCompY) > fabs(lamCompZ))
	  BUS_set_sens_v_reg_by_gen(bus,lamCompY);
	else
	  BUS_set_sens_v_reg_by_gen(bus,lamCompZ);
      }
    }	  
    
    // Update counted flag
    bus_counted[bus_index[k]] = TRUE;
  }
}

void CONSTR_REG_GEN_free(Constr* c) {
  // Nothing
}
