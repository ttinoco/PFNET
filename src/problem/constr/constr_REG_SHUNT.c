/** @file constr_REG_SHUNT.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_SHUNT.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_SHUNT.h>

void CONSTR_REG_SHUNT_init(Constr* c) {
  
  // Local variables
  Net* net;
  int num_Jconstr;
  
  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 4*NET_get_num_switched_shunts(CONSTR_get_network(c))*NET_get_num_periods(net);
  CONSTR_set_Hcounter(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_SHUNT_clear(Constr* c) {
  
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
   
  // Flags
  CONSTR_clear_bus_counted(c);  
}

void CONSTR_REG_SHUNT_count_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  int* Acounter;
  int* Jcounter;
  int* Aconstr_index;
  int* Jconstr_index;
  int* Hcounter;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  //******
  
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Shunts
      //*******
      
      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {
	
	// Linear
	//*******
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&     // b var
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var
	  
	  // A
	  (*Acounter)++; // b
	  (*Acounter)++; // y
	  (*Acounter)++; // z
	  
	  (*Aconstr_index)++;    
	}
    
	// Nonlinear constraints 1 (vmax,vmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var
      
	  // J
	  (*Jcounter)++; // dcompVmin/dy
	  (*Jcounter)++; // dcompVmax/dz
      
	  // H
	  Hcounter[*Jconstr_index]++;   // y and y (vmin)
	  Hcounter[*Jconstr_index+1]++; // z and z (vmax)
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	    Hcounter[*Jconstr_index]++;   // y and v (vmin)
	    Hcounter[*Jconstr_index+1]++; // z and v (vmax)
	  }
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	    Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	    Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
	  }
	}

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
	  // J
	  (*Jcounter)++; // dcompVmin/dv
	  (*Jcounter)++; // dcompVmax/dv
      
	  // H
	  Hcounter[*Jconstr_index]++;   // v and v (vmin)
	  Hcounter[*Jconstr_index+1]++; // v and v (vmax)
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	    Hcounter[*Jconstr_index]++;   // v and vl (vmin)
	    Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
	  }
	}
    
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var
	  
	  // J
	  (*Jcounter)++; // dcompVmin/dvl
	  (*Jcounter)++; // dcompVmax/dvh

	  // H 
	  Hcounter[*Jconstr_index]++;   // vl and vl (vmin)
	  Hcounter[*Jconstr_index+1]++; // vh and vh (vmax)
	}

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

	  // J
	  (*Jcounter)++; // dcompBmax/db
	  (*Jcounter)++; // dcompBmin/db

	  // H
	  Hcounter[*Jconstr_index+2]++; // b and b (bmax)
	  Hcounter[*Jconstr_index+3]++; // b and b (bmin)
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	    Hcounter[*Jconstr_index+2]++; // b and vl (bmax)
	    Hcounter[*Jconstr_index+3]++; // b and vh (bmin)
	  }
	}
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	  // J
	  (*Jcounter)++; // dcompBmax/dvl
	  (*Jcounter)++; // dcompBmin/dvh

	  // H 
	  Hcounter[*Jconstr_index+2]++; // vl and vl (bmax)
	  Hcounter[*Jconstr_index+3]++; // vh and vh (bmin)
	}

	// Inc J constr index
	(*Jconstr_index)++; // compVmin
	(*Jconstr_index)++; // compVmax
	(*Jconstr_index)++; // compBmax
	(*Jconstr_index)++; // compBmin
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_SHUNT_allocate(Constr* c) {

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
    MAT_set_data_array(H,(REAL*)malloc(Hcounter[i]*sizeof(REAL)));
    H_comb_nnz += Hcounter[i];
  }

  // H combined
  CONSTR_set_H_combined(c,MAT_new(num_vars,     // size1 (rows)
				  num_vars,     // size2 (cols)
				  H_comb_nnz)); // nnz
}

void CONSTR_REG_SHUNT_analyze_step(Constr* c, Branch* br, int t) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hvmin;
  Mat* Hvmax;
  Mat* Hbmin;
  Mat* Hbmax;
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
  int bus_index_t[2];
  int k;
  int m;
  int temp;
  int index_v;
  int index_vl;
  int index_vh;
  int index_b;
  int index_y;
  int index_z;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

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

  // Check pointers
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******
  
  // Buses
  //******
  
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet
 
      // Shunts
      //*******
      
      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {
    
	// Hessians (NOTE ORDER!!!)
	Hvmin = MAT_array_get(H_array,*Jconstr_index);
	Hvmax = MAT_array_get(H_array,*Jconstr_index+1);
	Hbmax = MAT_array_get(H_array,*Jconstr_index+2);
	Hbmin = MAT_array_get(H_array,*Jconstr_index+3);

	// Indices
	index_v = BUS_get_index_v_mag(bus,t);
	index_vl = BUS_get_index_vl(bus,t);
	index_vh = BUS_get_index_vh(bus,t);
	index_b = SHUNT_get_index_b(shunt,t);
	index_y = SHUNT_get_index_y(shunt,t);
	index_z = SHUNT_get_index_z(shunt,t);
    
	// Linear (b = b_0 + y - z)
	//*************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&     // b var
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var

	  // b
	  VEC_set(b,*Aconstr_index,SHUNT_get_b(shunt,t)); // current susceptance value
        
	  // A
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,index_b);
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; // b
	  
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,index_y);
	  MAT_set_d(A,*Acounter,-1.);
	  (*Acounter)++; // y
	  
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,index_z);
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++; // z
	  
	  (*Aconstr_index)++;    
	}
    
	// Nonlinear constraints 1 (vmin,vmax)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var
	  
	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index);
	  MAT_set_j(J,*Jcounter,index_y);
	  (*Jcounter)++; // dcompVmin/dy
	  
	  MAT_set_i(J,*Jcounter,*Jconstr_index+1);
	  MAT_set_j(J,*Jcounter,index_z);
	  (*Jcounter)++; // dcompVmax/dz
      
	  // H	
	  MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_y);
	  MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_y);
	  Hcounter[*Jconstr_index]++;   // y and y (vmin)
      
	  MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_z);
	  MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_z);
	  Hcounter[*Jconstr_index+1]++; // z and z (vmax)

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
	    MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_y);
	    MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_v);
	    Hcounter[*Jconstr_index]++;   // y and v (vmin)

	    MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_z);
	    MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	    Hcounter[*Jconstr_index+1]++; // z and v (vmax)
	  }

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	    MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_y);
	    MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_vl);
	    Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	
	    MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_z);
	    MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_vl);
	    Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
	  }
	}

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index);
	  MAT_set_j(J,*Jcounter,index_v);
	  (*Jcounter)++; // dcompVmin/dv

	  MAT_set_i(J,*Jcounter,*Jconstr_index+1);
	  MAT_set_j(J,*Jcounter,index_v);
	  (*Jcounter)++; // dcompVmax/dv
      
	  // H
	  MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_v);
	  MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_v);
	  Hcounter[*Jconstr_index]++;   // v and v (vmin)
      
	  MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	  MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	  Hcounter[*Jconstr_index+1]++; // v and v (vmax)

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	    
	    MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_v);
	    MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_vl);
	    Hcounter[*Jconstr_index]++;   // v and vl (vmin)
	    
	    MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	    MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_vh);
	    Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
	  }
	}
    
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index);
	  MAT_set_j(J,*Jcounter,index_vl);
	  (*Jcounter)++; // dcompVmin/dvl

	  MAT_set_i(J,*Jcounter,*Jconstr_index+1);
	  MAT_set_j(J,*Jcounter,index_vh);
	  (*Jcounter)++; // dcompVmax/dvh

	  // H 
	  MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_vl);
	  MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_vl);
	  Hcounter[*Jconstr_index]++;   // vl and vl (vmin)

	  MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_vh);
	  MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_vh);
	  Hcounter[*Jconstr_index+1]++; // vh and vh (vmax)
	}

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var
      
	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index+2);
	  MAT_set_j(J,*Jcounter,index_b);
	  (*Jcounter)++; // dcompBmax/db

	  MAT_set_i(J,*Jcounter,*Jconstr_index+3);
	  MAT_set_j(J,*Jcounter,index_b);
	  (*Jcounter)++; // dcompBmin/db
	  
	  // H
	  MAT_set_i(Hbmax,Hcounter[*Jconstr_index+2],index_b);
	  MAT_set_j(Hbmax,Hcounter[*Jconstr_index+2],index_b);
	  Hcounter[*Jconstr_index+2]++; // b and b (bmax)

	  MAT_set_i(Hbmin,Hcounter[*Jconstr_index+3],index_b);
	  MAT_set_j(Hbmin,Hcounter[*Jconstr_index+3],index_b);
	  Hcounter[*Jconstr_index+3]++; // b and b (bmin)

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	    MAT_set_i(Hbmax,Hcounter[*Jconstr_index+2],index_b);
	    MAT_set_j(Hbmax,Hcounter[*Jconstr_index+2],index_vl);
	    Hcounter[*Jconstr_index+2]++; // b and vl (bmax)

	    MAT_set_i(Hbmin,Hcounter[*Jconstr_index+3],index_b);
	    MAT_set_j(Hbmin,Hcounter[*Jconstr_index+3],index_vh);
	    Hcounter[*Jconstr_index+3]++; // b and vh (bmin)
	  }
	}

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	  // J
	  MAT_set_i(J,*Jcounter,*Jconstr_index+2);
	  MAT_set_j(J,*Jcounter,index_vl);
	  (*Jcounter)++; // dcompBmax/dvl

	  MAT_set_i(J,*Jcounter,*Jconstr_index+3);
	  MAT_set_j(J,*Jcounter,index_vh);
	  (*Jcounter)++; // dcompBmin/dvh

	  // H 
	  MAT_set_i(Hbmax,Hcounter[*Jconstr_index+2],index_vl);
	  MAT_set_j(Hbmax,Hcounter[*Jconstr_index+2],index_vl);
	  Hcounter[*Jconstr_index+2]++; // vl and vl (bmax)
	  
	  MAT_set_i(Hbmin,Hcounter[*Jconstr_index+3],index_vh);
	  MAT_set_j(Hbmin,Hcounter[*Jconstr_index+3],index_vh);
	  Hcounter[*Jconstr_index+3]++; // vh and vh (bmin)
	}

	// Inc J constr index
	(*Jconstr_index)++; // compVmin
	(*Jconstr_index)++; // compVmax
	(*Jconstr_index)++; // compBmax
	(*Jconstr_index)++; // compBmin
      }
    }
    
    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }

  // Done
  if ((t == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {
    
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

void CONSTR_REG_SHUNT_eval_step(Constr* c, Branch* br, int t, Vec* var_values) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* Hvmin;
  REAL* Hvmax;
  REAL* Hbmin;
  REAL* Hbmax;
  int* Jcounter;
  int* Jconstr_index;
  int* Hcounter;
  char* bus_counted;
  int bus_index_t[2];
  int k;
  REAL v;
  REAL vl;
  REAL vh;
  REAL vmin;
  REAL vmax;
  REAL b;
  REAL bmax;
  REAL bmin;
  REAL y;
  REAL z;  
  REAL sqrtermVmin;
  REAL sqrtermVmax;
  REAL sqrtermBmax;
  REAL sqrtermBmin;
  REAL norm = CONSTR_REG_SHUNT_NORM;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!f || !J || !Jcounter || !Jconstr_index || 
      !Hcounter || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Branch
  //*******

  // Buses
  //******
  
  for (k = 0; k < 2; k++) {
    
    bus = buses[k];
    
    if (!bus_counted[bus_index_t[k]]) { // not counted yet
      
      // Shunts
      //*******
      
      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {
    
	// v values
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
	  v = VEC_get(var_values,BUS_get_index_v_mag(bus,t));
	else
	  v = BUS_get_v_mag(bus,t);
	vmax = BUS_get_v_max(bus);
	vmin = BUS_get_v_min(bus);
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	  vl = VEC_get(var_values,BUS_get_index_vl(bus,t));
	  vh = VEC_get(var_values,BUS_get_index_vh(bus,t));
	}
	else {
	  vl = 0;
	  vh = 0;
	}
    
	// b values
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  b = VEC_get(var_values,SHUNT_get_index_b(shunt,t));
	}
	else
	  b = SHUNT_get_b(shunt,t);
	bmax = SHUNT_get_b_max(shunt);
	bmin = SHUNT_get_b_min(shunt);
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) {
	  y = VEC_get(var_values,SHUNT_get_index_y(shunt,t));
	  z = VEC_get(var_values,SHUNT_get_index_z(shunt,t));
	}
	else {
	  if (b > SHUNT_get_b(shunt,t)) {
	    y = b-SHUNT_get_b(shunt,t);
	    z = 0;
	  }
	  else {
	    y = 0;
	    z = SHUNT_get_b(shunt,t)-b;
	  }
	}

	// Terms
	sqrtermVmin = sqrt( (v+vl-vmin)*(v+vl-vmin) + y*y + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermVmax = sqrt( (vmax-v+vh)*(vmax-v+vh) + z*z + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermBmax = sqrt( (bmax-b)*(bmax-b) + vl*vl + 2*CONSTR_REG_SHUNT_PARAM );
	sqrtermBmin = sqrt( (b-bmin)*(b-bmin) + vh*vh + 2*CONSTR_REG_SHUNT_PARAM );
    
	// Hessians (NOTE ORDER!!!)
	Hvmin = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index));
	Hvmax = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+1));
	Hbmax = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+2));
	Hbmin = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+3));

	// f
	f[*Jconstr_index] = ((v+vl-vmin) + y - sqrtermVmin)*norm;   // vmin
	f[*Jconstr_index+1] = ((vmax-v+vh) + z - sqrtermVmax)*norm; // vmax
	f[*Jconstr_index+2] = ((bmax-b) + vl - sqrtermBmax)*norm;   // bmax
	f[*Jconstr_index+3] = ((b-bmin) + vh - sqrtermBmin)*norm;   // bmin
        
	// Nonlinear constraints 1 (vmin,vmax)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // yz var
      
	  // J
	  J[*Jcounter] = (1.-y/sqrtermVmin)*norm;
	  (*Jcounter)++; // dcompVmin/dy
      
	  J[*Jcounter] = (1.-z/sqrtermVmax)*norm;
	  (*Jcounter)++; // dcompVmax/dz
      
	  // H	
	  Hvmin[Hcounter[*Jconstr_index]] = -(((v+vl-vmin)*(v+vl-vmin)+2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	  Hcounter[*Jconstr_index]++;   // y and y (vmin)
      
	  Hvmax[Hcounter[*Jconstr_index+1]] = -(((vmax-v+vh)*(vmax-v+vh)+2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	  Hcounter[*Jconstr_index+1]++; // z and z (vmax)

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
	    Hvmin[Hcounter[*Jconstr_index]] = ((v+vl-vmin)*y/pow(sqrtermVmin,3.))*norm;
	    Hcounter[*Jconstr_index]++;   // y and v (vmin)

	    Hvmax[Hcounter[*Jconstr_index+1]] = -((vmax-v+vh)*z/pow(sqrtermVmax,3.))*norm;
	    Hcounter[*Jconstr_index+1]++; // z and v (vmax)
	  }

	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	    Hvmin[Hcounter[*Jconstr_index]] = ((v+vl-vmin)*y/pow(sqrtermVmin,3.))*norm;
	    Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	
	    Hvmax[Hcounter[*Jconstr_index+1]] = ((vmax-v+vh)*z/pow(sqrtermVmax,3.))*norm;
	    Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
	  }
	}

	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
	  // J
	  J[*Jcounter] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
	  (*Jcounter)++; // dcompVmin/dv

	  J[*Jcounter] = -((1.-(vmax-v+vh)/sqrtermVmax))*norm;
	  (*Jcounter)++; // dcompVmax/dv
	  
	  // H
	  Hvmin[Hcounter[*Jconstr_index]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	  Hcounter[*Jconstr_index]++;   // v and v (vmin)
	  
	  Hvmax[Hcounter[*Jconstr_index+1]] = -((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	  Hcounter[*Jconstr_index+1]++; // v and v (vmax)
	  
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	    Hvmin[Hcounter[*Jconstr_index]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	    Hcounter[*Jconstr_index]++;   // v and vl (vmin)

	    Hvmax[Hcounter[*Jconstr_index+1]] = ((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	    Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
	  }
	}
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	  // J
	  J[*Jcounter] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
	  (*Jcounter)++; // dcompVmin/dvl

	  J[*Jcounter] = (1.-(vmax-v+vh)/sqrtermVmax)*norm;
	  (*Jcounter)++; // dcompVmax/dvh

	  // H 
	  Hvmin[Hcounter[*Jconstr_index]] = -((y*y + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmin,3.))*norm;
	  Hcounter[*Jconstr_index]++;   // vl and vl (vmin)
	  
	  Hvmax[Hcounter[*Jconstr_index+1]] = -((z*z + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermVmax,3.))*norm;
	  Hcounter[*Jconstr_index+1]++; // vh and vh (vmax)
	}

	// Nonlinear constraints 2 (bmax,bmin)
	//************************************
	if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // t var
      
	  // J
	  J[*Jcounter] = -(1.-(bmax-b)/sqrtermBmax)*norm;
	  (*Jcounter)++; // dcompBmax/db
	  
	  J[*Jcounter] = (1.-(b-bmin)/sqrtermBmin)*norm;
	  (*Jcounter)++; // dcompBmin/db
	  
	  // H
	  Hbmax[Hcounter[*Jconstr_index+2]] = -((vl*vl + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmax,3.))*norm;
	  Hcounter[*Jconstr_index+2]++; // b and b (bmax)
	  
	  Hbmin[Hcounter[*Jconstr_index+3]] = -((vh*vh + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmin,3.))*norm;
	  Hcounter[*Jconstr_index+3]++; // b and b (bmin)
	  
	  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {

	    Hbmax[Hcounter[*Jconstr_index+2]] = -(vl*(bmax-b)/pow(sqrtermBmax,3.))*norm;
	    Hcounter[*Jconstr_index+2]++; // b and vl (bmax)

	    Hbmin[Hcounter[*Jconstr_index+3]] = (vh*(b-bmin)/pow(sqrtermBmin,3.))*norm;
	    Hcounter[*Jconstr_index+3]++; // b and vh (bmin)
	  }
	}
	
	if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

	  // J
	  J[*Jcounter] = (1.-vl/sqrtermBmax)*norm;
	  (*Jcounter)++; // dcompBmax/dvl

	  J[*Jcounter] = (1.-vh/sqrtermBmin)*norm;
	  (*Jcounter)++; // dcompBmin/dvh

	  // H 
	  Hbmax[Hcounter[*Jconstr_index+2]] = -(((bmax-b)*(bmax-b) + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmax,3.))*norm;
	  Hcounter[*Jconstr_index+2]++; // vl and vl (bmax)

	  Hbmin[Hcounter[*Jconstr_index+3]] = -(((b-bmin)*(b-bmin) + 2*CONSTR_REG_SHUNT_PARAM)/pow(sqrtermBmin,3.))*norm;
	  Hcounter[*Jconstr_index+3]++; // vh and vh (bmin)
	}

	// Inc J constr index
	(*Jconstr_index)++; // compVmin
	(*Jconstr_index)++; // compVmax
	(*Jconstr_index)++; // compBmax
	(*Jconstr_index)++; // compBmin
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }
}

void CONSTR_REG_SHUNT_store_sens_step(Constr* c, Branch* br, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Shunt* shunt;
  int* Jconstr_index;
  char* bus_counted;
  int bus_index_t[2];
  REAL lamCompVmin;
  REAL lamCompVmax;
  REAL lamCompBmax;
  REAL lamCompBmin;
  int k;
  int T;

  // Number of periods
  T = BRANCH_get_num_periods(br);
  
  // Constr data
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);

  // Check pointers
  if (!Jconstr_index || !bus_counted)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);
  for (k = 0; k < 2; k++)
    bus_index_t[k] = BUS_get_index(buses[k])*T+t;

  // Buses
  //******
  
  for (k = 0; k < 2; k++) {

    bus = buses[k];

    if (!bus_counted[bus_index_t[k]]) { // not counted yet

      // Shunts
      //*******
      
      for (shunt = BUS_get_reg_shunt(bus); shunt != NULL; shunt = SHUNT_get_reg_next(shunt)) {

	lamCompVmin = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // compVmin
	lamCompVmax = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // compVmax
	lamCompBmax = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // compBmax
	lamCompBmin = VEC_get(sf,*Jconstr_index);
	(*Jconstr_index)++; // compBmin

	if (fabs(lamCompVmin) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompVmin,t);
	if (fabs(lamCompVmax) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompVmax,t);
	if (fabs(lamCompBmax) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompBmax,t);
	if (fabs(lamCompBmin) > fabs(BUS_get_sens_v_reg_by_shunt(bus,t)))
	  BUS_set_sens_v_reg_by_shunt(bus,lamCompBmin,t);
      }
    }

    // Update counted flag
    bus_counted[bus_index_t[k]] = TRUE;
  }  
}

void CONSTR_REG_SHUNT_free(Constr* c) {
  // Nothing
}
