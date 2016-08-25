/** @file constr_REG_TRAN.c
 *  @brief This file defines the data structure and routines associated with the constraint of type REG_TRAN.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_REG_TRAN.h>

void CONSTR_REG_TRAN_init(Constr* c) {
  
  // Local variables
  Net* net;
  int num_Jconstr;
  
  // Init
  net = CONSTR_get_network(c);
  num_Jconstr = 4*NET_get_num_tap_changers_v(CONSTR_get_network(c))*NET_get_num_periods(net);;
  CONSTR_set_Hcounter(c,(int*)calloc(num_Jconstr,sizeof(int)),num_Jconstr);
  CONSTR_set_data(c,NULL);
}

void CONSTR_REG_TRAN_clear(Constr* c) {
  
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
}

void CONSTR_REG_TRAN_count_step(Constr* c, Branch* br, int tau) {
  
  // Local variables
  Bus* reg_bus;
  int* Acounter;
  int* Jcounter;
  int* Aconstr_index;
  int* Jconstr_index;
  int* Hcounter;
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);
  
  // Check pointers
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);

    // Linear
    //*******
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&     // t var
	BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var
        
      // A
      (*Acounter)++; // t
      (*Acounter)++; // y
      (*Acounter)++; // z
      
      (*Aconstr_index)++;    
    }
    
    // Nonlinear constraints 1 (vmax,vmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var
      
      // J
      (*Jcounter)++; // dcompVmin/dy
      (*Jcounter)++; // dcompVmax/dz
      
      // H
      Hcounter[*Jconstr_index]++;   // y and y (vmin)
      Hcounter[*Jconstr_index+1]++; // z and z (vmax)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	Hcounter[*Jconstr_index]++;   // y and v (vmin)
	Hcounter[*Jconstr_index+1]++; // z and v (vmax)
      }
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
      // J
      (*Jcounter)++; // dcompVmin/dv
      (*Jcounter)++; // dcompVmax/dv
      
      // H
      Hcounter[*Jconstr_index]++;   // v and v (vmin)
      Hcounter[*Jconstr_index+1]++; // v and v (vmax)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	Hcounter[*Jconstr_index]++;   // v and vl (vmin)
	Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
      }
    }
    
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      (*Jcounter)++; // dcompVmin/dvl
      (*Jcounter)++; // dcompVmax/dvh

      // H 
      Hcounter[*Jconstr_index]++;   // vl and vl (vmin)
      Hcounter[*Jconstr_index+1]++; // vh and vh (vmax)
    }

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var

      // J
      (*Jcounter)++; // dcompTmax/dt
      (*Jcounter)++; // dcompTmin/dt

      // H
      Hcounter[*Jconstr_index+2]++; // t and t (tmax)
      Hcounter[*Jconstr_index+3]++; // t and t (tmin)
      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	Hcounter[*Jconstr_index+2]++; // t and vl (tmax)
	Hcounter[*Jconstr_index+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      (*Jcounter)++; // dcompTmax/dvl
      (*Jcounter)++; // dcompTmin/dvh

      // H 
      Hcounter[*Jconstr_index+2]++; // vl and vl (tmax)
      Hcounter[*Jconstr_index+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*Jconstr_index)++; // compVmin
    (*Jconstr_index)++; // compVmax
    (*Jconstr_index)++; // compTmax
    (*Jconstr_index)++; // compTmin
  }
}

void CONSTR_REG_TRAN_allocate(Constr* c) {

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

void CONSTR_REG_TRAN_analyze_step(Constr* c, Branch* br, int tau) {
  
  // Local variables
  Bus* reg_bus;
  Vec* b;
  Mat* A;
  Mat* J;
  Mat* H_array;
  Mat* Hvmin;
  Mat* Hvmax;
  Mat* Htmin;
  Mat* Htmax;
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
  int k;
  int m;
  int temp;
  int index_yz_vmin;
  int index_yz_vmax;
  int index_vvio_tmax;
  int index_vvio_tmin;
  int index_v;
  int index_vl;
  int index_vh;
  int index_t;
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

  // Check pointers
  if (!Acounter || !Jcounter || !Aconstr_index ||
      !Jconstr_index || !Hcounter)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;
  
  if (BRANCH_is_tap_changer_v(br)) {
      
    reg_bus = BRANCH_get_reg_bus(br);
    
    // Hessians (NOTE ORDER!!!)
    Hvmin = MAT_array_get(H_array,*Jconstr_index);
    Hvmax = MAT_array_get(H_array,*Jconstr_index+1);
    Htmax = MAT_array_get(H_array,*Jconstr_index+2);
    Htmin = MAT_array_get(H_array,*Jconstr_index+3);

    // Indices
    if (BRANCH_has_pos_ratio_v_sens(br)) {
      index_yz_vmin = BRANCH_get_index_ratio_y(br,tau);
      index_yz_vmax = BRANCH_get_index_ratio_z(br,tau);
      index_vvio_tmax = BUS_get_index_vl(reg_bus,tau);
      index_vvio_tmin = BUS_get_index_vh(reg_bus,tau);
    }
    else {
      index_yz_vmin = BRANCH_get_index_ratio_z(br,tau);
      index_yz_vmax = BRANCH_get_index_ratio_y(br,tau);
      index_vvio_tmax = BUS_get_index_vh(reg_bus,tau);
      index_vvio_tmin = BUS_get_index_vl(reg_bus,tau);
    }
    index_v = BUS_get_index_v_mag(reg_bus,tau);
    index_vl = BUS_get_index_vl(reg_bus,tau);
    index_vh = BUS_get_index_vh(reg_bus,tau);
    index_t = BRANCH_get_index_ratio(br,tau);
    index_y = BRANCH_get_index_ratio_y(br,tau);
    index_z = BRANCH_get_index_ratio_z(br,tau);
    
    // Linear (t = t_0 + y - z)
    //*************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO) &&     // t var
	BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var

      // b
      VEC_set(b,*Aconstr_index,BRANCH_get_ratio(br,tau)); // current ratio value
        
      // A
      MAT_set_i(A,*Acounter,*Aconstr_index);
      MAT_set_j(A,*Acounter,index_t);
      MAT_set_d(A,*Acounter,1.);
      (*Acounter)++; // t

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
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var
      
      // J
      MAT_set_i(J,*Jcounter,*Jconstr_index);
      MAT_set_j(J,*Jcounter,index_yz_vmin);
      (*Jcounter)++; // dcompVmin/dy
      
      MAT_set_i(J,*Jcounter,*Jconstr_index+1);
      MAT_set_j(J,*Jcounter,index_yz_vmax);
      (*Jcounter)++; // dcompVmax/dz
      
      // H	
      MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_yz_vmin);
      MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_yz_vmin);
      Hcounter[*Jconstr_index]++;   // y and y (vmin)
      
      MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_yz_vmax);
      MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_yz_vmax);
      Hcounter[*Jconstr_index+1]++; // z and z (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
	MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_yz_vmin);
	MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_v);
	Hcounter[*Jconstr_index]++;   // y and v (vmin)

	MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_yz_vmax);
	MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	Hcounter[*Jconstr_index+1]++; // z and v (vmax)
      }

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_yz_vmin);
	MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_vl);
	Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	
	MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_yz_vmax);
	MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_vl);
	Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
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

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	MAT_set_i(Hvmin,Hcounter[*Jconstr_index],index_v);
	MAT_set_j(Hvmin,Hcounter[*Jconstr_index],index_vl);
	Hcounter[*Jconstr_index]++;   // v and vl (vmin)

	MAT_set_i(Hvmax,Hcounter[*Jconstr_index+1],index_v);
	MAT_set_j(Hvmax,Hcounter[*Jconstr_index+1],index_vh);
	Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
      }
    }
    
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

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

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
      
      // J
      MAT_set_i(J,*Jcounter,*Jconstr_index+2);
      MAT_set_j(J,*Jcounter,index_t);
      (*Jcounter)++; // dcompTmax/dt

      MAT_set_i(J,*Jcounter,*Jconstr_index+3);
      MAT_set_j(J,*Jcounter,index_t);
      (*Jcounter)++; // dcompTmin/dt

      // H
      MAT_set_i(Htmax,Hcounter[*Jconstr_index+2],index_t);
      MAT_set_j(Htmax,Hcounter[*Jconstr_index+2],index_t);
      Hcounter[*Jconstr_index+2]++; // t and t (tmax)

      MAT_set_i(Htmin,Hcounter[*Jconstr_index+3],index_t);
      MAT_set_j(Htmin,Hcounter[*Jconstr_index+3],index_t);
      Hcounter[*Jconstr_index+3]++; // t and t (tmin)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	MAT_set_i(Htmax,Hcounter[*Jconstr_index+2],index_t);
	MAT_set_j(Htmax,Hcounter[*Jconstr_index+2],index_vvio_tmax);
	Hcounter[*Jconstr_index+2]++; // t and vl (tmax)

	MAT_set_i(Htmin,Hcounter[*Jconstr_index+3],index_t);
	MAT_set_j(Htmin,Hcounter[*Jconstr_index+3],index_vvio_tmin);
	Hcounter[*Jconstr_index+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      MAT_set_i(J,*Jcounter,*Jconstr_index+2);
      MAT_set_j(J,*Jcounter,index_vvio_tmax);
      (*Jcounter)++; // dcompTmax/dvl

      MAT_set_i(J,*Jcounter,*Jconstr_index+3);
      MAT_set_j(J,*Jcounter,index_vvio_tmin);
      (*Jcounter)++; // dcompTmin/dvh

      // H 
      MAT_set_i(Htmax,Hcounter[*Jconstr_index+2],index_vvio_tmax);
      MAT_set_j(Htmax,Hcounter[*Jconstr_index+2],index_vvio_tmax);
      Hcounter[*Jconstr_index+2]++; // vl and vl (tmax)

      MAT_set_i(Htmin,Hcounter[*Jconstr_index+3],index_vvio_tmin);
      MAT_set_j(Htmin,Hcounter[*Jconstr_index+3],index_vvio_tmin);
      Hcounter[*Jconstr_index+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*Jconstr_index)++; // compVmin
    (*Jconstr_index)++; // compVmax
    (*Jconstr_index)++; // compTmax
    (*Jconstr_index)++; // compTmin
  }

  // Done
  if ((tau == T-1) && (BRANCH_get_index(br) == NET_get_num_branches(CONSTR_get_network(c))-1)) {
    
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

void CONSTR_REG_TRAN_eval_step(Constr* c, Branch* br, int tau, Vec* var_values) {
  
  // Local variables
  Bus* reg_bus;
  REAL* f;
  REAL* J;
  Mat* H_array;
  REAL* Hvmin;
  REAL* Hvmax;
  REAL* Htmin;
  REAL* Htmax;
  int* Jcounter;
  int* Jconstr_index;
  int* Hcounter;
  REAL v;
  REAL vl;
  REAL vh;
  REAL vmin;
  REAL vmax;
  REAL t;
  REAL tmax;
  REAL tmin;
  REAL y;
  REAL z;
  REAL yz_vmin;
  REAL yz_vmax;
  REAL vvio_tmax;
  REAL vvio_tmin;  
  REAL sqrtermVmin;
  REAL sqrtermVmax;
  REAL sqrtermTmax;
  REAL sqrtermTmin;
  REAL norm = CONSTR_REG_TRAN_NORM;

  // Constr data
  f = VEC_get_data(CONSTR_get_f(c));
  J = MAT_get_data_array(CONSTR_get_J(c));
  H_array = CONSTR_get_H_array(c);
  Jcounter = CONSTR_get_Jcounter_ptr(c);
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);
  Hcounter = CONSTR_get_Hcounter(c);

  // Check pointers
  if (!f || !J || !Jcounter || 
      !Jconstr_index || !Hcounter)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;  
  
  if (BRANCH_is_tap_changer_v(br)) {
    
    reg_bus = BRANCH_get_reg_bus(br);

    // v values
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG))
      v = VEC_get(var_values,BUS_get_index_v_mag(reg_bus,tau));
    else
      v = BUS_get_v_mag(reg_bus,tau);
    vmax = BUS_get_v_max(reg_bus);
    vmin = BUS_get_v_min(reg_bus);
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
      vl = VEC_get(var_values,BUS_get_index_vl(reg_bus,tau));
      vh = VEC_get(var_values,BUS_get_index_vh(reg_bus,tau));
    }
    else {
      vl = 0;
      vh = 0;
    }
    
    // t values
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
      t = VEC_get(var_values,BRANCH_get_index_ratio(br,tau));
    }
    else
      t = BRANCH_get_ratio(br,tau);
    tmax = BRANCH_get_ratio_max(br);
    tmin = BRANCH_get_ratio_min(br);
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) {
      y = VEC_get(var_values,BRANCH_get_index_ratio_y(br,tau));
      z = VEC_get(var_values,BRANCH_get_index_ratio_z(br,tau));
    }
    else {
      if (t > BRANCH_get_ratio(br,tau)) {
	y = t-BRANCH_get_ratio(br,tau);
	z = 0;
      }
      else {
	y = 0;
	z = BRANCH_get_ratio(br,tau)-t;
      }
    }

    // values that depend on sensitivity
    if (BRANCH_has_pos_ratio_v_sens(br)) {
      yz_vmin = y;
      yz_vmax = z;
      vvio_tmax = vl;
      vvio_tmin = vh;
    }
    else {
      yz_vmin = z;
      yz_vmax = y;
      vvio_tmax = vh;
      vvio_tmin = vl;
    }

    // Terms
    sqrtermVmin = sqrt( (v+vl-vmin)*(v+vl-vmin) + yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermVmax = sqrt( (vmax-v+vh)*(vmax-v+vh) + yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermTmax = sqrt( (tmax-t)*(tmax-t) + vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM );
    sqrtermTmin = sqrt( (t-tmin)*(t-tmin) + vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM );
    
    // Hessians (NOTE ORDER!!!)
    Hvmin = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index));
    Hvmax = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+1));
    Htmax = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+2));
    Htmin = MAT_get_data_array(MAT_array_get(H_array,*Jconstr_index+3));

    // f
    f[*Jconstr_index] = ((v+vl-vmin) + yz_vmin - sqrtermVmin)*norm;   // vmin
    f[*Jconstr_index+1] = ((vmax-v+vh) + yz_vmax - sqrtermVmax)*norm; // vmax
    f[*Jconstr_index+2] = ((tmax-t) + vvio_tmax - sqrtermTmax)*norm;  // tmax
    f[*Jconstr_index+3] = ((t-tmin) + vvio_tmin - sqrtermTmin)*norm;  // tmin
        
    // Nonlinear constraints 1 (vmin,vmax)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // yz var
      
      // J
      J[*Jcounter] = (1.-yz_vmin/sqrtermVmin)*norm;
      (*Jcounter)++; // dcompVmin/dy
      
      J[*Jcounter] = (1.-yz_vmax/sqrtermVmax)*norm;
      (*Jcounter)++; // dcompVmax/dz
      
      // H	
      Hvmin[Hcounter[*Jconstr_index]] = -(((v+vl-vmin)*(v+vl-vmin)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      Hcounter[*Jconstr_index]++;   // y and y (vmin)
      
      Hvmax[Hcounter[*Jconstr_index+1]] = -(((vmax-v+vh)*(vmax-v+vh)+2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      Hcounter[*Jconstr_index+1]++; // z and z (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) {
	
	Hvmin[Hcounter[*Jconstr_index]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
	Hcounter[*Jconstr_index]++;   // y and v (vmin)

	Hvmax[Hcounter[*Jconstr_index+1]] = -((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
	Hcounter[*Jconstr_index+1]++; // z and v (vmax)
      }

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	Hvmin[Hcounter[*Jconstr_index]] = ((v+vl-vmin)*yz_vmin/pow(sqrtermVmin,3.))*norm;
	Hcounter[*Jconstr_index]++;   // y and vl (vmin)
	
	Hvmax[Hcounter[*Jconstr_index+1]] = ((vmax-v+vh)*yz_vmax/pow(sqrtermVmax,3.))*norm;
	Hcounter[*Jconstr_index+1]++; // z and vh (vmax)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VMAG)) { // v var
      
      // J
      J[*Jcounter] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
      (*Jcounter)++; // dcompVmin/dv

      J[*Jcounter] = -((1.-(vmax-v+vh)/sqrtermVmax))*norm;
      (*Jcounter)++; // dcompVmax/dv
      
      // H
      Hvmin[Hcounter[*Jconstr_index]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      Hcounter[*Jconstr_index]++;   // v and v (vmin)
      
      Hvmax[Hcounter[*Jconstr_index+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      Hcounter[*Jconstr_index+1]++; // v and v (vmax)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {
	
	Hvmin[Hcounter[*Jconstr_index]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
	Hcounter[*Jconstr_index]++;   // v and vl (vmin)

	Hvmax[Hcounter[*Jconstr_index+1]] = ((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
	Hcounter[*Jconstr_index+1]++; // v and vh (vmax)
      }
    }
    
    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      J[*Jcounter] = (1.-(v+vl-vmin)/sqrtermVmin)*norm;
      (*Jcounter)++; // dcompVmin/dvl

      J[*Jcounter] = (1.-(vmax-v+vh)/sqrtermVmax)*norm;
      (*Jcounter)++; // dcompVmax/dvh

      // H 
      Hvmin[Hcounter[*Jconstr_index]] = -((yz_vmin*yz_vmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmin,3.))*norm;
      Hcounter[*Jconstr_index]++;   // vl and vl (vmin)

      Hvmax[Hcounter[*Jconstr_index+1]] = -((yz_vmax*yz_vmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermVmax,3.))*norm;
      Hcounter[*Jconstr_index+1]++; // vh and vh (vmax)
    }

    // Nonlinear constraints 2 (tmax,tmin)
    //************************************
    if (BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) { // t var
      
      // J
      J[*Jcounter] = -(1.-(tmax-t)/sqrtermTmax)*norm;
      (*Jcounter)++; // dcompTmax/dt
      
      J[*Jcounter] = (1.-(t-tmin)/sqrtermTmin)*norm;
      (*Jcounter)++; // dcompTmin/dt

      // H
      Htmax[Hcounter[*Jconstr_index+2]] = -((vvio_tmax*vvio_tmax + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
      Hcounter[*Jconstr_index+2]++; // t and t (tmax)

      Htmin[Hcounter[*Jconstr_index+3]] = -((vvio_tmin*vvio_tmin + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
      Hcounter[*Jconstr_index+3]++; // t and t (tmin)

      if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) {

	Htmax[Hcounter[*Jconstr_index+2]] = -(vvio_tmax*(tmax-t)/pow(sqrtermTmax,3.))*norm;
	Hcounter[*Jconstr_index+2]++; // t and vl (tmax)

	Htmin[Hcounter[*Jconstr_index+3]] = (vvio_tmin*(t-tmin)/pow(sqrtermTmin,3.))*norm;
	Hcounter[*Jconstr_index+3]++; // t and vh (tmin)
      }
    }

    if (BUS_has_flags(reg_bus,FLAG_VARS,BUS_VAR_VVIO)) { // vl and vh var

      // J
      J[*Jcounter] = (1.-vvio_tmax/sqrtermTmax)*norm;
      (*Jcounter)++; // dcompTmax/dvl

      J[*Jcounter] = (1.-vvio_tmin/sqrtermTmin)*norm;
      (*Jcounter)++; // dcompTmin/dvh

      // H 
      Htmax[Hcounter[*Jconstr_index+2]] = -(((tmax-t)*(tmax-t) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmax,3.))*norm;
      Hcounter[*Jconstr_index+2]++; // vl and vl (tmax)

      Htmin[Hcounter[*Jconstr_index+3]] = -(((t-tmin)*(t-tmin) + 2*CONSTR_REG_TRAN_PARAM)/pow(sqrtermTmin,3.))*norm;
      Hcounter[*Jconstr_index+3]++; // vh and vh (tmin)
    }

    // Inc J constr index
    (*Jconstr_index)++; // compVmin
    (*Jconstr_index)++; // compVmax
    (*Jconstr_index)++; // compTmax
    (*Jconstr_index)++; // compTmin
  }
}

void CONSTR_REG_TRAN_store_sens_step(Constr* c, Branch* br, int tau, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {

  // Local variables
  Bus* reg_bus;
  int* Jconstr_index;
  REAL lamCompVmin;
  REAL lamCompVmax;
  REAL lamCompTmax;
  REAL lamCompTmin;
  
  // Constr data
  Jconstr_index = CONSTR_get_Jconstr_index_ptr(c);

  // Check pointer
  if (!Jconstr_index)
    return;

  // Check outage
  if (BRANCH_is_on_outage(br))
    return;

  if (BRANCH_is_tap_changer_v(br)) {

    reg_bus = BRANCH_get_reg_bus(br);
    
    lamCompVmin = VEC_get(sf,*Jconstr_index);
    (*Jconstr_index)++; // compVmin
    lamCompVmax = VEC_get(sf,*Jconstr_index);
    (*Jconstr_index)++; // compVmax
    lamCompTmax = VEC_get(sf,*Jconstr_index);
    (*Jconstr_index)++; // compTmax
    lamCompTmin = VEC_get(sf,*Jconstr_index);
    (*Jconstr_index)++; // compTmin

    if (fabs(lamCompVmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmin,tau);
    if (fabs(lamCompVmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompVmax,tau);
    if (fabs(lamCompTmax) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmax,tau);
    if (fabs(lamCompTmin) > fabs(BUS_get_sens_v_reg_by_tran(reg_bus,tau)))
      BUS_set_sens_v_reg_by_tran(reg_bus,lamCompTmin,tau);
  }
}

void CONSTR_REG_TRAN_free(Constr* c) {
  // Nothing
}
