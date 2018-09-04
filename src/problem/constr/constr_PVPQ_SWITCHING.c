/** @file constr_PVPQ_SWITCHING.c
 *  @brief This file defines the data structure and routines associated with the constraint of type PVPQ_SWITCHING.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/constr_PVPQ_SWITCHING.h>

struct Constr_PVPQ_SWITCHING_Data {

  char* fix_flag; 
};

Constr* CONSTR_PVPQ_SWITCHING_new(Net* net) {
  Constr* c = CONSTR_new(net);
  CONSTR_set_func_count_step(c,&CONSTR_PVPQ_SWITCHING_count_step);
  CONSTR_set_func_allocate(c,&CONSTR_PVPQ_SWITCHING_allocate);
  CONSTR_set_func_analyze_step(c,&CONSTR_PVPQ_SWITCHING_analyze_step);
  CONSTR_set_func_eval_step(c,&CONSTR_PVPQ_SWITCHING_eval_step);
  CONSTR_set_func_store_sens_step(c,&CONSTR_PVPQ_SWITCHING_store_sens_step);
  CONSTR_set_func_free(c,&CONSTR_PVPQ_SWITCHING_free);
  CONSTR_set_name(c,"PVPQ switching");
  return c;
}

void CONSTR_PVPQ_SWITCHING_count_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Gen* gen;
  int* A_nnz;
  int* A_row;
  int num;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;

  // Regulated bus (includes slack)
  if (BUS_is_regulated_by_gen(bus)) {
    
    num  = 0;
    
    // v
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
      num += 1;
    
    // Q
    for (gen = BUS_get_reg_gen(bus); gen != NULL; gen = GEN_get_reg_next(gen)) {
      if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage(gen))
        num += 1;
    }
    
    if (num > 1) {
      (*A_nnz) += num*(num-1);
      (*A_row) += num-1;
    }
  }
}

void CONSTR_PVPQ_SWITCHING_allocate(Constr* c) {

  // Local variables
  Constr_PVPQ_SWITCHING_Data* data;
  int num_vars;
  Bus* bus;
  Net* net;
  int i;
  int t;

  net = CONSTR_get_network(c);
  num_vars = NET_get_num_vars(net);

  // Data (var-dependent)
  CONSTR_PVPQ_SWITCHING_free(c);
  data = (Constr_PVPQ_SWITCHING_Data*)malloc(sizeof(Constr_PVPQ_SWITCHING_Data));
  ARRAY_zalloc(data->fix_flag,char,num_vars);
  for (i = 0; i < NET_get_num_buses(net); i++) {
    bus = NET_get_bus(net,i);
    if (BUS_is_regulated_by_gen(bus)) { // includes slack bus
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
        for (t = 0; t < NET_get_num_periods(net); t++)
          data->fix_flag[BUS_get_index_v_mag(bus,t)] = TRUE; // reg v starts fixed
      }
    }
  }
  CONSTR_set_data(c,(void*)data);
}

void CONSTR_PVPQ_SWITCHING_analyze_step(Constr* c, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Gen* gen1;
  Gen* gen2;
  Gen* gen3;
  int* A_nnz;
  int* A_row;
  Vec* b;
  Mat* A;
  REAL alpha1;
  REAL alpha2;
  REAL Q;
  REAL Q_min;
  REAL Q_max;
  BOOL Q_available;
  Constr_PVPQ_SWITCHING_Data* data;

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  data = (Constr_PVPQ_SWITCHING_Data*)CONSTR_get_data(c);

  // Check pointer
  if (!A_nnz || !A_row || !data || !bus)
    return;

  // Regulated bus (includes slack)
  if (BUS_is_regulated_by_gen(bus)) {

    // Q available
    Q_available = FALSE;
    for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
      if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage(gen1)) {
        Q_available = TRUE;
        break;
      }
    }
	  
    // v var and fixed
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) && Q_available &&
        data->fix_flag[BUS_get_index_v_mag(bus,t)]) {
      
      VEC_set(b,*A_row,BUS_get_v_set(bus,t));
      
      // v
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
      
      // Q
      for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
        if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage(gen1)) {
          MAT_set_i(A,*A_nnz,*A_row);
          MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen1,t));
          MAT_set_d(A,*A_nnz,0.);
          (*A_nnz)++;
        }
      }
      
      (*A_row)++;
    }
    
    // Q var and fixed
    for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
      if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) &&
          !GEN_is_on_outage(gen1) &&
          data->fix_flag[GEN_get_index_Q(gen1,t)]) {
	
        Q = GEN_get_Q(gen1,t);
        Q_max = GEN_get_Q_max(gen1);
        Q_min = GEN_get_Q_min(gen1);
	
        if (fabs(Q-Q_min) < fabs(Q-Q_max))
          VEC_set(b,*A_row,Q_min);
        else
          VEC_set(b,*A_row,Q_max);
	
        // v
        if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
          MAT_set_i(A,*A_nnz,*A_row);
          MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
          MAT_set_d(A,*A_nnz,0.);
          (*A_nnz)++;
        }
	
        // Q
        for (gen2 = BUS_get_reg_gen(bus); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
          if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage(gen2)) {
            MAT_set_i(A,*A_nnz,*A_row);
            MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen2,t));
            if (gen2 == gen1)
              MAT_set_d(A,*A_nnz,1.);
            else
              MAT_set_d(A,*A_nnz,0.);
            (*A_nnz)++;
          }
        }
	
        (*A_row)++;
      }
    }
    
    // Q var and free pairs
    gen1 = BUS_get_reg_gen(bus);
    while(gen1) {
      
      // Candidate 1
      if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) &&
          !GEN_is_on_outage(gen1) &&
          !data->fix_flag[GEN_get_index_Q(gen1,t)]) {
	
        for (gen2 = GEN_get_reg_next(gen1); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
	  
          // Candidate 2
          if (GEN_has_flags(gen2,FLAG_VARS,GEN_VAR_Q) &&
              !GEN_is_on_outage(gen2) &&
              !data->fix_flag[GEN_get_index_Q(gen2,t)]) {
	    
            VEC_set(b,*A_row,0.);
	    
            alpha1 = GEN_get_Q_par(gen1);
            if (alpha1 < CONSTR_PVPQ_SWITCHING_PARAM)
              alpha1 = CONSTR_PVPQ_SWITCHING_PARAM;
	    
            alpha2 = GEN_get_Q_par(gen2);
            if (alpha2 < CONSTR_PVPQ_SWITCHING_PARAM)
              alpha2 = CONSTR_PVPQ_SWITCHING_PARAM;
	    
            // v
            if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
              MAT_set_i(A,*A_nnz,*A_row);
              MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
              MAT_set_d(A,*A_nnz,0.);
              (*A_nnz)++;
            }
	    
            // Q
            for (gen3 = BUS_get_reg_gen(bus); gen3 != NULL; gen3 = GEN_get_reg_next(gen3)) {
              if (GEN_has_flags(gen3,FLAG_VARS,GEN_VAR_Q) && !GEN_is_on_outage(gen3)) {
                MAT_set_i(A,*A_nnz,*A_row);
                MAT_set_j(A,*A_nnz,GEN_get_index_Q(gen3,t));
                if (gen3 == gen1)
                  MAT_set_d(A,*A_nnz,alpha2);
                else if (gen3 == gen2)
                  MAT_set_d(A,*A_nnz,-alpha1);
                else
                  MAT_set_d(A,*A_nnz,0.);
                (*A_nnz)++;
              }
            }
	    
            (*A_row)++;
            break;
          }
        }
	
        // Move forward
        gen1 = gen2;
      }
      else {
	
        // Move forward
        gen1 = GEN_get_reg_next(gen1);
      }
    }
  }
}

void CONSTR_PVPQ_SWITCHING_eval_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* values, Vec* values_extra) {
  // Nothing to do
}

void CONSTR_PVPQ_SWITCHING_store_sens_step(Constr* c, Bus* bus, BusDC* busdc, int t, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_PVPQ_SWITCHING_free(Constr* c) {

  // Local variables
  Constr_PVPQ_SWITCHING_Data* data = (Constr_PVPQ_SWITCHING_Data*)CONSTR_get_data(c);

  // Free
  if (data) {
    if (data->fix_flag)
      free(data->fix_flag);
    free(data);
  }

  // Clear
  CONSTR_set_data(c,NULL);
}

char* CONSTR_PVPQ_SWITCHING_get_flags(Constr* c) {

  // Local variables
  Constr_PVPQ_SWITCHING_Data* data = (Constr_PVPQ_SWITCHING_Data*)CONSTR_get_data(c);

  // Check
  if (!data)
    return NULL;

  // Return
  return data->fix_flag;
}
