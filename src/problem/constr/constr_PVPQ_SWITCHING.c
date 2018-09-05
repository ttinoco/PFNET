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
#include <pfnet/reg_obj.h>
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
  int* A_nnz;
  int* A_row;
  int num;

  // Constr data
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);

  // Check pointer
  if (!A_nnz || !A_row || !bus)
    return;
      
  // Regulated bus
  if (BUS_is_v_set_regulated(bus)) {
    
    num = 0;
    
    // v
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG))
      num += 1;
    
    // Q
    num += REG_OBJ_count_candidates(bus);
    
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
    if (BUS_is_v_set_regulated(bus)) {
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
  int* A_nnz;
  int* A_row;
  Vec* b;
  Mat* A;
  REAL alpha1;
  REAL alpha2;
  REAL Q;
  REAL Q_min;
  REAL Q_max;
  Constr_PVPQ_SWITCHING_Data* data;

  void* obj1;
  char obj1_type;
  void* obj2;
  char obj2_type;
  void* obj3;
  char obj3_type;

  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  A_nnz = CONSTR_get_A_nnz_ptr(c);
  A_row = CONSTR_get_A_row_ptr(c);
  data = (Constr_PVPQ_SWITCHING_Data*)CONSTR_get_data(c);

  // Check pointer
  if (!A_nnz || !A_row || !data || !bus)
    return;
  
  // Regulated bus
  if (BUS_is_v_set_regulated(bus)) {
	
    // v var and fixed
    //****************
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
        REG_OBJ_count_candidates(bus) > 0 &&
        data->fix_flag[BUS_get_index_v_mag(bus,t)]) {
      
      VEC_set(b,*A_row,BUS_get_v_set(bus,t));
          
      // v
      MAT_set_i(A,*A_nnz,*A_row);
      MAT_set_j(A,*A_nnz,BUS_get_index_v_mag(bus,t));
      MAT_set_d(A,*A_nnz,1.);
      (*A_nnz)++;
          
      // Q
      for(REG_OBJ_init(&obj1_type,&obj1,bus); obj1 != NULL; REG_OBJ_next(&obj1_type,&obj1,bus)) {
        if (REG_OBJ_is_candidate(obj1_type,obj1)) {
          MAT_set_i(A,*A_nnz,*A_row);
          MAT_set_j(A,*A_nnz,REG_OBJ_get_index_Q(obj1_type,obj1,t));
          MAT_set_d(A,*A_nnz,0.);
          (*A_nnz)++;
        }
      }
          
      (*A_row)++;
    }
        
    // Q var and fixed
    //****************	
    for(REG_OBJ_init(&obj1_type,&obj1,bus); obj1 != NULL; REG_OBJ_next(&obj1_type,&obj1,bus)) {
      if (REG_OBJ_is_candidate(obj1_type,obj1) && data->fix_flag[REG_OBJ_get_index_Q(obj1_type,obj1,t)]) {
            
        Q = REG_OBJ_get_Q(obj1_type,obj1,t);
        Q_max = REG_OBJ_get_Q_max(obj1_type,obj1);
        Q_min = REG_OBJ_get_Q_min(obj1_type,obj1);
            
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
        for(REG_OBJ_init(&obj2_type,&obj2,bus); obj2 != NULL; REG_OBJ_next(&obj2_type,&obj2,bus)) {
          if (REG_OBJ_is_candidate(obj2_type,obj2)) {
            MAT_set_i(A,*A_nnz,*A_row);
            MAT_set_j(A,*A_nnz,REG_OBJ_get_index_Q(obj2_type,obj2,t));
            if (obj2 == obj1)
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
    //*********************
    REG_OBJ_init(&obj1_type,&obj1,bus);
    while(obj1) {
          
      // Candidate 1
      if (REG_OBJ_is_candidate(obj1_type,obj1) &&
          !data->fix_flag[REG_OBJ_get_index_Q(obj1_type,obj1,t)]) {
            
        obj2_type = obj1_type;
        obj2 = obj1;
        for(REG_OBJ_next(&obj2_type,&obj2,bus); obj2 != NULL; REG_OBJ_next(&obj2_type,&obj2,bus)) {
              
          // Candidate 2
          if (REG_OBJ_is_candidate(obj2_type,obj2) &&
              !data->fix_flag[REG_OBJ_get_index_Q(obj2_type,obj2,t)]) {
                
            VEC_set(b,*A_row,0.);
                
            alpha1 = REG_OBJ_get_Q_par(obj1_type,obj1);
            if (alpha1 < CONSTR_PVPQ_SWITCHING_PARAM)
              alpha1 = CONSTR_PVPQ_SWITCHING_PARAM;
                
            alpha2 = REG_OBJ_get_Q_par(obj2_type,obj2);
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
            for(REG_OBJ_init(&obj3_type,&obj3,bus); obj3 != NULL; REG_OBJ_next(&obj3_type,&obj3,bus)) {
              if (REG_OBJ_is_candidate(obj3_type,obj3)) {
                MAT_set_i(A,*A_nnz,*A_row);
                MAT_set_j(A,*A_nnz,REG_OBJ_get_index_Q(obj3_type,obj3,t));
                if (obj3 == obj1)
                  MAT_set_d(A,*A_nnz,alpha2);
                else if (obj3 == obj2)
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
        obj1_type = obj2_type;
        obj1 = obj2;
      }
      else {
            
        // Move forward
        REG_OBJ_next(&obj1_type,&obj1,bus);
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
