/** @file func_REG_VAR.c
 *  @brief This file defines the data structure and routines associated with the function of type REG_VAR.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/array.h>
#include <pfnet/func_REG_VAR.h>

struct Func_REG_VAR_Data {
  
  REAL* x0; // center
  REAL* w;  // weights
  int num_vars;
  unsigned long int state_tag;
};

Func* FUNC_REG_VAR_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_init(f,&FUNC_REG_VAR_init);
  FUNC_set_func_count_step(f,&FUNC_REG_VAR_count_step);
  FUNC_set_func_allocate(f,&FUNC_REG_VAR_allocate);
  FUNC_set_func_analyze_step(f,&FUNC_REG_VAR_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_REG_VAR_eval_step);
  FUNC_set_func_free(f,&FUNC_REG_VAR_free);
  FUNC_set_func_set_parameter(f,&FUNC_REG_VAR_set_parameter);
  FUNC_set_name(f,"variable regularization");
  FUNC_init(f);
  return f;
}

void FUNC_REG_VAR_init(Func* f) {

  // Local variables
  Func_REG_VAR_Data* data;
  int num_vars = NET_get_num_vars(FUNC_get_network(f));
  
  // Init
  data = (Func_REG_VAR_Data*)malloc(sizeof(Func_REG_VAR_Data));
  ARRAY_zalloc(data->x0,REAL,num_vars);
  ARRAY_zalloc(data->w,REAL,num_vars);
  data->num_vars = num_vars;
  data->state_tag = FUNC_get_state_tag(f);
  FUNC_set_data(f,data);
}

void FUNC_REG_VAR_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Func_REG_VAR_Data* data;
  
  // Func data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  data = (Func_REG_VAR_Data*)FUNC_get_data(f);

  // Check pointers
  if (!Hphi_nnz || !data)
    return;
  
  // Count
  if (BUS_get_index(bus) == 0 && !busdc && t == 0)
    (*Hphi_nnz) += NET_get_num_vars(FUNC_get_network(f));
}

void FUNC_REG_VAR_allocate(Func* f) {
  
  // Local vars
  Func_REG_VAR_Data* data = (Func_REG_VAR_Data*)FUNC_get_data(f);
  int num_vars = NET_get_num_vars(FUNC_get_network(f));

  // Check
  if (!data)
    return;

  // Alloc
  if (data->num_vars != num_vars ||
      data->state_tag != FUNC_get_state_tag(f)) {
    free(data->x0);
    free(data->w);
    ARRAY_zalloc(data->x0,REAL,num_vars);
    ARRAY_zalloc(data->w,REAL,num_vars);
    data->num_vars = num_vars;
    data->state_tag = FUNC_get_state_tag(f);
  }  
}

void FUNC_REG_VAR_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Mat* Hphi;
  Func_REG_VAR_Data* data;
  int i;
  
  // Func data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  data = (Func_REG_VAR_Data*)FUNC_get_data(f);

  // Check pointers
  if (!Hphi_nnz || !data || !(data->w) || !Hphi)
    return;

  // Analyze
  if (BUS_get_index(bus) == 0 && !busdc && t == 0) {
    for (i = 0; i < NET_get_num_vars(FUNC_get_network(f)); i++) {
      MAT_set_i(Hphi,*Hphi_nnz,i);
      MAT_set_j(Hphi,*Hphi_nnz,i);
      (*Hphi_nnz)++;
    }
  }
}


void FUNC_REG_VAR_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Func_REG_VAR_Data* data;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int i;
  REAL x;
  
  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);
  data = (Func_REG_VAR_Data*)FUNC_get_data(f);

  // Check pointers
  if (!phi || !gphi || !data || !(data->w) || !(data->x0) || !Hphi || !Hphi_nnz)
    return;

  // Eval
  if (BUS_get_index(bus) == 0 && !busdc && t == 0) {
    for (i = 0; i < NET_get_num_vars(FUNC_get_network(f)); i++) {
      x = VEC_get(var_values,i);
      (*phi) += data->w[i]*pow(x-data->x0[i],2.);
      gphi[i] = 2.*data->w[i]*(x-data->x0[i]);
      Hphi[*Hphi_nnz] = 2.*data->w[i];
      (*Hphi_nnz)++;
    }
  }
}

void FUNC_REG_VAR_free(Func* f) {

  // Local variables
  Func_REG_VAR_Data* data;

  // Get data
  data = (Func_REG_VAR_Data*)FUNC_get_data(f);

  // Free
  if (data) {
    if (data->x0)
      free(data->x0);
    if (data->w)
      free(data->w);
    free(data);
  }

  // Set data
  FUNC_set_data(f,NULL);
}

void FUNC_REG_VAR_set_parameter(Func* f, char* key, void* value) {

  // Local variables
  REAL* w;
  REAL* x0;
  Net* net = FUNC_get_network(f);
  Func_REG_VAR_Data* data = (Func_REG_VAR_Data*)FUNC_get_data(f);
  int num_vars = NET_get_num_vars(net);

  // Check
  if (!data)
    return;

  // Bad num vars
  if (data->num_vars != num_vars) {
    FUNC_set_error(f,"network num vars changed");
    return;
  }

  // Bad network state
  if (data->state_tag != FUNC_get_state_tag(f)) {
    FUNC_set_error(f,"network has changed");
    return;
  }

  // Set 
  if (strcmp(key,"w") == 0) { // w
    w = (REAL*)value;
    memcpy(data->w,w,num_vars*sizeof(REAL));
  }
	
  else if (strcmp(key,"x0") == 0) { // x0
    x0 = (REAL*)value;
    memcpy(data->x0,x0,num_vars*sizeof(REAL));
  }
  
  else // unknown
    FUNC_set_error(f,"invalid parameter");
}
