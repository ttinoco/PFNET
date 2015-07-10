/** @file problem.c
 *  @brief This file defines the Prob data structure and its associated methods.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/problem.h>

struct Prob {

  // Error
  BOOL error_flag;                     /**< @brief Error flags */
  char error_string[PROB_BUFFER_SIZE]; /**< @brief Error string */

  // Constraints
  Constr* constr; /**< @brief List of constraints */

  // Functions
  Func* func;     /**< @brief List of objective functions */

  // Heurisics
  Heur* heur;     /**< @brief List of heuristics */
  
  // Network
  Net* net;       /**< @brief Power flow network */

  // Point
  Vec* point;     /**< @brief Vector of problem variables */

  // Objective function
  REAL phi;       /**< @brief Combined objective value */
  Vec* gphi;      /**< @brief Gradient of combined objective function */
  Mat* Hphi;      /**< @brief Hessian of combined objective function */

  // Linear equality constraints
  Vec* b;         /**< @brief Combined right-hand side of linear constraints */
  Mat* A;         /**< @brief Combined linear constraint matrix */
  Mat* Z;         /**< @brief Null-space matrix of linear constraints */
  
  // Nonlinear equality constraints
  Vec* f;           /**< @brief Nonlinear constraint function values */
  Mat* J;           /**< @brief Jacobian matrix of nonlinear constraints */
  Mat* H_combined;  /**< @brief Combined Hessians of nonlinear constraints */
};

void PROB_add_constr(Prob* p, int type) {
  if (p) {
    if (!PROB_find_constr(p,type))
      p->constr = CONSTR_list_add(p->constr,CONSTR_new(type,p->net));
  }
}

void PROB_add_func(Prob* p, int type, REAL weight) {
  if (p)
    p->func = FUNC_list_add(p->func,FUNC_new(type,weight,p->net));
}

void PROB_add_heur(Prob* p, int type) {
  if (p)
    p->heur = HEUR_list_add(p->heur,
			    HEUR_new(type,p->net));
}

void PROB_analyze(Prob* p) {

  // Local variables
  Branch* br;
  Constr* c;
  Func* f;
  int Arow;
  int Annz;
  int Jrow;
  int Jnnz;
  int Hphinnz;
  int Hcombnnz;
  int num_vars;
  int i;
  
  if (!p)
    return;

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  
  // Count
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_count_branch(p->constr,br);
    FUNC_list_count_branch(p->func,br);
  }
    
  // Allocate
  CONSTR_list_allocate(p->constr);
  FUNC_list_allocate(p->func);

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);

  // Analyze
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_analyze_branch(p->constr,br);
    FUNC_list_analyze_branch(p->func,br);
  }

  // Allocate problem matvec
  Annz = 0;
  Arow = 0;
  Jrow = 0;
  Jnnz = 0;
  Hphinnz = 0;
  Hcombnnz = 0;
  num_vars = NET_get_num_vars(p->net);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    Arow += MAT_get_size1(CONSTR_get_A(c));
    Annz += MAT_get_nnz(CONSTR_get_A(c));
    Jrow += MAT_get_size1(CONSTR_get_J(c));
    Jnnz += MAT_get_nnz(CONSTR_get_J(c));
    Hcombnnz += MAT_get_nnz(CONSTR_get_H_combined(c));
  }
  for (f = p->func; f != NULL; f = FUNC_get_next(f))
    Hphinnz += MAT_get_nnz(FUNC_get_Hphi(f));
  p->phi = 0;
  p->gphi = VEC_new(num_vars);
  p->Hphi = MAT_new(num_vars,num_vars,Hphinnz);
  p->b = VEC_new(Arow);
  p->A = MAT_new(Arow,num_vars,Annz);
  p->f = VEC_new(Jrow);
  p->J = MAT_new(Jrow,num_vars,Jnnz);
  p->H_combined = MAT_new(num_vars,num_vars,Hcombnnz);

  // Update 
  PROB_update_lin(p);
  PROB_update_nonlin_struc(p);

  // Construct Z
  PROB_construct_Z(p);
}

void PROB_apply_heuristics(Prob* p, Vec* point) {

  // Local variables
  Branch* br;
  Constr* c;
  int i;
  
  if (!p)
    return;

  // Clear
  HEUR_list_clear(p->heur,p->net);

  // Apply
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    HEUR_list_apply_to_branch(p->heur,p->constr,p->net,br,point);
  }

  // Udpate A and b
  PROB_update_lin(p);
}

void PROB_eval(Prob* p, Vec* point) {

  // Local variables
  Branch* br;
  int i;
  
  if (!p)
    return;

  // Clear
  CONSTR_list_clear(p->constr);
  FUNC_list_clear(p->func);
  NET_clear_properties(p->net);

  // Eval
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_eval_branch(p->constr,br,point);
    FUNC_list_eval_branch(p->func,br,point);
    NET_update_properties_branch(p->net,br,point);
  }

  // Update 
  PROB_update_nonlin_data(p);
}

void PROB_store_sens(Prob* p, Vec* sens) {

  // Local variables
  Branch* br;
  int i;
  
  if (!p || !sens)
    return;

  // Check size
  if (VEC_get_size(sens) != VEC_get_size(p->f)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }

  // Clear
  CONSTR_list_clear(p->constr);

  // Store sens
  for (i = 0; i < NET_get_num_branches(p->net); i++) {
    br = NET_get_branch(p->net,i);
    CONSTR_list_store_sens_branch(p->constr,br,sens);
  }
}

void PROB_del(Prob* p) {
  if (p) {
    PROB_clear(p);
    free(p);
    p = NULL;
  }
}

void PROB_clear(Prob* p) {
  if (p) {
    
    // Free data
    CONSTR_list_del(p->constr);
    FUNC_list_del(p->func);
    HEUR_list_del(p->heur);
    VEC_del(p->point);
    VEC_del(p->b);
    MAT_del(p->A);
    MAT_del(p->Z);
    VEC_del(p->f);
    MAT_del(p->J);
    MAT_del(p->H_combined);
    VEC_del(p->gphi);
    MAT_del(p->Hphi);    

    // Re-initialize
    PROB_init(p);
  }
}

void PROB_combine_H(Prob* p, Vec* coeff, BOOL ensure_psd) {
  
  // Local variables
  Constr* c;
  int Hcombnnz;
  REAL* Hcomb;
  REAL* Hcomb_constr;
  int i;
  
  if (!p || !coeff)
    return;

  // Check size
  if (VEC_get_size(coeff) != VEC_get_size(p->f)) {
    sprintf(p->error_string,"invalid vector size");
    p->error_flag = TRUE;
    return;
  }
  
  // Combine
  CONSTR_list_combine_H(p->constr,coeff,ensure_psd);

  // Combine and update
  Hcombnnz = 0;
  Hcomb = MAT_get_data_array(p->H_combined);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    Hcomb_constr = MAT_get_data_array(CONSTR_get_H_combined(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_H_combined(c)); i++) {
      Hcomb[Hcombnnz] = Hcomb_constr[i];
      Hcombnnz++;
    }
  }
}

void PROB_construct_Z(Prob* p) {

  // Local variables
  Constr* c;
  Bus* bus;
  Gen* gen;
  Gen* gen1;
  Gen* gen2;
  Shunt* shunt;
  Branch* branch;
  BOOL has_REG_GEN = FALSE;
  BOOL has_REG_TRAN = FALSE;
  BOOL has_REG_SHUNT = FALSE;
  BOOL has_FIX = FALSE;
  BOOL has_PAR = FALSE;
  int Znnz;
  int Zcol;
  int i;
  int num_free_P;
  int num_free_Q;
  int n;
  int m;
  int nz;
  REAL dQ;
  REAL max_dQ;
  
  if (!p || !p->A || !p->net)
    return;

  // Set flags
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
    if (CONSTR_get_type(c) == CONSTR_TYPE_FIX)
      has_FIX = TRUE;
    if (CONSTR_get_type(c) == CONSTR_TYPE_REG_GEN)
      has_REG_GEN = TRUE;
    if (CONSTR_get_type(c) == CONSTR_TYPE_PAR_GEN)
      has_PAR = TRUE;
    if (CONSTR_get_type(c) == CONSTR_TYPE_REG_TRAN)
      has_REG_TRAN = TRUE;
    if (CONSTR_get_type(c) == CONSTR_TYPE_REG_SHUNT)
      has_REG_SHUNT = TRUE;
  }

  // Count
  //******
  Znnz = 0;

  // Buses
  for (i = 0; i < NET_get_num_buses(p->net); i++) { 
    
    bus = NET_get_bus(p->net,i);

    // w
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG) &&
	(!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG))) // var and not fixed
      Znnz++; // w free

    // v
    if (has_REG_GEN &&
	BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	BUS_is_regulated_by_gen(bus) &&
	!BUS_is_slack(bus) &&
	BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // reg gen constraint
      
      if (!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG)) { // v not fixed
	Znnz++;
	Znnz++; // v and y coupled
	
	Znnz++;
	Znnz++; // v and z coupled
      }
      else { // v fixed
	Znnz++;
	Znnz++; // y and z coupled
      }
    }
    else { // no reg gen constraint
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	  (!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG))) // v var and not fixed
	Znnz++;  // v free
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	Znnz++;  // y free
	Znnz++;  // z free
      }
    }
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
      Znnz++; // vl free
      Znnz++; // vh free
    }

    // slack gens P
    if (has_PAR && BUS_is_slack(bus)) {

      // number of free
      num_free_P = 0;
      for (gen1 = BUS_get_gen(bus); gen1 != NULL; gen1 = GEN_get_next(gen1)) {
	if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P) &&
	    (!has_FIX || !GEN_has_flags(gen1,FLAG_FIXED,GEN_VAR_P))) // var and not fixed
	  num_free_P++;
      }

      gen1 = BUS_get_gen(bus);
      if (num_free_P == GEN_list_len(gen1)) { // all var and not fixed	
	for (gen2 = BUS_get_gen(bus); gen2 != NULL; gen2 = GEN_get_next(gen2))
	  Znnz++; // P coupled
      }
    }
    
    // reg gen Q
    if (has_PAR && BUS_is_regulated_by_gen(bus)) {

      // number of free
      num_free_Q = 0;
      for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
	if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) &&
	    (!has_FIX || !GEN_has_flags(gen1,FLAG_FIXED,GEN_VAR_Q))) // var and not fixed
	  num_free_Q++;
      }

      gen1 = BUS_get_reg_gen(bus);
      if (num_free_Q == GEN_list_reg_len(gen1)) { // all var and not fixed
	for (gen2 = BUS_get_reg_gen(bus); gen2 != NULL; gen2 = GEN_get_reg_next(gen2))
	  Znnz++; // Q coupled
      }
    }    
  }
  
  // Other gens
  for (i = 0; i < NET_get_num_gens(p->net); i++) {
    
    gen = NET_get_gen(p->net,i);

    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) &&
	(!has_FIX || !GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P)) &&
	(!has_PAR || !GEN_is_slack(gen)))
      Znnz++; // P free
    
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q) &&
	(!has_FIX || !GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q)) &&
	(!has_PAR || !GEN_is_regulator(gen)))
      Znnz++; // Q free
  }

  // Branches
  for (i = 0; i < NET_get_num_branches(p->net); i++) {

    branch = NET_get_branch(p->net,i);
    
    if (has_REG_TRAN &&
	BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO) &&
	BRANCH_is_tap_changer_v(branch) &&
	BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // reg tran constraint
      
      if (!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_RATIO)) { // ratio not fixed
	Znnz++; 
	Znnz++; // ratio and ratio_y coupled

	Znnz++; 
	Znnz++; // ratio and ratio_z coupled
      }
      else { // ratio fixed
	Znnz++; 
	Znnz++; // ratio_y and ratio_z coupled
      }
    }
    else { // no reg tran constraints
      if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO) &&
	  (!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_RATIO))) // ratio var and not fixed
	Znnz++;  // ratio free
      if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) {
	Znnz++; // ratio_y free
	Znnz++; // ratio_z free
      }
    }
    
    if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_PHASE) &&
	(!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_PHASE)))
      Znnz++; // phase free
  }

  // Shunts
  for (i = 0; i < NET_get_num_shunts(p->net); i++) {
    
    shunt = NET_get_shunt(p->net,i);

    if (has_REG_SHUNT &&
	SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	SHUNT_is_switched_v(shunt) &&
	SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // reg shunt constraint
      
      if (!has_FIX || !SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC)) { // susc not fixed
	Znnz++; 
	Znnz++; // susc and susc_y coupled
	
	Znnz++; 
	Znnz++; // susc and susc_z coupled
      }
      else { // susc fixed
	Znnz++; 
	Znnz++; // susc_y and susc_z coupled
      }
    }
    else { // no reg shunt constraint
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	  (!has_FIX || !SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC))) // susc var and not fixed
	Znnz++; // susceptance free
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // susc deviations var and not fixed
	Znnz++; // susc_y free
	Znnz++; // susc_z free
      }
    }
  }

  // Allocate
  //*********
  m = MAT_get_size1(p->A); // constraints
  n = MAT_get_size2(p->A); // variables
  nz = n-m >= 0 ? n-m : 0; // rank of null space (assuming A is full rank)
  p->Z = MAT_new(n,nz,Znnz);
		 
  // Fill
  //*****
  Zcol = 0;
  Znnz = 0;

  // Buses
  for (i = 0; i < NET_get_num_buses(p->net); i++) {
    
    bus = NET_get_bus(p->net,i);

    // w
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG) &&
	(!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG))) { // var and not fixed
      MAT_set_i(p->Z,Znnz,BUS_get_index_v_ang(bus));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // w free
      Zcol++;
    }

    // v	
    if (has_REG_GEN &&                    
	BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	BUS_is_regulated_by_gen(bus) &&
	!BUS_is_slack(bus) &&
	BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) { // reg gen constraint
      
      if (!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG)) { // v not fixed
	MAT_set_i(p->Z,Znnz,BUS_get_index_v_mag(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++;
	MAT_set_i(p->Z,Znnz,BUS_get_index_y(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // v y coupled
	Zcol++;

	MAT_set_i(p->Z,Znnz,BUS_get_index_v_mag(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++;
	MAT_set_i(p->Z,Znnz,BUS_get_index_z(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,-1.);
	Znnz++; // v z coupled
	Zcol++;
      }
      else { // v fixed
	MAT_set_i(p->Z,Znnz,BUS_get_index_y(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++;
	MAT_set_i(p->Z,Znnz,BUS_get_index_z(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // y z coupled
	Zcol++;
      }
    }
    else { // no reg gen constraint

      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG) &&
	  (!has_FIX || !BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG))) { // v var and not fixed
	MAT_set_i(p->Z,Znnz,BUS_get_index_v_mag(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // v free
	Zcol++;
      } 
      if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VDEV)) {
	MAT_set_i(p->Z,Znnz,BUS_get_index_y(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // y free
	Zcol++;

	MAT_set_i(p->Z,Znnz,BUS_get_index_z(bus));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // z free
	Zcol++;
      }
    }
    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VVIO)) {
      MAT_set_i(p->Z,Znnz,BUS_get_index_vl(bus));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // vl free
      Zcol++;

      MAT_set_i(p->Z,Znnz,BUS_get_index_vh(bus));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // vh free
      Zcol++;
    }

    // slack gens P
    if (has_PAR && BUS_is_slack(bus)) {

      // number of free
      num_free_P = 0;
      for (gen1 = BUS_get_gen(bus); gen1 != NULL; gen1 = GEN_get_next(gen1)) {
	if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_P) &&
	    (!has_FIX || !GEN_has_flags(gen1,FLAG_FIXED,GEN_VAR_P))) // var and not fixed
	  num_free_P++;
      }

      gen1 = BUS_get_gen(bus);
      if (num_free_P == GEN_list_len(gen1)) { // all var and not fixed	
	for (gen2 = BUS_get_gen(bus); gen2 != NULL; gen2 = GEN_get_next(gen2)) {
	  MAT_set_i(p->Z,Znnz,GEN_get_index_P(gen2));
	  MAT_set_j(p->Z,Znnz,Zcol);
	  MAT_set_d(p->Z,Znnz,1.);
	  Znnz++; // P coupled
	}
	Zcol++;
      }
    }
    
    // reg gen Q
    if (has_PAR && BUS_is_regulated_by_gen(bus)) {

      // number of free
      max_dQ = 0;
      num_free_Q = 0;
      for (gen1 = BUS_get_reg_gen(bus); gen1 != NULL; gen1 = GEN_get_reg_next(gen1)) {
	if (GEN_has_flags(gen1,FLAG_VARS,GEN_VAR_Q) &&
	    (!has_FIX || !GEN_has_flags(gen1,FLAG_FIXED,GEN_VAR_Q))) { // var and not fixed
	  num_free_Q++;
	  dQ = GEN_get_Q_max(gen1)-GEN_get_Q_min(gen1);
	  max_dQ = (dQ > max_dQ) ? dQ : max_dQ;
	}
      }

      gen1 = BUS_get_reg_gen(bus);
      if (num_free_Q == GEN_list_reg_len(gen1)) { // all var and not fixed
	for (gen2 = BUS_get_reg_gen(bus); gen2 != NULL; gen2 = GEN_get_reg_next(gen2)) {
	  MAT_set_i(p->Z,Znnz,GEN_get_index_Q(gen2));
	  MAT_set_j(p->Z,Znnz,Zcol);
	  MAT_set_d(p->Z,Znnz,(GEN_get_Q_max(gen2)-GEN_get_Q_min(gen2))/max_dQ);
	  Znnz++; // Q coupled
	}
	Zcol++;
      }
    }    
  }
  
  // Other gens
  for (i = 0; i < NET_get_num_gens(p->net); i++) {
    
    gen = NET_get_gen(p->net,i);

    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P) &&
	(!has_FIX || !GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P)) &&
	(!has_PAR || !GEN_is_slack(gen))) {
      MAT_set_i(p->Z,Znnz,GEN_get_index_P(gen));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // P free
      Zcol++;
    }
    
    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q) &&
	(!has_FIX || !GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q)) &&
	(!has_PAR || !GEN_is_regulator(gen))) {
      MAT_set_i(p->Z,Znnz,GEN_get_index_Q(gen));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // Q free
      Zcol++;
    }
  }

  // Branches
  for (i = 0; i < NET_get_num_branches(p->net); i++) {

    branch = NET_get_branch(p->net,i);

    if (has_REG_TRAN &&
	BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO) &&
	BRANCH_is_tap_changer_v(branch) &&
	BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) { // reg tran constraint
      
      if (!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_RATIO)) { // ratio not fixed
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_y(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // ratio and ratio_y coupled
	Zcol++;

	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_z(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,-1.);
	Znnz++; // ratio and ratio_z coupled
	Zcol++;
      }
      else { // ratio fixed
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_y(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_z(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // ratio_y and ratio_z coupled
	Zcol++;
      }
    }
    else { // no reg tran constraint
      if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO) &&
	  (!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_RATIO))) { // ratio var and not fixed
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++;  // ratio free
	Zcol++;
      }
      if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_RATIO_DEV)) {
	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_y(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // ratio_y free
	Zcol++;

	MAT_set_i(p->Z,Znnz,BRANCH_get_index_ratio_z(branch));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // ratio_z free
	Zcol++;
      }
    }
        
    if (BRANCH_has_flags(branch,FLAG_VARS,BRANCH_VAR_PHASE) &&
	(!has_FIX || !BRANCH_has_flags(branch,FLAG_FIXED,BRANCH_VAR_PHASE))) {
      MAT_set_i(p->Z,Znnz,BRANCH_get_index_phase(branch));
      MAT_set_j(p->Z,Znnz,Zcol);
      MAT_set_d(p->Z,Znnz,1.);
      Znnz++; // phase free
      Zcol++;
    }
  }
  
  // Shunts
  for (i = 0; i < NET_get_num_shunts(p->net); i++) {
    
    shunt = NET_get_shunt(p->net,i);

    if (has_REG_SHUNT &&
	SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	SHUNT_is_switched_v(shunt) &&
	SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) { // reg shunt constraint
      
      if (!has_FIX || !SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC)) { // susc not fixed
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_b(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_y(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // susc and susc_y coupled
	Zcol++;

	MAT_set_i(p->Z,Znnz,SHUNT_get_index_b(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_z(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,-1.);
	Znnz++; // susc and susc_z coupled
	Zcol++;
      }
      else { // susc fixed
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_y(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; 
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_z(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // susc_y and susc_z coupled
	Zcol++;
      }
    }
    else {
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
	  (!has_FIX || !SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC))) {
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_b(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // susceptance free
	Zcol++;
      }
      if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC_DEV)) {
	MAT_set_i(p->Z,Znnz,SHUNT_get_index_y(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // susc_y free
	Zcol++;

	MAT_set_i(p->Z,Znnz,SHUNT_get_index_z(shunt));
	MAT_set_j(p->Z,Znnz,Zcol);
	MAT_set_d(p->Z,Znnz,1.);
	Znnz++; // susc_z free
	Zcol++;
      }
    }
  }

  // Check counters
  if (Znnz != MAT_get_nnz(p->Z) || Zcol != nz) {
    sprintf(p->error_string,"error while constructing Z");
    p->error_flag = TRUE;
  }
}

Constr* PROB_find_constr(Prob* p, int constr_type) {
  Constr* c;
  if (p) {
    for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {
      if (CONSTR_get_type(c) == constr_type)
	return c;
    }
    return NULL;
  }
  else
    return NULL;
}

Constr* PROB_get_constr(Prob* p) {
  if (p)
    return p->constr;
  else
    return NULL;
}

char* PROB_get_error_string(Prob* p) {
  if (!p)
    return NULL;
  else
    return p->error_string;
}

Func* PROB_get_func(Prob* p) {
  if (p)
    return p->func;
  else
    return NULL;
}

Heur* PROB_get_heur(Prob* p) {
  if (p)
    return p->heur;
  else
    return NULL;
}

Vec* PROB_get_init_point(Prob* p) {
  if (p) {
    p->point = NET_get_var_values(p->net);
    return p->point;
  } 
  return NULL;
}

Net* PROB_get_network(Prob* p) {
  if (p)
    return p->net;
  else
    return NULL;
}

REAL PROB_get_phi(Prob* p) {
  if (p)
    return p->phi;
  else
    return 0;
}

Vec* PROB_get_gphi(Prob* p) {
  if (p)
    return p->gphi;
  else
    return NULL;
}

Mat* PROB_get_Hphi(Prob* p) {
  if (p)
    return p->Hphi;
  else
    return NULL;
}

Mat* PROB_get_A(Prob* p) {
  if (p)
    return p->A;
  else 
    return NULL;
}

Mat* PROB_get_Z(Prob* p) {
  if (p)
    return p->Z;
  else 
    return NULL;
}

Vec* PROB_get_b(Prob* p) {
  if (p)
    return p->b;
  else 
    return NULL;
}

Mat* PROB_get_J(Prob* p) {
  if (p)
    return p->J;
  else 
    return NULL;
}

Vec* PROB_get_f(Prob* p) {
  if (p)
    return p->f;
  else 
    return NULL;
}

Mat* PROB_get_H_combined(Prob* p) {
  if (p)
    return p->H_combined;
  else
    return NULL;
}

BOOL PROB_has_error(Prob* p) {
  if (!p)
    return FALSE;
  else
    return p->error_flag;
}

void PROB_init(Prob* p) {
  if (p) {

    // Error
    p->error_flag = FALSE;
    strcpy(p->error_string,"");
    
    p->constr = NULL;
    p->func = NULL;
    p->heur = NULL;
    p->point = NULL;
    p->net = NULL;
    p->phi = 0;
    p->gphi = NULL;
    p->Hphi = NULL;
    p->b = NULL;
    p->A = NULL;
    p->Z = NULL;
    p->f = NULL;
    p->J = NULL;
    p->H_combined = NULL;
  }
}

Prob* PROB_new(void) {
  Prob* p = (Prob*)malloc(sizeof(Prob));
  PROB_init(p);
  return p;
}

void PROB_set_network(Prob* p, Net* net) {
  if (p)
    p->net = net;
}

void PROB_show(Prob* p) {
  if (p) {
    printf("\nProblem\n");
    printf("functions  : %d\n",FUNC_list_len(p->func));
    printf("constraints: %d\n",CONSTR_list_len(p->constr));  
    printf("heuristics : %d\n",HEUR_list_len(p->heur));
  }
}

void PROB_update_nonlin_struc(Prob* p) {
  
  // Local variables
  Func* f;
  Constr* c;
  int* Hphi_i;
  int* Hphi_j;
  int* Hphi_i_func;
  int* Hphi_j_func;
  int Hphinnz;
  int* Ji;
  int* Jj;
  int* Ji_constr;
  int* Jj_constr;
  int Jrow;
  int Jnnz;
  int* Hcomb_i;
  int* Hcomb_j;
  int* Hcomb_i_constr;
  int* Hcomb_j_constr;
  int Hcombnnz;
  int i;

  // Hphi
  Hphinnz = 0;
  Hphi_i = MAT_get_row_array(p->Hphi);
  Hphi_j = MAT_get_col_array(p->Hphi);
  for (f = p->func; f != NULL; f = FUNC_get_next(f)) {
    Hphi_i_func = MAT_get_row_array(FUNC_get_Hphi(f));
    Hphi_j_func = MAT_get_col_array(FUNC_get_Hphi(f));
    for (i = 0; i < MAT_get_nnz(FUNC_get_Hphi(f)); i++) {
      Hphi_i[Hphinnz] = Hphi_i_func[i];
      Hphi_j[Hphinnz] = Hphi_j_func[i];
      Hphinnz++;
    }     
  }

  // J and Hcomb
  Jnnz = 0;
  Jrow = 0;
  Ji = MAT_get_row_array(p->J);
  Jj = MAT_get_col_array(p->J);
  Hcombnnz = 0;
  Hcomb_i = MAT_get_row_array(p->H_combined);
  Hcomb_j = MAT_get_col_array(p->H_combined);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // J
    Ji_constr = MAT_get_row_array(CONSTR_get_J(c));
    Jj_constr = MAT_get_col_array(CONSTR_get_J(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_J(c)); i++) {
      Ji[Jnnz] = Ji_constr[i]+Jrow;
      Jj[Jnnz] = Jj_constr[i];
      Jnnz++;
    }
    Jrow += MAT_get_size1(CONSTR_get_J(c));

    // H comb
    Hcomb_i_constr = MAT_get_row_array(CONSTR_get_H_combined(c));
    Hcomb_j_constr = MAT_get_col_array(CONSTR_get_H_combined(c));
    for (i = 0; i < MAT_get_nnz(CONSTR_get_H_combined(c)); i++) {
      Hcomb_i[Hcombnnz] = Hcomb_i_constr[i];
      Hcomb_j[Hcombnnz] = Hcomb_j_constr[i];
      Hcombnnz++;
    }
  }
}

void PROB_update_nonlin_data(Prob* p) {
  
  // Local variables
  Func* func;
  Constr* c;
  REAL* gphi;
  REAL* gphi_func;
  REAL* Hphi;
  REAL* Hphi_func;
  REAL weight;
  int Hphinnz;
  REAL* f;
  REAL* f_constr;
  REAL* J;
  REAL* J_constr;
  int Jnnz;
  int Jrow;
  int i;
  
  // phi and derivatives
  p->phi = 0;
  Hphinnz = 0;
  VEC_set_zero(p->gphi);
  gphi = VEC_get_data(p->gphi);
  Hphi = MAT_get_data_array(p->Hphi);
  for (func = p->func; func != NULL; func = FUNC_get_next(func)) {

    // Weight
    weight = FUNC_get_weight(func);

    // phi
    p->phi += weight*FUNC_get_phi(func);

    //gphi
    gphi_func = VEC_get_data(FUNC_get_gphi(func));
    for (i = 0; i < NET_get_num_vars(p->net); i++)
      gphi[i] += weight*gphi_func[i];

    // Hphi
    Hphi_func = MAT_get_data_array(FUNC_get_Hphi(func));
    for (i = 0; i < MAT_get_nnz(FUNC_get_Hphi(func)); i++) {
      Hphi[Hphinnz] = weight*Hphi_func[i];
      Hphinnz++;
    }     
  }

  // f and derivatives
  Jnnz = 0;
  Jrow = 0;
  f = VEC_get_data(p->f);
  J = MAT_get_data_array(p->J);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // J and f of constraint
    f_constr = VEC_get_data(CONSTR_get_f(c));
    J_constr = MAT_get_data_array(CONSTR_get_J(c));

    // Update J
    for (i = 0; i < MAT_get_nnz(CONSTR_get_J(c)); i++) {
      J[Jnnz] = J_constr[i];
      Jnnz++;
    }

    // Update f
    for (i = 0; i < MAT_get_size1(CONSTR_get_J(c)); i++) {
      f[Jrow] = f_constr[i];
      Jrow++;
    }
  }
}

void PROB_update_lin(Prob* p) {
  
  // Local variables
  Constr* c;
  REAL* b;
  REAL* b_constr;
  int* Ai;
  int* Aj;
  REAL* Ad;
  int* Ai_constr;
  int* Aj_constr;
  REAL* Ad_constr;
  int Annz;
  int Arow;
  int i;
  
  if (!p)
    return;

  // Update
  Annz = 0;
  Arow = 0;
  b = VEC_get_data(p->b);
  Ai = MAT_get_row_array(p->A);
  Aj = MAT_get_col_array(p->A);
  Ad = MAT_get_data_array(p->A);
  for (c = p->constr; c != NULL; c = CONSTR_get_next(c)) {

    // A and b of constraint
    b_constr = VEC_get_data(CONSTR_get_b(c));
    Ai_constr = MAT_get_row_array(CONSTR_get_A(c));
    Aj_constr = MAT_get_col_array(CONSTR_get_A(c));
    Ad_constr = MAT_get_data_array(CONSTR_get_A(c));

    // Update A
    for (i = 0; i < MAT_get_nnz(CONSTR_get_A(c)); i++) {
      Ai[Annz] = Ai_constr[i]+Arow;
      Aj[Annz] = Aj_constr[i];
      Ad[Annz] = Ad_constr[i];
      Annz++;
    }

    // Update b
    for (i = 0; i < MAT_get_size1(CONSTR_get_A(c)); i++) {
      b[Arow] = b_constr[i];
      Arow++;
    }
  }
}

