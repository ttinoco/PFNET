/** @file constr_FIX.c
 *  @brief This file defines the data structure and routines associated with the constraint of type FIX.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015-2016, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/constr_FIX.h>

void CONSTR_FIX_init(Constr* c) {

  // Init
  CONSTR_set_data(c,NULL);
}

void CONSTR_FIX_clear(Constr* c) {
  
  // Counters
  CONSTR_set_Acounter(c,0);
  CONSTR_set_Aconstr_index(c,0);

  // Flags
  CONSTR_clear_bus_counted(c);
}

void CONSTR_FIX_count_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  int* Acounter;
  int* Aconstr_index;
  char* bus_counted;
  int i;
  
  // Constr data
  Acounter = CONSTR_get_Acounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !Aconstr_index || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Tap ratio 
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    (*Acounter)++;
    (*Aconstr_index)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    (*Acounter)++;
    (*Aconstr_index)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) && 
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	(*Acounter)++;
	
	// Extra nz for regulating generator (for PV-PQ switching?)
	if (BUS_is_regulated_by_gen(bus) && 
	    GEN_has_flags(BUS_get_reg_gen(bus),FLAG_VARS,GEN_VAR_Q))
	  (*Acounter)++;
	
	(*Aconstr_index)++;
      }
	
      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) && 
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	(*Acounter)++;
	(*Aconstr_index)++;
      }
      
      // Generators	  	  
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
	
	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }
      
      // Variable generators	  	  
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) && 
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
	
	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) && 
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }
    }
      
    // Update counted flag
    bus_counted[BUS_get_index(bus)] = TRUE;    
  }
}

void CONSTR_FIX_allocate(Constr *c) {
  
  // Local variables
  int num_constr;
  int num_vars;
  int Acounter;
  
  num_vars = NET_get_num_vars(CONSTR_get_network(c));
  num_constr = CONSTR_get_Aconstr_index(c);
  Acounter = CONSTR_get_Acounter(c);

  // J f
  CONSTR_set_J(c,MAT_new(0,num_vars,0));
  CONSTR_set_f(c,VEC_new(0));

  // G l u
  CONSTR_set_G(c,MAT_new(0,num_vars,0));
  CONSTR_set_l(c,VEC_new(0));
  CONSTR_set_u(c,VEC_new(0));

  // b
  CONSTR_set_b(c,VEC_new(num_constr));

  // A
  CONSTR_set_A(c,MAT_new(num_constr, // size1 (rows)
			 num_vars,   // size2 (cols)
			 Acounter)); // nnz
}

void CONSTR_FIX_analyze_branch(Constr* c, Branch* br) {
  
  // Local variables
  Bus* buses[2];
  Bus* bus;
  Gen* gen;
  Vargen* vargen;
  Shunt* shunt;
  int* Acounter;
  int* Aconstr_index;
  char* bus_counted;
  Vec* b;
  Mat* A;
  int i;
  
  // Cosntr data
  b = CONSTR_get_b(c);
  A = CONSTR_get_A(c);
  Acounter = CONSTR_get_Acounter_ptr(c);
  Aconstr_index = CONSTR_get_Aconstr_index_ptr(c);
  bus_counted = CONSTR_get_bus_counted(c);
  if (!Acounter || !Aconstr_index || !bus_counted)
    return;

  // Bus data
  buses[0] = BRANCH_get_bus_from(br);
  buses[1] = BRANCH_get_bus_to(br);

  // Tap ratio 
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_RATIO) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_RATIO)) {
    VEC_set(b,*Aconstr_index,BRANCH_get_ratio(br));
    MAT_set_i(A,*Acounter,*Aconstr_index);
    MAT_set_j(A,*Acounter,BRANCH_get_index_ratio(br));
    MAT_set_d(A,*Acounter,1.);
    (*Acounter)++;
    (*Aconstr_index)++;
  }

  // Phase shift
  if (BRANCH_has_flags(br,FLAG_FIXED,BRANCH_VAR_PHASE) && 
      BRANCH_has_flags(br,FLAG_VARS,BRANCH_VAR_PHASE)) {
    VEC_set(b,*Aconstr_index,BRANCH_get_phase(br));
    MAT_set_i(A,*Acounter,*Aconstr_index);
    MAT_set_j(A,*Acounter,BRANCH_get_index_phase(br));
    MAT_set_d(A,*Acounter,1.);
    (*Acounter)++;
    (*Aconstr_index)++;
  }

  // Buses
  for (i = 0; i < 2; i++) {
    
    bus = buses[i];
    
    if (!bus_counted[BUS_get_index(bus)]) {

      // Voltage magnitude (V_MAG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VMAG) && 
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {
	if (BUS_is_regulated_by_gen(bus))
	  VEC_set(b,*Aconstr_index,BUS_get_v_set(bus));
	else
	  VEC_set(b,*Aconstr_index,BUS_get_v_mag(bus));
	MAT_set_i(A,*Acounter,*Aconstr_index);
	MAT_set_j(A,*Acounter,BUS_get_index_v_mag(bus));
	MAT_set_d(A,*Acounter,1.);
	(*Acounter)++;

	// Extra nz for regulating generator (for PV-PQ switching?)
	if (BUS_is_regulated_by_gen(bus) && 
	    GEN_has_flags(BUS_get_reg_gen(bus),FLAG_VARS,GEN_VAR_Q)) {
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,GEN_get_index_Q(BUS_get_reg_gen(bus)));
	  MAT_set_d(A,*Acounter,0.);
	  (*Acounter)++;
	}
	
	(*Aconstr_index)++;
      }
      
      // Voltage angle (V_ANG)
      if (BUS_has_flags(bus,FLAG_FIXED,BUS_VAR_VANG) && 
	  BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VANG)) {
	VEC_set(b,*Aconstr_index,BUS_get_v_ang(bus));
	MAT_set_i(A,*Acounter,*Aconstr_index);
	MAT_set_j(A,*Acounter,BUS_get_index_v_ang(bus));
	MAT_set_d(A,*Acounter,1.);
	(*Acounter)++;
	(*Aconstr_index)++;
      }
      
      // Generators
      for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {
	
	// Active power (P)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_P) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_P)) {
	  VEC_set(b,*Aconstr_index,GEN_get_P(gen));
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,GEN_get_index_P(gen));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (GEN_has_flags(gen,FLAG_FIXED,GEN_VAR_Q) && 
	    GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) {
	  VEC_set(b,*Aconstr_index,GEN_get_Q(gen));
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,GEN_get_index_Q(gen));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }

      // Variable generators
      for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)) {
	
	// Active power (P)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_P) && 
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_P)) {
	  VEC_set(b,*Aconstr_index,VARGEN_get_P(vargen));
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,VARGEN_get_index_P(vargen));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}

	// Reactive power (Q)
	if (VARGEN_has_flags(vargen,FLAG_FIXED,VARGEN_VAR_Q) && 
	    VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) {
	  VEC_set(b,*Aconstr_index,VARGEN_get_Q(vargen));
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,VARGEN_get_index_Q(vargen));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
	
	// Susceptance (b)
	if (SHUNT_has_flags(shunt,FLAG_FIXED,SHUNT_VAR_SUSC) && 
	    SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {
	  VEC_set(b,*Aconstr_index,SHUNT_get_b(shunt));
	  MAT_set_i(A,*Acounter,*Aconstr_index);
	  MAT_set_j(A,*Acounter,SHUNT_get_index_b(shunt));
	  MAT_set_d(A,*Acounter,1.);
	  (*Acounter)++;
	  (*Aconstr_index)++;
	}
      }
    }

    // Update counted flag
    bus_counted[BUS_get_index(bus)] = TRUE;    
  } 
}

void CONSTR_FIX_eval_branch(Constr* c, Branch *br, Vec* var_values) {
  // Nothing to do
}

void CONSTR_FIX_store_sens_branch(Constr* c, Branch *b, Vec* sA, Vec* sf, Vec* sGu, Vec* sGl) {
  // Nothing
}

void CONSTR_FIX_free(Constr* c) {
  // Nothing to do
}
