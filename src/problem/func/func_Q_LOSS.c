/** @file func_Q_LOSS.c
 *  @brief This file defines the data structure and routines associated with the function of type Q_LOSS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_Q_LOSS.h>

Func* FUNC_Q_LOSS_new(REAL weight, Net* net) {
  Func* f = FUNC_new(weight,net);
  FUNC_set_func_count_step(f,&FUNC_Q_LOSS_count_step);
  FUNC_set_func_analyze_step(f,&FUNC_Q_LOSS_analyze_step);
  FUNC_set_func_eval_step(f,&FUNC_Q_LOSS_eval_step);
  FUNC_set_name(f,"reactive power loss");
  return f;
}

void FUNC_Q_LOSS_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  int* Hphi_nnz;
  Shunt* shunt;

  // Func data
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // Bus
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {
        if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC) &&
            SHUNT_is_in_service(shunt))

          // b & v var
          (*Hphi_nnz)++;
      }

      // v var
      (*Hphi_nnz)++;
    
    }
}

void FUNC_Q_LOSS_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

  // Local variables
  Shunt* shunt;
  int* Hphi_nnz;
  Mat* Hphi;

  // Constr data
  Hphi = FUNC_get_Hphi(f);
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!Hphi_nnz || !Hphi || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // v var
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

      // Shunts
      for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

        if (!SHUNT_is_in_service(shunt))
          continue;

        if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) {

          // b & v
          MAT_set_i(Hphi,*Hphi_nnz,SHUNT_get_index_b(shunt,t));
          MAT_set_j(Hphi,*Hphi_nnz,BUS_get_index_v_mag(bus,t));
          (*Hphi_nnz)++;
        }
      }

      // Hphi  >> v & v
      MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
      MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
      (*Hphi_nnz)++;

    }
}

void FUNC_Q_LOSS_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

  // Local variables
  Branch* br;
  Shunt* shunt;
  Gen* gen;
  Load* load;
  Vargen* vargen;
  REAL* phi;
  REAL* gphi;
  REAL* Hphi;
  int* Hphi_nnz;
  int index_v_mag;
  REAL v;
  REAL tot_b = 0.0;
  REAL gphi_vmag = 0.0;

  // Func data
  phi = FUNC_get_phi_ptr(f);
  gphi = VEC_get_data(FUNC_get_gphi(f));
  Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
  Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

  // Check pointers
  if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
    return;

  // Out of service
  if (!BUS_is_in_service(bus))
    return;

  // Generators
  for (gen = BUS_get_gen(bus); gen != NULL; gen = GEN_get_next(gen)) {

    // Out of service
    if (!GEN_is_in_service(gen))
      continue;

    if (GEN_has_flags(gen,FLAG_VARS,GEN_VAR_Q)) { // Q var

      // phi
      (*phi) += VEC_get(var_values,GEN_get_index_Q(gen,t));

      // gphi
      gphi[GEN_get_index_Q(gen,t)] = 1.0;
    }
    else {

      // phi
      (*phi) += GEN_get_Q(gen,t);
    }
  }

  // Variable generators
  for (vargen = BUS_get_vargen(bus); vargen != NULL; vargen = VARGEN_get_next(vargen)){

    // Out of service
    if (!VARGEN_is_in_service(vargen))
      continue;

    if (VARGEN_has_flags(vargen,FLAG_VARS,VARGEN_VAR_Q)) { // Q var

      // phi
      (*phi) += VEC_get(var_values,VARGEN_get_index_Q(vargen,t));

      // gphi
      gphi[VARGEN_get_index_Q(vargen,t)] = 1.0;
    }
    else {

      // phi
      (*phi) += VARGEN_get_Q(vargen,t);
    }

  }

  // Loads
    for (load = BUS_get_load(bus); load != NULL; load = LOAD_get_next(load)) {

    // Out of service
    if (!LOAD_is_in_service(load))
      continue;

    // Variable
    if (LOAD_has_flags(load,FLAG_VARS,LOAD_VAR_Q)) {

      // phi
      (*phi) -= VEC_get(var_values,LOAD_get_index_Q(load,t));

      // gphi
      gphi[LOAD_get_index_Q(load,t)] = -1.0;

    }

    // Constant
    else {

      // phi
      (*phi) -= LOAD_get_Q(load,t);
    }
  }

  // Bus voltage
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

      index_v_mag = BUS_get_index_v_mag(bus,t);
      v = VEC_get(var_values,index_v_mag);
  }
  else {
      v = BUS_get_v_mag(bus,t);
  }

  // Branches_k
  for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;

    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

      // phi
      (*phi) += BRANCH_get_b_k(br) * pow(v, 2.);

      // gphi
      gphi_vmag += 2.* BRANCH_get_b_k(br) * v;

      tot_b += BRANCH_get_b_k(br);
      }

    else {

      // phi
      (*phi) += BRANCH_get_b_k(br) * pow(v, 2.);

      }
  }

  // Branches_m
  for (br = BUS_get_branch_m(bus); br != NULL; br = BRANCH_get_next_m(br)) {

    // Out of service
    if (!BRANCH_is_in_service(br))
      continue;

    if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {

      // phi
      (*phi) += BRANCH_get_b_m(br) * pow(v, 2.);

      // gphi
      gphi_vmag += 2.* BRANCH_get_b_m(br) * v;

      tot_b += BRANCH_get_b_m(br);
    }
    else {

      // phi
      (*phi) += BRANCH_get_b_m(br) * pow(v, 2.);
    }
  }

  // Shunts
  for (shunt = BUS_get_shunt(bus); shunt != NULL; shunt = SHUNT_get_next(shunt)) {

    if (!SHUNT_is_in_service(shunt))
      continue;

    if (SHUNT_has_flags(shunt,FLAG_VARS,SHUNT_VAR_SUSC)) { // b var

        // phi
        (*phi) += VEC_get(var_values,SHUNT_get_index_b(shunt,t)) * pow(v,2.);

        // gphi > shunt
        gphi[SHUNT_get_index_b(shunt,t)] = pow(v,2.);


            if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {  // V var

                 // gphi > vmag
                gphi_vmag += 2. * VEC_get(var_values,SHUNT_get_index_b(shunt,t)) * v;

                // Hphi > v & b
                Hphi[*Hphi_nnz] = 2. * v;
                (*Hphi_nnz)++;

                tot_b += VEC_get(var_values,SHUNT_get_index_b(shunt,t));
            }

    }
    else{  // b const

        // phi
        (*phi) += SHUNT_get_b(shunt,t) * pow(v,2.);

        if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)) {  // V var

            // gphi
            gphi_vmag += 2. * SHUNT_get_b(shunt,t) * v;

            tot_b += SHUNT_get_b(shunt,t);
        }
    }
  }

  // Hessian for bus
  if (BUS_has_flags(bus,FLAG_VARS,BUS_VAR_VMAG)){

    // gphi
    gphi[index_v_mag] = gphi_vmag;

    // Hphi
    Hphi[*Hphi_nnz] = 2. * tot_b;
    (*Hphi_nnz)++;
  }
}
