/** @file func_BRN_QLOSS.c
 *  @brief This file defines the data structure and routines associated with the function of type BRN_QLOSS.
 *
 * This file is part of PFNET.
 *
 * Copyright (c) 2015, Tomas Tinoco De Rubira.
 *
 * PFNET is released under the BSD 2-clause license.
 */

#include <pfnet/func_BRN_QLOSS.h>

Func* FUNC_BRN_QLOSS_new(REAL weight, Net* net) {
    Func* f = FUNC_new(weight, net);
    FUNC_set_func_count_step(f, &FUNC_BRN_QLOSS_count_step);
    FUNC_set_func_analyze_step(f, &FUNC_BRN_QLOSS_analyze_step);
    FUNC_set_func_eval_step(f, &FUNC_BRN_QLOSS_eval_step);
    FUNC_set_name(f, "branch Q loss");
    return f;
}

void FUNC_BRN_QLOSS_count_step(Func* f, Bus* bus, BusDC* busdc, int t) {

    // Local variables
    int* Hphi_nnz;
    Branch* br;
    Bus* bus_m;

    // Func data
    Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

    // Check pointers
    if (!Hphi_nnz || !bus)
        return;

    // Out of service or bus not in subsys
    if (!BUS_is_in_service(bus) || !BUS_is_in_subsys(bus))
        return;
    
    // Branches_k
    for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

        // Out of service
        if (!BRANCH_is_in_service(br))
            continue;

        bus_m = BRANCH_get_bus_m(br);

        // Bus_m not in subsys or not in service
        if (!BUS_is_in_subsys(bus_m) || !BUS_is_in_service(bus_m))
            continue;

        // Hphi
        // vk is variable
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {
            
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {

                // vk vm
                (*Hphi_nnz)++;
            }
            
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {

                // vk wm
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

                // vk a_km
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                // vk phsh
                (*Hphi_nnz)++;
            }
        }

        // wk
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {
            
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {

                // wk vm
                (*Hphi_nnz)++;
            }
            
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {

                // wk wm
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

                // wk a_km
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                // wk phsh
                (*Hphi_nnz)++;
            }
        }

        // vm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {

            // vm vm
            (*Hphi_nnz)++;
            
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {

                // vm wm
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

                // vm a_km
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                // vm phsh
                (*Hphi_nnz)++;
            }
        }

        // wm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {

            // wm wm
            (*Hphi_nnz)++;
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

                // wm a_km
                (*Hphi_nnz)++;
            }
            
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                // wm phsh
                (*Hphi_nnz)++;
            }
        }

        // a_km 
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

            // a_km a_km
            (*Hphi_nnz)++;
           
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                // a_km phsh
                (*Hphi_nnz)++;
            }
        }

        // phsh
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

            // phsh phsh
            (*Hphi_nnz)++;
        }
    }

    // vk is variable
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {

        // vk vk
        (*Hphi_nnz)++;
        
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {

            // vk wk
            (*Hphi_nnz)++;
        }
    }

    // wk
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {

        // wk wk
        (*Hphi_nnz)++;
    }

}

void FUNC_BRN_QLOSS_analyze_step(Func* f, Bus* bus, BusDC* busdc, int t) {

    // Local variables
    Branch* br;
    Bus* bus_m;
    int* Hphi_nnz;
    Mat* Hphi;

    // Constr data
    Hphi = FUNC_get_Hphi(f);
    Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

    // Check pointers
    if (!Hphi_nnz || !Hphi || !bus)
        return;

    // Out of service or bus not in subsys
    if (!BUS_is_in_service(bus) || !BUS_is_in_subsys(bus))
        return;

    // Branches_k
    for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

        // Out of service
        if (!BRANCH_is_in_service(br))
            continue;

        bus_m = BRANCH_get_bus_m(br);

        // Bus_m not in subsys or not in service
        if (!BUS_is_in_subsys(bus_m) || !BUS_is_in_service(bus_m))
            continue;

        // Hphi
        // vk is variable
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {

            // vk vm
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
                (*Hphi_nnz)++;
            }

            // vk wm
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
                (*Hphi_nnz)++;
            }

            // vk a_km
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
                (*Hphi_nnz)++;
            }

            // vk phsh
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
                (*Hphi_nnz)++;
            }
        }

        // wk
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {

            // wk vm
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
                (*Hphi_nnz)++;
            }

            // wk wm
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
                (*Hphi_nnz)++;
            }

            // wk a_km
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
                (*Hphi_nnz)++;
            }

            // wk phsh
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
                (*Hphi_nnz)++;
            }
        }

        // vm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {

            // vm vm
            MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
            MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
            (*Hphi_nnz)++;   

            // vm wm
            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { 
                
                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
                MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
                (*Hphi_nnz)++;
            }

            // vm a_km
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
                (*Hphi_nnz)++;
            }

            // vm phsh
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus_m, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
                (*Hphi_nnz)++;
            }
        }

        // wm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {
            
            // wm wm
            MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
            MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
            (*Hphi_nnz)++;   

            // wm a_km
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { 

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
                (*Hphi_nnz)++;
            }

            // wm phsh
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

                MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus_m, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
                (*Hphi_nnz)++;
            }
        }

        // a_km 
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

            // a_km a_km
            MAT_set_i(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
            MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
            (*Hphi_nnz)++;   

            // a_km phsh
            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {
                
                MAT_set_i(Hphi, *Hphi_nnz, BRANCH_get_index_ratio(br, t));
                MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
                (*Hphi_nnz)++;
            }
        }

        // phsh
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {
            
            // phsh phsh
            MAT_set_i(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
            MAT_set_j(Hphi, *Hphi_nnz, BRANCH_get_index_phase(br, t));
            (*Hphi_nnz)++;    
        }
    }

    // vk is variable
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {

        // vk vk
        MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
        MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
        (*Hphi_nnz)++;  

        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {
            
            // vk wk
            MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_mag(bus, t));
            MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
            (*Hphi_nnz)++;
        }
    }

    // wk
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {

        // wk wk
        MAT_set_i(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
        MAT_set_j(Hphi, *Hphi_nnz, BUS_get_index_v_ang(bus, t));
        (*Hphi_nnz)++;  

    }
}

void FUNC_BRN_QLOSS_eval_step(Func* f, Bus* bus, BusDC* busdc, int t, Vec* var_values) {

    // Local variables
    Branch* br;
    REAL* phi;
    REAL* gphi;
    REAL* Hphi;
    Bus* bus_m;
    int* Hphi_nnz;

    REAL vk;
    REAL vm;
    REAL wk;
    REAL wm;
    REAL theta;   // bransh admittance angle
    REAL y;       // bransh admittance magnitude |brn.g + j brn.b|
    REAL phsh;    // phase shift angle
    REAL a_km;    // tap ratio from k
    REAL b;       // branch susceptance

    REAL sin_w_k_m_ph_th;  // sin(wk - wm - phsh - theta)
    REAL cos_w_k_m_ph_th;  // cos(wk - wm - phsh - theta)
    REAL sin_w_m_k_ph_th;  // sin(wm - wk + phsh - theta)
    REAL cos_w_m_k_ph_th;  // cos(wm - wk + phsh - theta)

    REAL hphi_vk = 0;
    REAL hphi_wk = 0;
    REAL hphi_vkwk = 0;

    int index_vk;
    int index_wk;
    int index_vm;
    int index_wm;
    int index_a_km;
    int index_phsh;

    // Func data
    phi = FUNC_get_phi_ptr(f);
    gphi = VEC_get_data(FUNC_get_gphi(f));
    Hphi = MAT_get_data_array(FUNC_get_Hphi(f));
    Hphi_nnz = FUNC_get_Hphi_nnz_ptr(f);

    // Check pointers
    if (!phi || !gphi || !Hphi || !Hphi_nnz || !bus)
        return;

    // Out of service or bus not in subsys
    if (!BUS_is_in_service(bus) || !BUS_is_in_subsys(bus))
        return;

    // Bus voltage mag
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {

        index_vk = BUS_get_index_v_mag(bus, t);
        vk = VEC_get(var_values, index_vk);
    }
    else {
        vk = BUS_get_v_mag(bus, t);
    }

    // Bus voltage angle
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {

        index_wk = BUS_get_index_v_ang(bus, t);
        wk = VEC_get(var_values, index_wk);
    }
    else {
        wk = BUS_get_v_ang(bus, t);
    }

    // Branches_k
    for (br = BUS_get_branch_k(bus); br != NULL; br = BRANCH_get_next_k(br)) {

        // Out of service
        if (!BRANCH_is_in_service(br))
            continue;

        bus_m = BRANCH_get_bus_m(br);

        // Bus_m not in subsys or not in service
        if (!BUS_is_in_subsys(bus_m) || !BUS_is_in_service(bus_m))
            continue;

        // Bus_m voltage mag
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {

            index_vm = BUS_get_index_v_mag(bus_m, t);
            vm = VEC_get(var_values, index_vm);
        }
        else {
            vm = BUS_get_v_mag(bus_m, t);
        }

        // Bus_m voltage angle
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {

            index_wm = BUS_get_index_v_ang(bus_m, t);
            wm = VEC_get(var_values, index_wm);
        }
        else {
            wm = BUS_get_v_ang(bus_m, t);
        }

        // Branch Data
        y = sqrt(pow(BRANCH_get_b(br), 2.) + pow(BRANCH_get_g(br), 2.));
        // atan2(y/x)
        theta = atan2(BRANCH_get_b(br), BRANCH_get_g(br));
        b = BRANCH_get_b(br);

        // Get phase shift angle, this should be subtracted on side "m" (to)
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {

            index_phsh = BRANCH_get_index_phase(br, t);
            phsh = VEC_get(var_values, index_phsh);
        }
        else {
            phsh = BRANCH_get_phase(br, t);
        }

        // Get tap ratio from k, a_mk = 1 always
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {

            index_a_km = BRANCH_get_index_ratio(br, t);
            a_km = VEC_get(var_values, index_a_km);
        }
        else {
            a_km = BRANCH_get_ratio(br, t);
        }

        // Calculate handy values
        sin_w_k_m_ph_th = sin(wk - wm - phsh - theta);
        cos_w_k_m_ph_th = cos(wk - wm - phsh - theta);
        sin_w_m_k_ph_th = sin(wm - wk + phsh - theta);
        cos_w_m_k_ph_th = cos(wm - wk + phsh - theta);

        // phi
        // Equation for reference:
        // I^2 x = - a_km^2*vk^2*b - vm^2*b - a_km*vm*vk*y(sin(wk-wm-phsh-theta_brn) + sin(wm-wk+phsh-theta_brn))
        // 
        (*phi) += - pow(a_km,2)*pow(vk,2)*b - pow(vm,2)*b - a_km*vk*vm*y*(sin_w_k_m_ph_th + sin_w_m_k_ph_th);

        // gphi and Hphi
        //          vk      wk      vm      wm     a_km    phsh 
        // vk       * 
        // wk       *       *       
        // vm       *       *       *              
        // wm       *       *       *       *       
        // a_km     *       *       *       *       *       
        // phsh     *       *       *       *       *       *       

        // vk
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {
            gphi[index_vk] += -2 * pow(a_km, 2) * vk * b - a_km * vm * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
            hphi_vk += -2 * pow(a_km, 2) * b;  // vk&vk

            if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) { // vk wk
                hphi_vkwk += -a_km * vm * y * (cos_w_k_m_ph_th - cos_w_m_k_ph_th);
            }

            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) { // vk vm
                Hphi[*Hphi_nnz] += -a_km * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { // vk wm
                Hphi[*Hphi_nnz] += -a_km * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { // vk a_km
                Hphi[*Hphi_nnz] += -4 * a_km * vk * b - vm * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) { // vk phsh
                Hphi[*Hphi_nnz] += -a_km * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }
        }

        // wk
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {
            gphi[index_wk] += -a_km * vk * vm * y * (cos_w_k_m_ph_th - cos_w_m_k_ph_th);
            hphi_wk += -a_km * vk * vm * y * (-sin_w_k_m_ph_th - sin_w_m_k_ph_th);  // wk&wk


            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) { // wk vm
                Hphi[*Hphi_nnz] += -a_km * vk * y * (cos_w_k_m_ph_th - cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { // wk wm
                Hphi[*Hphi_nnz] += -a_km * vk * vm * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { // wk a_km
                Hphi[*Hphi_nnz] += -vk * vm * y * (cos_w_k_m_ph_th - cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) { // wk phsh
                Hphi[*Hphi_nnz] += -a_km * vk * vm * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }
        }

        // vm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VMAG)) {
            gphi[index_vm] += -2 * vm * b - a_km * vk * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
            Hphi[*Hphi_nnz] += -2 * b;  // vm&vm
            (*Hphi_nnz)++;

            if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) { // vm wm
                Hphi[*Hphi_nnz] += -a_km * vk * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { // vm a_km
                Hphi[*Hphi_nnz] += -vk * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) { // vm phsh
                Hphi[*Hphi_nnz] += -a_km * vk * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }
        }

        // wm
        if (BUS_has_flags(bus_m, FLAG_VARS, BUS_VAR_VANG)) {
            gphi[index_wm] += -a_km * vk * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
            Hphi[*Hphi_nnz] += -a_km * vk * vm * y * (-sin_w_k_m_ph_th - sin_w_m_k_ph_th);  // wm&wm
            (*Hphi_nnz)++;

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) { // wm a_km
                Hphi[*Hphi_nnz] += -vk * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) { // wm phsh
                Hphi[*Hphi_nnz] += -a_km * vk * vm * y * (-sin_w_k_m_ph_th - sin_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }
        }

        // a_km 
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_RATIO)) {
            gphi[index_a_km] += -2 * a_km * pow(vk, 2.) * b - vk * vm * y * (sin_w_k_m_ph_th + sin_w_m_k_ph_th);
            Hphi[*Hphi_nnz] += -2 * pow(vk, 2.) * b;  // a_km a_km
            (*Hphi_nnz)++;

            if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) { // a_km phsh
                Hphi[*Hphi_nnz] += -vk * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
                (*Hphi_nnz)++;
            }
        }

        // phsh
        if (BRANCH_has_flags(br, FLAG_VARS, BRANCH_VAR_PHASE)) {
            gphi[index_phsh] += -a_km * vk * vm * y * (-cos_w_k_m_ph_th + cos_w_m_k_ph_th);
            Hphi[*Hphi_nnz] += -a_km * vk * vm * y * (-sin_w_k_m_ph_th - sin_w_m_k_ph_th);   // phsh phsh
            (*Hphi_nnz)++;
        }
    }

    // vk vk 
    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VMAG)) {
        Hphi[*Hphi_nnz] += hphi_vk;  // vk&vk
        (*Hphi_nnz)++;
        
        if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) { // vk wk
            Hphi[*Hphi_nnz] += hphi_vkwk;
            (*Hphi_nnz)++;
        }
    }

    if (BUS_has_flags(bus, FLAG_VARS, BUS_VAR_VANG)) {
        Hphi[*Hphi_nnz] += hphi_wk;  // wk&wk
        (*Hphi_nnz)++;
    }
}
