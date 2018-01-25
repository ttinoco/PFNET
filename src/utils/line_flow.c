/* Copyright (c) 2017, ETH Zurich, Dmitry Shchetinin
 * https://github.com/dmitry-shchetinin/LineFlow
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions 
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright 
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright 
 * notice, this list of conditions and the following disclaimer in the 
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its 
 * contributors may be used to endorse or promote products derived from 
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "line_flow.h"

/*define custom types*/
struct LF_Options { //contains parameters of the algorithm
	int computation_mode;
	int iter_Max;
	int N_constraints_max;
	int N_adjustments;
	int tr_model_type;
	double delta_max_user;
	double eps_tolerance;
	double error_max;
	double max_error_change;
	double ratio_threshold;
	ApproximationType approximation;
};


struct LF_Results { //contains results of the linearization
	double* A_matrix;
	double* b_vector;
	LF_ResultFlag flag_result;
	int N_created_constraints;
	double max_error_est;
	char message[80];
};


struct Point3D { //point in a 3D space
	double V_i;
	double V_j;
	double delta;
};


struct IntersectionLine { //contains information related to intersection lines
	Point3D point_begin;
	Point3D point_end;
	double slope; //in (V_i, delta) plane
	double tangent_slope_begin; //ddelta/dV_i_begin
	double tangent_slope_end; //ddelta/dV_i_begin
	double offset;
};


struct FeasibleRegion { //contains information related to feasible region
	BOOL top_left_corner_feasible; //shows whether point (V_i_min, V_j_max) is feasible
	BOOL bottom_right_corner_feasible; //shows whether point (V_i_max, V_j_min) is feasible
	BOOL convex; //shows whether the region is convex
	BOOL approximation_needed; //shows whether linear approximation is needed
	double I_at_Vbox_corners[4]; //values of I at corners of (V_i, V_j) box and delta=delta_max_user
	double offset_min; //max(offset of lower boundary line, offset of || to ridge line that passes through bottom right corner)
	double offset_max; //min(offset of upper boundary line, offset of || to ridge line that passes through top left corner)
};


struct VoltageLimits { //describes min and max values of voltage magnitudes at the beginning (V_i) and end (V_j) of the line
	double V_i_min;
	double V_i_max;
	double V_j_min;
	double V_j_max;
};


struct LF_Workspace { //contains intermediate variables used to construct approximation
	IntersectionLine* lines;
	double* plane_normals;
	double* b_planes;
	double* errors;
	VoltageLimits V_limits;
	FeasibleRegion feas_reg;
	double max_error_estimate;
	double branch_data[10];
	double I_max_user; //after transformer
	double delta_max;
	double slope_all_lines;
	int constraint_counter;
	BranchType branch_type;
	BOOL all_points_reachable;
};


struct LineSegmentEdge { //contains parameters needed to establish convexity of feasible region for line between gen and load
	TangentLineType tangent_line_type;
	double V_fixed;
	double V_min;
	double V_max;
	double delta_V_min;
	double delta_V_max;
	double d_delta_V_min;
	double d_delta_V_max;
	double line_slope;
	double line_offset;
	BOOL convex;
};




/*Outer interface function. Check input branch branch_data, initialize workspace and output structures,
and construct linear approximation*/
LF_Results* LF_construct(LF_Branch* branch, int flow_side, LF_Options* options_in) {
	LF_Options* options;
	LF_Results* results;
	LF_Workspace* workspace;
	BOOL use_default_options = LF_FALSE;

	if (!options_in) {
		options = LF_get_default_options();
		use_default_options = LF_TRUE;
	}
	else
		options = options_in;

	//initialize output structure
	results = LF_initialize_results(options, flow_side);

	//check algorithm's options
	if (!LF_check_options(options, results, flow_side))
		return results; //if options are incorrect, linearization is not performed, error message is recorded

	//check branch branch_data
	if (!LF_check_branch_data(branch, results))
		return results; //if branch branch_data are incorrect, linearization is not performed, error message is recorded

	//initialize LF_Workspace variables
	workspace = LF_initizalize_workspace(options, flow_side);

	//CONSTRUCT LINEAR APPROXIMATION
	LF_linearize_one_line(branch, flow_side, options, workspace, results);

	//free memory
	LF_free_workspace(workspace);
	if (use_default_options) {
		free(options);
		options = NULL;
	}

	return results;
}


/*create default algorithm options*/
LF_Options* LF_get_default_options() {
	LF_Options* options;
	options = (LF_Options *)malloc(sizeof(LF_Options));
	options->computation_mode = 2;
	options->iter_Max = 25;
	options->N_constraints_max = 15;
	options->N_adjustments = 4;
	options->tr_model_type = 0;
	options->delta_max_user = 85.0 * LF_PI / 180.0;
	options->eps_tolerance = 0.0001;
	options->error_max = 5.0;
	options->max_error_change = 0.1;
	options->ratio_threshold = 0.9;
	options->approximation = conservative;
	return options;
}


/*initialize output structure*/
LF_Results* LF_initialize_results(LF_Options* options, int flow_side) {
	LF_Results* results = (LF_Results *)malloc(sizeof(LF_Results));
	int coef;
	if (flow_side == 3)
		coef = 2;
	else
		coef = 1;

	if (options->N_constraints_max > 0) {
		results->A_matrix = (double *)malloc(options->N_constraints_max * 6 * coef * sizeof(double));
		results->b_vector = (double *)malloc(options->N_constraints_max * 2 * coef * sizeof(double));
	}
	else {
		results->A_matrix = NULL;
		results->b_vector = NULL;
	}
	results->flag_result = non_binding;
	results->max_error_est = 0.0;
	results->N_created_constraints = 0;
	return results;
}


/*free memory taken by output structure*/
void LF_free_results(LF_Results* results) {
	free(results->A_matrix);
	free(results->b_vector);
	free(results);
	results = NULL;
}


/*check correctness of algorithm options*/
BOOL LF_check_options(LF_Options* options, LF_Results* results, int flow_side) {
	//check computation mode
	if (options->computation_mode < 1 || options->computation_mode>2) {
		strcpy(results->message, "Input error: computation mode can be either 1 or 2.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check number of iterations
	if (options->iter_Max < 25) {
		strcpy(results->message, "Input error: must have iter_Max>=25.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check number of adjustments
	if (options->N_adjustments < 0) {
		strcpy(results->message, "Input error: must have N_adjustments>=0.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check transformer model type
	if (options->tr_model_type < 0 || options->tr_model_type>1) {
		strcpy(results->message, "Input error: transformer model type can be 0 or 1.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check maximum angle 
	if (options->delta_max_user<LF_PI / 4.0 || options->delta_max_user>88.0 * LF_PI / 180.0) {
		strcpy(results->message, "Input error: must have 45<=delta_max_user<=88 degrees.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check tolerance
	if (options->eps_tolerance <= 0.0) {
		strcpy(results->message, "Input error: eps_tolerance must have a positive value.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check max error change
	if (options->max_error_change <= 0.0) {
		strcpy(results->message, "Input error: max_error_change must have a positive value.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check ratio threshold
	if (options->ratio_threshold <= 0.0 || options->ratio_threshold >= 1.0) {
		strcpy(results->message, "Input error: must have 0<ratio_threshold<1.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check that max number of constraints is positive
	if (options->N_constraints_max <= 0) {
		strcpy(results->message, "Input error: must have N_constraints_max>=1.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check that maximum error is positive
	if (options->error_max <= 0.0) {
		strcpy(results->message, "Input error: maximum error must have a positive value.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check that flow side is 1, 2 or 3
	if (flow_side < 1 || flow_side>3) {
		strcpy(results->message, "Input error: flow_side can be 1, 2, or 3.");
		results->flag_result = error_options;
		return LF_FALSE;
	}

	//check that approximation is either conservative or relaxed
	if (options->approximation != conservative && options->approximation != relaxed) {
		strcpy(results->message, "Input error: approximation must be either conservative or relaxed.");
		results->flag_result = error_options;
		return LF_FALSE;
	}
	return LF_TRUE;
}


/*check correctness of line branch_data*/
BOOL LF_check_branch_data(LF_Branch* branch, LF_Results* results) {
	//check lower bounds on voltage magnitudes
	if (branch->V_i_min <= 0.0 || branch->V_j_min <= 0.0) {
		strcpy(results->message, "Input error: lower bounds on voltage magnitudes must be positive.");
		results->flag_result = error_branch_data;
		return LF_FALSE;
	}

	//check consistency of bounds on voltage magnitudes
	if (branch->V_i_min > branch->V_i_max || branch->V_j_min > branch->V_j_max) {
		strcpy(results->message, "Input error: must have V_min<=V_max.");
		results->flag_result = error_branch_data;
		return LF_FALSE;
	}

	//check that branch has non-zero series parameters
	if ((branch->g == 0.0 && branch->b == 0.0) || fabs(branch->g) > 1e6 || fabs(branch->b) > 1e6) {
		strcpy(results->message, "Input error: g and b cannot be >1e6 or simultaneously zero.");
		results->flag_result = error_branch_data;
		return LF_FALSE;
	}

	//check that transformer's ratio is nonnegative
	if (branch->t_ratio < 0.0) {
		strcpy(results->message, "Input error: transformer's ratio must be nonnegative.");
		results->flag_result = error_branch_data;
		return LF_FALSE;
	}

	//check that thermal limit is positive
	if (branch->I_max < 0.0) {
		strcpy(results->message, "Input error: value of maximum current must be positive.");
		results->flag_result = error_branch_data;
		return LF_FALSE;
	}
	else if (branch->I_max == 0.0) {
		strcpy(results->message, "Thermal limit is zero => no approximation was constructed.");
		results->flag_result = zero_limit;
		return LF_FALSE;
	}
	return LF_TRUE;
}


/*allocate memory for workspace*/
LF_Workspace* LF_initizalize_workspace(LF_Options* options, int flow_side) {
	int num_lines, coef;
	LF_Workspace* workspace = (LF_Workspace *)malloc(sizeof(LF_Workspace));

	if (flow_side == 3)
		coef = 2;
	else
		coef = 1;

	//get maximum number of intersection lines and linear constraints on one side of the surface
	num_lines = max_number_of_intersection_lines(options->N_constraints_max);

	//allocate memory for arrays and structures
	workspace->lines = (IntersectionLine*)malloc(num_lines*sizeof(IntersectionLine));
	workspace->plane_normals = (double *)malloc(6 * coef * options->N_constraints_max * sizeof(double));
	workspace->b_planes = (double *)malloc(2 * coef * options->N_constraints_max * sizeof(double));
	workspace->errors = (double *)malloc((options->N_constraints_max + 1) * sizeof(double));
	return workspace;
}


/*free memory taken by workspace*/
void LF_free_workspace(LF_Workspace* workspace) {
	free(workspace->lines);
	free(workspace->plane_normals);
	free(workspace->b_planes);
	free(workspace->errors);
	free(workspace);
	workspace = NULL;
}


/*construct linear approximation using provided workspace and results structures*/
void LF_linearize_one_line(LF_Branch* branch, int flow_side, LF_Options* options, LF_Workspace* workspace, LF_Results* results) {
	int N_constraints_end = 0;
	workspace->constraint_counter = 0;
	workspace->max_error_estimate = 0;
	if (flow_side == 3 && branch->b_sh == 0.0) {
		if (branch->t_ratio <= 1.0)
			flow_side = 1;
		else
			flow_side = 2;
	}

	if (flow_side != 3) {
		//compute values derived from branch parameters that are necessary to construct approximation
		compute_branch_parameters(branch, flow_side, options, workspace);

		//construct approximation for the desired end of the line
		linearize_one_line_end(flow_side, options, workspace, results);
		if (results->flag_result == infeasible)
			return;

		//record output parameters
		record_all_constraints(workspace, branch, results);
	}
	else {
		//construct approximation for the end of the line
		compute_branch_parameters(branch, 2, options, workspace);
		linearize_one_line_end(2, options, workspace, results);
		if (results->flag_result == infeasible)
			return;
		N_constraints_end = workspace->constraint_counter;

		//construct approximation for the beginning of the line
		compute_branch_parameters(branch, 1, options, workspace);
		linearize_one_line_end(1, options, workspace, results);
		if (results->flag_result == infeasible)
			return;

		//record relevant output parameters
		if (N_constraints_end >= 2 && workspace->constraint_counter >= N_constraints_end + 2 &&
			(branch->t_ratio == 0.0 || (branch->t_ratio == 1.0 && options->tr_model_type == 0)))
			record_relevant_constraints(workspace, branch, N_constraints_end, results);
		else
			record_all_constraints(workspace, branch, results);
	}
}


/*construct linear approximation for one end of the line*/
void linearize_one_line_end(int flow_side, LF_Options* options, LF_Workspace* workspace, LF_Results* results) {
	VoltageLimits V_limits_initial;
	int N_constraints_for_upper_part = workspace->constraint_counter;
	
	//check feasibility of the nonlinear constraint
	if (is_feasible_region_empty(workspace)) {
		strcpy(results->message, "Nonlinear constraint is infeasible.");
		results->flag_result = infeasible;
		results->N_created_constraints = 0;
		return;
	}

	V_limits_initial = workspace->V_limits;

	//approximate upper part of the surface
	workspace->branch_data[9] = 1.0;
	workspace->delta_max = options->delta_max_user;
	compute_parameters_of_feasible_region(workspace, options);
	if (workspace->feas_reg.approximation_needed)
		compute_parameters_of_approximation(workspace, options);

	//approximate lower part of the surface
	if (is_approximation_symmetric(workspace)) {
		N_constraints_for_upper_part = workspace->constraint_counter - N_constraints_for_upper_part;
		reflect_approximation(workspace, N_constraints_for_upper_part);
	}
	else {
		workspace->branch_data[9] = -1.0;
		workspace->delta_max = -options->delta_max_user;
		workspace->V_limits = V_limits_initial;
		compute_parameters_of_feasible_region(workspace, options);
		if (workspace->feas_reg.approximation_needed)
			compute_parameters_of_approximation(workspace, options);
	}
}





/*compute values derived from branch parameters that will be extensively used throughout the algorithm*/
void compute_branch_parameters(LF_Branch* branch, int flow_side, LF_Options* options, LF_Workspace* workspace) {
	double b_sh, ratio;

	//take transformer into account
	workspace->I_max_user = branch->I_max;
	if (branch->t_ratio == 0.0) {
		b_sh = branch->b_sh / 2.0;
		ratio = 1.0;
	}
	else { //compute parameters depending on the transformer's model
		if (options->tr_model_type == 0) //if Matpower convention is used
			b_sh = branch->b_sh / 2.0;
		else { //if Russian convention is used
			if (flow_side == 1)
				b_sh = branch->b_sh;
			else
				b_sh = 0.0;
		}
		ratio = branch->t_ratio;
		if (flow_side == 1)
			workspace->I_max_user = workspace->I_max_user*ratio;
	}

	//fill out voltage limits after the transformer
	workspace->V_limits.V_i_min = branch->V_i_min / ratio;
	workspace->V_limits.V_j_min = branch->V_j_min;
	workspace->V_limits.V_i_max = branch->V_i_max / ratio;
	workspace->V_limits.V_j_max = branch->V_j_max;

	//compute frequenly used values derived from g, b, b_sh
	if (flow_side == 1) {
		workspace->branch_data[0] = branch->g *branch->g + (branch->b + b_sh) *(branch->b + b_sh);
		workspace->branch_data[1] = branch->g * branch->g + branch->b * branch->b;
		workspace->branch_data[2] = branch->g * b_sh;
	}
	else {
		workspace->branch_data[0] = branch->g * branch->g + branch->b * branch->b;
		workspace->branch_data[1] = branch->g *branch->g + (branch->b + b_sh) *(branch->b + b_sh);
		workspace->branch_data[2] = -branch->g * b_sh;
	}
	workspace->branch_data[3] = branch->g * branch->g + branch->b * branch->b + branch->b*b_sh;
	workspace->branch_data[4] = workspace->I_max_user * workspace->I_max_user;
	workspace->branch_data[5] = 4.0 * workspace->branch_data[0] * workspace->branch_data[1];
	workspace->branch_data[6] = -workspace->branch_data[2] / (2.0 * workspace->branch_data[0] * workspace->branch_data[1]);
	workspace->branch_data[7] = fabs(workspace->branch_data[3]) / (2.0 * workspace->branch_data[0] * workspace->branch_data[1]);
	workspace->branch_data[8] = sqrt(workspace->branch_data[0] / workspace->branch_data[1]);

	//check what kind of branch we have
	if (branch->V_i_min == branch->V_i_max && branch->V_j_min == branch->V_j_max)
		workspace->branch_type = gen_gen;
	else if (branch->V_i_min != branch->V_i_max && branch->V_j_min != branch->V_j_max)
		workspace->branch_type = load_load;
	else
		workspace->branch_type = gen_load;
}


/*check if intersection of ViVj box with feasible region of line flow constraint is empty*/
BOOL is_feasible_region_empty(LF_Workspace* workspace) {
	//check if lower boundary line is above (V_i_min, V_j_max) point or upper boundary line is below (V_i_max, V_j_min) point
	if (workspace->V_limits.V_j_max<workspace->branch_data[8] * workspace->V_limits.V_i_min -
		sqrt(workspace->branch_data[4] / workspace->branch_data[1]) ||
		workspace->V_limits.V_j_min>workspace->branch_data[8] * workspace->V_limits.V_i_max + 
		sqrt(workspace->branch_data[4] / workspace->branch_data[1]))
		return LF_TRUE;
	else
		return LF_FALSE;
}


/*extract the information on the feasible region of line flow constraint*/
void compute_parameters_of_feasible_region(LF_Workspace* workspace, LF_Options* options) {

	//initialize values of flags that describe properties of feasible region to their default values
	set_default_flags_for_feasible_region(&workspace->feas_reg);

	//compute maximum possible values of I at the corners of the ViVj box
	compute_I_max_at_Vbox_corners(workspace);

	//check reachability of points in ViVj box
	if (max_element_in_array(workspace->feas_reg.I_at_Vbox_corners, 4) <= workspace->I_max_user) {
		workspace->feas_reg.approximation_needed = LF_FALSE; //all corners are not reachable
		return;
	}
	else if (min_element_in_array(workspace->feas_reg.I_at_Vbox_corners, 4) >= workspace->I_max_user)
		workspace->all_points_reachable = LF_TRUE; //all corners are reachable
	else {
		workspace->all_points_reachable = LF_FALSE; //there are reachable and unreachable corners
	    //tighten limits on V to make sure that at least three corners of ViVj box are reachable
		tighten_V_limits_for_reachability(workspace);
	}

	//check and record whether the whole ViVj box is feasible
	establish_feasibility_of_Vbox_corners(workspace);

	//compute min amd max offsets of a line that is || to the ridge and contains at least one feasible point
	compute_offsets_of_outermost_lines(workspace);

	//perform extra steps if the line is between gen and load
	if (workspace->branch_type == gen_load) {
		//make sure all points within V limits are feasible
		tighten_V_limits_for_feasibility(workspace);
		//check and record if the region is convex and if it non-convex, retain only the part that causes nonconvexity
		establish_convexity_of_feasible_region(workspace, options);
	}

	//recompute maximum possible values of I at the corners of the ViVj box
	if (!workspace->all_points_reachable)
		compute_I_max_at_Vbox_corners(workspace);
}


/*initialize parameters of feasible region to their default values*/
void set_default_flags_for_feasible_region(FeasibleRegion* feas_reg) {
	feas_reg->top_left_corner_feasible = LF_TRUE;
	feas_reg->bottom_right_corner_feasible = LF_TRUE;
	feas_reg->convex = LF_TRUE;
	feas_reg->approximation_needed = LF_TRUE;
}


/*compute values of current at box corners and delta=+-delta_max_user-Kt_shift*/
void compute_I_max_at_Vbox_corners(LF_Workspace* workspace) {
	workspace->feas_reg.I_at_Vbox_corners[0] = current_magnitude(workspace->V_limits.V_i_min, workspace->V_limits.V_j_min,
		workspace->delta_max, workspace->branch_data);
	workspace->feas_reg.I_at_Vbox_corners[1] = current_magnitude(workspace->V_limits.V_i_min, workspace->V_limits.V_j_max,
		workspace->delta_max, workspace->branch_data);
	workspace->feas_reg.I_at_Vbox_corners[2] = current_magnitude(workspace->V_limits.V_i_max, workspace->V_limits.V_j_min,
		workspace->delta_max, workspace->branch_data);
	workspace->feas_reg.I_at_Vbox_corners[3] = current_magnitude(workspace->V_limits.V_i_max, workspace->V_limits.V_j_max,
		workspace->delta_max, workspace->branch_data);
}


/*tighten voltage limits to make sure that at least three corners of the box can have I>=I_max_user for delta<=delta_max_user*/
void tighten_V_limits_for_reachability(LF_Workspace* workspace) {
	if (workspace->feas_reg.I_at_Vbox_corners[3] >= workspace->I_max_user) {
		//if (V_i_max, V_j_max) corner is reachable, we check two adjacent corners
		if (workspace->feas_reg.I_at_Vbox_corners[1] < workspace->I_max_user) {
			//then(V_i_min, V_j_max) corner is unreachable => shrinking the box along V_i dimension
			workspace->V_limits.V_i_min = tightened_V_limit(workspace, workspace->V_limits.V_j_max, 2);
		}
		if (workspace->feas_reg.I_at_Vbox_corners[2] < workspace->I_max_user) {
			//then(V_i_max, V_j_min) corner is unreachable => shrinking the box along V_j dimension
			workspace->V_limits.V_j_min = tightened_V_limit(workspace, workspace->V_limits.V_i_max, 1);
		}
	}
	else { //This is a rare case that happens for high values of b_sh. Check which corner is reachable
		if (workspace->feas_reg.I_at_Vbox_corners[1] >= workspace->I_max_user) { //then(V_i_min, V_j_max) is the only reachable corner
			//shrink the box along V_i_dimension
			workspace->V_limits.V_i_max = tightened_V_limit(workspace, workspace->V_limits.V_j_max, 2);
			//shrink the box along V_j dimension
			workspace->V_limits.V_j_min = tightened_V_limit(workspace, workspace->V_limits.V_i_min, 1);
		}
		else { //then(V_i_max, V_j_min) is the only reachable corner
			//shrink the box along V_i_dimension
			workspace->V_limits.V_i_min = tightened_V_limit(workspace, workspace->V_limits.V_j_min, 2);
			//shrink the box along V_j dimension
			workspace->V_limits.V_j_max = tightened_V_limit(workspace, workspace->V_limits.V_i_max, 1);
		}
	}
}


/*tighten limit on V_i or V_j to make sure that at this limit I=I_max_user for delta=delta_max_user*/
double tightened_V_limit(LF_Workspace* workspace, double V_value, int flag_V_fixed) {
	double V_value_opt = 0.0, a6, D, V_root1, V_root2, buf1, buf2, V_min, V_max;
	//UPDATE MATH DESCRIPTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (flag_V_fixed == 1) {
		//V_i is fixed => shrink the box along V_j dimension
		buf1 = workspace->branch_data[1];
		buf2 = workspace->branch_data[0];
		V_min = workspace->V_limits.V_j_min;
		V_max = workspace->V_limits.V_j_max;
	}
	else {
		//V_j is fixed => shrink the box along V_i dimension
		buf1 = workspace->branch_data[0];
		buf2 = workspace->branch_data[1];
		V_min = workspace->V_limits.V_i_min;
		V_max = workspace->V_limits.V_i_max;
	}

	if (V_min == V_max)
		return V_min;

	a6 = -workspace->branch_data[3] * cos(workspace->delta_max) + workspace->branch_data[2] * sin(workspace->delta_max);
	D = (V_value*V_value*a6*a6) - buf1 *
		(buf2 * V_value*V_value - workspace->branch_data[4]);
	if (D < 0)
		D = 0;
	else
		D = sqrt(D);
	V_root1 = (-V_value*a6 - D) / buf1;
	V_root2 = (-V_value*a6 + D) / buf1;
	if (V_root1 > V_min && V_root1 < V_max)
		V_value_opt = V_root1;
	else if (V_root2 > V_min && V_root2 < V_max)
		V_value_opt = V_root2;
	return V_value_opt;
}


/*check and record whether top-left and bottom-right box corners are feasible*/
void establish_feasibility_of_Vbox_corners(LF_Workspace* workspace) {
	double offset = sqrt(workspace->branch_data[4] / workspace->branch_data[1]); //abs value of the offset of boundary lines
	//top_left_corner
	if (workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_min > offset)
		workspace->feas_reg.top_left_corner_feasible = LF_FALSE;
	//bottom-right corner
	if (fabs(workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_max) > offset)
		workspace->feas_reg.bottom_right_corner_feasible = LF_FALSE;
}


/*compute min amd max offsets of a line that is || to the ridge and contains at least one feasible point*/
void compute_offsets_of_outermost_lines(LF_Workspace* workspace) {
	double offset = sqrt(workspace->branch_data[4] / workspace->branch_data[1]); //abs value of the offset of boundary lines

	//top_left_corner
	if (workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_min > offset)
		workspace->feas_reg.offset_max = offset;
	else {
		if (workspace->feas_reg.I_at_Vbox_corners[3] >= workspace->I_max_user ||
			workspace->feas_reg.I_at_Vbox_corners[1] >= workspace->I_max_user)
			workspace->feas_reg.offset_max = workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_min;
		else
			workspace->feas_reg.offset_max = fmax(workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_max,
				workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_min);
	}

	//bottom-right corner
	if (fabs(workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_max) > offset)
		workspace->feas_reg.offset_min = -offset;
	else {
		if (workspace->feas_reg.I_at_Vbox_corners[3] >= workspace->I_max_user ||
			workspace->feas_reg.I_at_Vbox_corners[2] >= workspace->I_max_user)
			workspace->feas_reg.offset_min = workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_max;
		else
			workspace->feas_reg.offset_min = fmin(workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_max,
				workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_min);
	}
}


/*tighten limits on V for line between gen and load to make sure all V points within limits are feasible*/
void tighten_V_limits_for_feasibility(LF_Workspace* workspace) {
	double  offset1, offset2, coef_ViVj;
	offset1 = sqrt(workspace->branch_data[4] / workspace->branch_data[1]);
	offset2 = sqrt(workspace->branch_data[4] / workspace->branch_data[0]);
	coef_ViVj = workspace->branch_data[8];

	if (workspace->V_limits.V_i_min == workspace->V_limits.V_i_max) {
		coef_ViVj = workspace->branch_data[8];
		if (workspace->V_limits.V_j_max >= coef_ViVj * workspace->V_limits.V_i_min + offset1)
			workspace->V_limits.V_j_max = coef_ViVj * workspace->V_limits.V_i_min + offset1;
		if (workspace->V_limits.V_j_min <= coef_ViVj * workspace->V_limits.V_i_min - offset1)
			workspace->V_limits.V_j_min = coef_ViVj * workspace->V_limits.V_i_min - offset1;
	}
	else {
		coef_ViVj = 1.0 / workspace->branch_data[8];
		if (workspace->V_limits.V_i_max >= coef_ViVj*workspace->V_limits.V_j_min + offset2)
			workspace->V_limits.V_i_max = coef_ViVj*workspace->V_limits.V_j_min + offset2;
		if (workspace->V_limits.V_i_min <= coef_ViVj*workspace->V_limits.V_j_min - offset2)
			workspace->V_limits.V_i_min = coef_ViVj*workspace->V_limits.V_j_min - offset2;
	}
}


/*check if the feasible region is convex for line between gen and load and if it is non-convex, retain only the part that causes nonconvexity*/
void establish_convexity_of_feasible_region(LF_Workspace* workspace, LF_Options* options) {
	LineSegmentEdge segment = { .convex = LF_FALSE };
	double V_inflection, sign=workspace->branch_data[9];

	if (workspace->feas_reg.bottom_right_corner_feasible == LF_FALSE ||
		workspace->feas_reg.top_left_corner_feasible == LF_FALSE) {
		workspace->feas_reg.convex = LF_TRUE;  //nonconvexity can only happen if both box corners are feasible
		return;
	}

	compute_parameters_of_edge_line_segment(workspace->branch_data, &workspace->V_limits, &segment);

	if (sign*segment.line_slope > sign*segment.d_delta_V_min || sign*segment.line_slope < sign*segment.d_delta_V_max)
		workspace->feas_reg.convex = LF_FALSE;
	else {
		workspace->feas_reg.convex = LF_TRUE;
		return;
	}

	//in case of non-convex segment, retain only the part of the interval on which the curve creates non-convexity of the segment
	if ((segment.line_slope<segment.d_delta_V_min && segment.line_slope<segment.d_delta_V_max) ||
		(segment.line_slope>segment.d_delta_V_min && segment.line_slope>segment.d_delta_V_max)) { //curve changes its curvature on the interval
		//find inflection point (the point at which the curve changes its curvature)
		V_inflection = compute_value_by_bisection(segment.V_min, segment.V_max, options->eps_tolerance, workspace->branch_data,
			bisection_check_for_inflection_point, &segment);

		//tighten bounds on V (retain only the part of interval on which the curve creates non-convexity of the segment)
		if (sign*segment.line_slope>sign*segment.d_delta_V_min) {
			if (segment.tangent_line_type == fixed_V_i)
				workspace->V_limits.V_j_max = V_inflection;
			else
				workspace->V_limits.V_i_max = V_inflection;
		}
		else {
			if (segment.tangent_line_type == fixed_V_i)
				workspace->V_limits.V_j_min = V_inflection;
			else
				workspace->V_limits.V_i_min = V_inflection;
		}
	}
}


/*compute specific parameters of edge line segment for line between load and generator*/
void compute_parameters_of_edge_line_segment(double *branch_data, VoltageLimits* Vbox, LineSegmentEdge* segment) {
	if (fabs(Vbox->V_i_min - Vbox->V_i_max) < fabs(Vbox->V_j_min - Vbox->V_j_max)) {
		segment->tangent_line_type = fixed_V_i;
		segment->V_fixed = Vbox->V_i_min;
		segment->V_min = Vbox->V_j_min;
		segment->V_max = Vbox->V_j_max;
	}
	else {
		segment->tangent_line_type = fixed_V_j;
		segment->V_fixed = Vbox->V_j_min;
		segment->V_min = Vbox->V_i_min;
		segment->V_max = Vbox->V_i_max;
	}

	//compute surface angles at the ends of feasible interval
	segment->delta_V_min = max_angle_for_given_current(Vbox->V_i_min, Vbox->V_j_min, branch_data);
	segment->delta_V_max = max_angle_for_given_current(Vbox->V_i_max, Vbox->V_j_max, branch_data);

	//compute slope and offset of initial approximation line
	segment->line_slope = (segment->delta_V_max - segment->delta_V_min) / (segment->V_max - segment->V_min);
	segment->line_offset = segment->delta_V_min - segment->line_slope*segment->V_min;

	//compute the slope of the tangent line at the ends of feasible interval
	segment->d_delta_V_min = slope_of_tangent_line(Vbox->V_i_min, Vbox->V_j_min, branch_data, segment->tangent_line_type);
	segment->d_delta_V_max = slope_of_tangent_line(Vbox->V_i_max, Vbox->V_j_max, branch_data, segment->tangent_line_type);
}


/*Construct linear approximation in a nested loop. The outer loop increases the number of constraints and
the inner loop refines the offsets of intersection lines to reduce the approxation error*/
void compute_parameters_of_approximation(LF_Workspace* workspace, LF_Options* options) {
	BOOL finish_outer_loop = LF_FALSE, finish_inner_loop;
	int N_constructed_constraints, N_previously_constructed_constraints = workspace->constraint_counter, 
		iter_inner, N_constraints_max, N_new_constraints = 0, N_constraints_min;
	double max_error = 0.0, min_error = 0.0, max_error_previous = 1e10, width_relative;

	//make sure only one iteration is performed in case of non-convex region for line between load and gen
	if (!workspace->feas_reg.convex || workspace->branch_type == gen_gen)
		N_constraints_max = 1;
	else
		N_constraints_max = options->N_constraints_max;

	//choose the minimum number of linear constraints that will be considered
	if (options->computation_mode == 1)
		N_constructed_constraints = N_constraints_max;
	else if ((!workspace->feas_reg.top_left_corner_feasible || !workspace->feas_reg.bottom_right_corner_feasible) 
		&& N_constraints_max>=2)
		N_constructed_constraints = 2;
	else
		N_constructed_constraints = 1;
	N_constraints_min = N_constructed_constraints;

	//outer loop goes over the desired range of number of linear constraints
	workspace->errors[0] = 0.0;
	while (!finish_outer_loop) {
		compute_offsets_of_new_intersection_lines(workspace, N_constructed_constraints, N_new_constraints);
		N_constructed_constraints += N_new_constraints;

		//inner loop refines offsets of intersection lines for a given number of linear constraints to make the error smaller
		finish_inner_loop = LF_FALSE;
		iter_inner = 0;
		while (!finish_inner_loop) {
			//reset constraints counter
			workspace->constraint_counter = N_previously_constructed_constraints;

			//construct linear constraints
			construct_linear_constraints(workspace, options, N_constructed_constraints);

			//estimate errors
			if (options->approximation == conservative) {
				estimate_conservative_approximation_errors(workspace, N_constructed_constraints);
				min_error = min_element_in_array(workspace->errors, N_constructed_constraints + 1);
			}
			else {
				estimate_relaxed_approximation_errors(workspace, N_constructed_constraints);
				min_error = min_error_relaxed_approximation(workspace->errors, N_constructed_constraints);
			}
			max_error = max_element_in_array(workspace->errors, N_constructed_constraints + 1);

			//update counter of created constraints
			workspace->constraint_counter = workspace->constraint_counter + N_constructed_constraints;

			//check the stopping criteria of inner loop
			if (max_error < options->error_max || min_error / max_error >= options->ratio_threshold ||
				iter_inner == options->N_adjustments)
				finish_inner_loop = LF_TRUE;
			else {
				refine_offsets_of_existing_interseciton_lines(workspace, N_constructed_constraints, options->approximation);
				iter_inner++;
			}
		}

		//compute the relative 'width' of constraint in the middle
		width_relative = compute_width_of_middle_constraint(workspace, N_constructed_constraints);

		//check the stopping criteria of outer loop
		if (max_error < options->error_max || N_constructed_constraints == N_constraints_max ||
			(max_error_previous > max_error && max_error_previous - max_error < N_new_constraints*options->max_error_change) ||
			width_relative<0.45)
			finish_outer_loop = LF_TRUE;
		else if (max_error_previous < max_error && N_constructed_constraints>2) {
			N_constraints_max = N_constructed_constraints - N_new_constraints; //to make sure we exit after the next iteration
			N_new_constraints = -N_new_constraints;
		}
		else if (max_error_previous - max_error < max_error - options->error_max &&
			N_constraints_max - N_constructed_constraints>1 && N_constructed_constraints>N_constraints_min)
			N_new_constraints = 2;
		else
			N_new_constraints = 1;
		max_error_previous = max_error;
	}

	//record the estimate of the maximum error
	if (max_error > workspace->max_error_estimate)
		workspace->max_error_estimate = max_error;
}


/*coumpute the relative width of the constraint in the middle*/
double compute_width_of_middle_constraint(LF_Workspace* workspace, int N_constructed_constraints) {
	int index = N_constructed_constraints / 2;
	IntersectionLine line1, line2;
	if (N_constructed_constraints < 7)
		return 1.0;
	
	obtain_actual_intersection_line(workspace, index-1, index, &line1);
	obtain_actual_intersection_line(workspace, index, index + 1, &line2);

	return ((line2.offset - line1.offset)*N_constructed_constraints /
			(workspace->feas_reg.offset_max - workspace->feas_reg.offset_min));
}


/*compute the values of offsets of intersection lines given that their number has changed*/
void compute_offsets_of_new_intersection_lines(LF_Workspace* workspace, int N_constructed_constraints, int N_new_constraints) {
	double dif_offset = workspace->feas_reg.offset_max - workspace->feas_reg.offset_min;
	int N_constraints_to_be_constructed = N_constructed_constraints + N_new_constraints;
	IntersectionLine* lines = workspace->lines;
	
	//check if smart splitting is applicable
	if (N_constraints_to_be_constructed <= 2) {
		lines[0].offset = workspace->feas_reg.offset_min;
		lines[1].offset = (workspace->feas_reg.offset_min + workspace->feas_reg.offset_max) / 2.0;
		lines[2].offset = workspace->feas_reg.offset_max;
		return;
	}

	if (workspace->errors[0] == 0.0) {
		distribute_offsets_evenly(workspace, N_constructed_constraints, dif_offset);
		return;
	}

	if (N_new_constraints > 0)
		for (int i = N_constructed_constraints; i < N_constraints_to_be_constructed; i++)
			insert_intersection_line(lines, dif_offset, i);
	else
		for (int i = N_constructed_constraints; i > N_constraints_to_be_constructed; i--)
			remove_intersection_line(lines, dif_offset, i);
}


/*insert new intersection line based on the offsets of already constructed intersection lines*/
void insert_intersection_line(IntersectionLine* lines, double dif_offset, int N_constructed_constraints) {
	double buf1, buf2, buf3, coef, length_new;
	int N_constraints_to_be_constructed = N_constructed_constraints + 1,
		index_of_new_constraint = N_constraints_to_be_constructed / 2 - 1; //index of newly created constraint in a list of constraints

	//compute the 'length' of the constraint to be inserted
	if ((N_constraints_to_be_constructed - 1) % 2 == 0) //if it is even
		length_new = (lines[index_of_new_constraint + 2].offset - lines[index_of_new_constraint].offset) / 2.0;
	else //if it is odd
		length_new = lines[index_of_new_constraint + 1].offset - lines[index_of_new_constraint].offset;

	//get normalizing coefficient
	coef = dif_offset / (dif_offset + length_new);

	//compute new values of offsets
	//for constraints that are located in the list before the newly constructed constraint
	buf1 = lines[0].offset;
	for (int i = 0; i <= index_of_new_constraint; i++) {
		buf2 = lines[i + 1].offset;
		lines[i + 1].offset = lines[i].offset + (buf2 - buf1)*coef;
		buf1 = buf2;
	}
	//for the newly constructed constraint
	buf2 = lines[index_of_new_constraint + 2].offset;
	lines[index_of_new_constraint + 2].offset = lines[index_of_new_constraint + 1].offset + length_new*coef;
	//for constraints that are located in the list after the newly constructed constraint
	for (int i = index_of_new_constraint + 2; i < N_constraints_to_be_constructed; i++) {
		buf3 = lines[i + 1].offset;
		lines[i + 1].offset = lines[i].offset + (buf2 - buf1)*coef;
		buf1 = buf2;
		buf2 = buf3;
	}
}


/*remove intersection lines based on the offsets of already constructed intersection lines*/
void remove_intersection_line(IntersectionLine* lines, double dif_offset, int N_constructed_constraints) {
	double length_to_remove, coef, buf1, buf2;
	int index_constraint_to_remove = (N_constructed_constraints - 1) / 2;

	//compute the 'length' of constraints to be removed
	length_to_remove = lines[index_constraint_to_remove + 1].offset - lines[index_constraint_to_remove].offset;

	//get normalizing coefficient
	coef = dif_offset / (dif_offset - length_to_remove);

	//compute new values of offsets
	//for constraints that are located in the list before the removed constraint
	buf1 = lines[0].offset;
	for (int i = 0; i < index_constraint_to_remove; i++) {
		buf2 = lines[i + 1].offset;
		lines[i + 1].offset = lines[i].offset + (buf2 - buf1)*coef;
		buf1 = buf2;
	}
	//for constraints that are located in the list after the removed constraint
	buf1 = lines[index_constraint_to_remove + 1].offset;
	for (int i = index_constraint_to_remove + 1; i < N_constructed_constraints; i++) {
		buf2 = lines[i + 1].offset;
		lines[i].offset = lines[i - 1].offset + (buf2 - buf1)*coef;
		buf1 = buf2;
	}
}

/*refine the values of offsets of existing intersection lines in order to reduce the approximation error*/
void refine_offsets_of_existing_interseciton_lines(LF_Workspace* workspace, 
	int N_constructed_constraints, ApproximationType approximation) {
	double mean_error, buf1, buf2, buf3, dif_alpha;
	double *errors = workspace->errors;;
	IntersectionLine* lines = workspace->lines;

	dif_alpha = workspace->feas_reg.offset_max - workspace->feas_reg.offset_min;

	//check if smart splitting is applicable
	if (N_constructed_constraints == 1) {
		lines[0].offset = workspace->feas_reg.offset_min;
		lines[1].offset = (workspace->feas_reg.offset_min + workspace->feas_reg.offset_max) / 2.0;
		lines[2].offset = workspace->feas_reg.offset_max;
		return;
	}

	if (errors[0] == 0.0) {
		distribute_offsets_evenly(workspace, N_constructed_constraints, dif_alpha);
		return;
	}

	//refining offsets is slightly different for computing conservative and relaxed approximation due to different error representation
	if (approximation == conservative) {
		//compute the mean of errors
		buf1 = 0.0;
		for (int i = 0; i < N_constructed_constraints; i++)
			buf1 = buf1 + errors[i];
		mean_error = buf1 / N_constructed_constraints;
		//compute normalization coefficient
		buf1 = 0.0;
		for (int i = 0; i < N_constructed_constraints; i++)
			buf1 = buf1 + sqrt(mean_error / errors[i])*(lines[i + 1].offset - lines[i].offset);
		//refine values of offsets
		buf2 = lines[0].offset;
		for (int i = 0; i < N_constructed_constraints; i++) {
			buf3 = lines[i + 1].offset;
			lines[i + 1].offset = lines[i].offset + dif_alpha*sqrt(mean_error / errors[i])*(buf3 - buf2) / buf1;
			buf2 = buf3;
		}
	}
	else {
		//compute the mean of errors
		buf1 = (errors[0] + errors[N_constructed_constraints]) / 2;
		for (int i = 1; i < N_constructed_constraints; i++)
			buf1 = buf1 + errors[i];
		mean_error = buf1 / N_constructed_constraints;
		//compute normalization coefficient
		buf1 = 0.0;
		for (int i = 0; i < N_constructed_constraints; i++)
			buf1 = buf1 + sqrt(2 * mean_error / (errors[i] + errors[i + 1]))*(lines[i + 1].offset - lines[i].offset);
		//refine values of offsets
		buf2 = lines[0].offset;
		for (int i = 0; i < N_constructed_constraints; i++) {
			buf3 = lines[i + 1].offset;
			lines[i + 1].offset = lines[i].offset + dif_alpha*sqrt(2 * mean_error / (errors[i] + errors[i + 1]))*(buf3 - buf2) / buf1;
			buf2 = buf3;
		}
	}
}


/*compute the offsets of intersection lines such that the difference between two adjacent offsets is the same for all lines*/
void distribute_offsets_evenly(LF_Workspace* workspace, int N_constructed_constraints, double dif_alpha) {
	IntersectionLine* lines = workspace->lines;
	double alpha_min = workspace->feas_reg.offset_min;
	for (int i = 0; i <= N_constructed_constraints; i++)
		lines[i].offset = alpha_min + dif_alpha * i / N_constructed_constraints;
}


/*determine parameters of linear constraints for a given number of constraints*/
void construct_linear_constraints(LF_Workspace* workspace, LF_Options* options, int N_constraints) {
	if (workspace->branch_type == gen_gen) { //special case when line is between two gens
		workspace->errors[0] = 0.0;
		workspace->errors[1] = 0.0;
		workspace->plane_normals[workspace->constraint_counter * 3] = 0;
		workspace->plane_normals[workspace->constraint_counter * 3 + 1] = 0;
		workspace->plane_normals[workspace->constraint_counter * 3 + 2] = workspace->branch_data[9];
		workspace->b_planes[workspace->constraint_counter] = max_angle_for_given_current(workspace->V_limits.V_i_min,
			workspace->V_limits.V_j_min, workspace->branch_data) * workspace->plane_normals[workspace->constraint_counter * 3 + 2];
	}
	else if (!workspace->feas_reg.convex) { //special case when line is between gen and load and the feasible region is non-convex
		construct_approximation_to_nonconvex_region(workspace, options);
		estimate_approximation_error_for_nonconvex_region(workspace, options);
	}
	else {
		construct_intersection_lines(workspace, options, N_constraints);

		//shift_intersection_lines(workspace, options, N_constraints);

		compute_plane_normals(workspace, N_constraints, options->approximation);

		compute_plane_offsets(workspace, N_constraints, options);

		if (options->approximation == conservative && !workspace->feas_reg.bottom_right_corner_feasible && N_constraints > 1)
			update_plane_for_infeasible_box_corner(workspace, 0, N_constraints, options); 
		if (options->approximation == conservative && !workspace->feas_reg.top_left_corner_feasible && N_constraints > 1)
			update_plane_for_infeasible_box_corner(workspace, N_constraints - 1, N_constraints, options);
	}
}


/*construct conservative linear approximation for line between gen and load*/
void construct_approximation_to_nonconvex_region(LF_Workspace* workspace, LF_Options* options) {
	LineSegmentEdge segment = { .convex = LF_FALSE };
	double V_temp, delta_curve, sign = workspace->branch_data[9];
	int index_variable_V;

	compute_parameters_of_edge_line_segment(workspace->branch_data, &workspace->V_limits, &segment);

	if (options->approximation == conservative) {
		//use bisection to find the most distant pointon on the curve from the line
		V_temp = compute_value_by_bisection(segment.V_min, segment.V_max, options->eps_tolerance, workspace->branch_data,
			bisection_check_for_line_with_fixed_V, &segment);

		//compute the value of delta on the curve at the obtained point
		if (segment.tangent_line_type == fixed_V_i)
			delta_curve = max_angle_for_given_current(segment.V_fixed, V_temp, workspace->branch_data);
		else
			delta_curve = max_angle_for_given_current(V_temp, segment.V_fixed, workspace->branch_data);
	}
	else {
		V_temp = segment.V_min;
		delta_curve = segment.delta_V_min;
	}

	if (segment.tangent_line_type == fixed_V_i)
		index_variable_V = 1;
	else
		index_variable_V = 0;

	//compute parameters of approximation
	workspace->plane_normals[3 * workspace->constraint_counter + index_variable_V] = -sign*segment.line_slope;
	workspace->plane_normals[3 * workspace->constraint_counter + 1 - index_variable_V] = 0.0; //coefficient is zero for fixed V
	workspace->plane_normals[3 * workspace->constraint_counter + 2] = sign;
	workspace->b_planes[workspace->constraint_counter] = sign*(delta_curve - segment.line_slope*V_temp);
}


/*estimate approximation error for line between load and generator in case of non-convex feasible region*/
void estimate_approximation_error_for_nonconvex_region(LF_Workspace* workspace, LF_Options* options) {
	LineSegmentEdge segment = { .convex = LF_FALSE };
	double error_begin, error_end, V_temp;

	//fo conservative approximation we will only look at end points
	if (options->approximation == conservative) {
		error_begin = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_min, workspace->V_limits.V_j_min, 0);
		error_end = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_max, workspace->V_limits.V_j_max, 0);
		workspace->errors[0] = fmax(error_begin, error_end);
		return;
	}

	//for the relaxed case we find the most distant point from the approximation line and measure error there
	compute_parameters_of_edge_line_segment(workspace->branch_data, &workspace->V_limits, &segment);

	//use bisection to find the most distant pointon on the curve from the line
	V_temp = compute_value_by_bisection(segment.V_min, segment.V_max, options->eps_tolerance, workspace->branch_data,
		bisection_check_for_line_with_fixed_V, &segment);

	if (segment.tangent_line_type == fixed_V_i)
		workspace->errors[0] = error_value_in_ViVj_point(workspace, segment.V_fixed, V_temp, 0);
	else
		workspace->errors[0] = error_value_in_ViVj_point(workspace, V_temp, segment.V_fixed, 0);
}


/*compute all parameters related to intersection lines given their offsets*/
void construct_intersection_lines(LF_Workspace* workspace, LF_Options* options, int N_constraints) {
	double delta_begin, delta_end, delta_max = fabs(workspace->delta_max),
		*branch_data = workspace->branch_data;
	IntersectionLine* lines = workspace->lines;
	int n_lines_total = max_number_of_intersection_lines(N_constraints);

	for (int i = 0; i < n_lines_total; i++) {
		compute_intersections_of_line_with_Vbox(&lines[i], branch_data[8], &workspace->V_limits);

		//compute the values of delta at both ends of relevant line segment
		delta_begin = max_angle_for_given_current_slow(lines[i].point_begin.V_i, lines[i].point_begin.V_j, branch_data);
		delta_end = max_angle_for_given_current_slow(lines[i].point_end.V_i, lines[i].point_end.V_j, branch_data);

		//check if entire line segment is relevant for approximation
		if (fabs(delta_begin) > delta_max || fabs(delta_end) > delta_max) {
			shrink_intersection_line_segment(delta_begin, delta_end, workspace, &lines[i]);
			//recompute the values of delta at both ends of relevant line segment
			delta_begin = max_angle_for_given_current_slow(lines[i].point_begin.V_i, lines[i].point_begin.V_j, branch_data);
			delta_end = max_angle_for_given_current_slow(lines[i].point_end.V_i, lines[i].point_end.V_j, branch_data);
		}

		//compute the slope of relaxed (outer) intersection line
		if (lines[i].point_begin.V_i != lines[i].point_end.V_i)
			lines[i].slope = (delta_end - delta_begin) / (lines[i].point_end.V_i - lines[i].point_begin.V_i);
		else
			lines[i].slope = 0;

		//compute the value of ddelta/dV_i at V_i_begin
		lines[i].tangent_slope_begin = slope_of_tangent_line(lines[i].point_begin.V_i, lines[i].point_begin.V_j, branch_data, intersection_line);
		//compute the value of ddelta/dV_i at V_i_end
		lines[i].tangent_slope_end = slope_of_tangent_line(lines[i].point_end.V_i, lines[i].point_end.V_j, branch_data, intersection_line);
	}

	workspace->slope_all_lines = slope_of_intersection_lines(lines, n_lines_total);
	compute_deltas_at_ends_of_intersection_lines(workspace, n_lines_total, options);
}


/*compute parameters of points at which a line intersects with Vbox*/
void compute_intersections_of_line_with_Vbox(IntersectionLine* line, double line_slope, VoltageLimits* Vbox) {
	//beginning of the line segment

	if (line_slope * Vbox->V_i_min + line->offset >= Vbox->V_j_min) {
		line->point_begin.V_i = Vbox->V_i_min;
		line->point_begin.V_j = line_slope * Vbox->V_i_min + line->offset;
	}
	else {
		line->point_begin.V_i = (Vbox->V_j_min - line->offset) / line_slope;
		line->point_begin.V_j = Vbox->V_j_min;
	}
	//end of the line segment
	if (line_slope * Vbox->V_i_max + line->offset <= Vbox->V_j_max) {
		line->point_end.V_i = Vbox->V_i_max;
		line->point_end.V_j = line_slope * Vbox->V_i_max + line->offset;
	}
	else {
		line->point_end.V_i = (Vbox->V_j_max - line->offset) / line_slope;
		line->point_end.V_j = Vbox->V_j_max;
	}
}


/*update values of end points of the intersection line to make sure that all points between them can have I>=I_max for delta<=delta_max*/
void shrink_intersection_line_segment(double delta_begin, double delta_end, LF_Workspace* workspace, IntersectionLine* line) {
	double V_i_out;

	if (fabs(delta_begin) > fabs(workspace->delta_max) &&
		fabs(delta_end) > fabs(workspace->delta_max)) { //then the entire line is irrelevant for determining the slope
		line->slope = 0.0;
		return;
	}

	//find point(V_i, V_j) along the line at which delta on the surface is the maximum allowed
	V_i_out = update_end_of_intersection_line_segment(workspace, line->offset);
	if (fabs(delta_begin) > fabs(workspace->delta_max)) { //then only second part of the line is relevant
		line->point_begin.V_i = V_i_out;
		line->point_begin.V_j = V_i_out*workspace->branch_data[8] + line->offset;
	}
	else { //then only second part of the line is relevant
		line->point_end.V_i = V_i_out;
		line->point_end.V_j = V_i_out*workspace->branch_data[8] + line->offset;
	}
}


/*compute the slope that all intersection lines will have in (V_i, delta) plane*/
double slope_of_intersection_lines(IntersectionLine* lines, int n_lines_total) {
	double slope_all_lines = 0.0, length_squared, sum_of_line_length_squares = 0.0;

	//different 'outer' intersection lines will have different weights of their slopes, which are made proportional to squares of line lengths
	for (int i = 1; i < n_lines_total - 1; i++) {
		if (lines[i].slope != 0.0) {
			length_squared = (lines[i].point_end.V_i - lines[i].point_begin.V_i)*(lines[i].point_end.V_i - lines[i].point_begin.V_i) +
				(lines[i].point_end.V_j - lines[i].point_begin.V_j)*(lines[i].point_end.V_j - lines[i].point_begin.V_j);
			slope_all_lines = slope_all_lines + length_squared*lines[i].slope;
			sum_of_line_length_squares = sum_of_line_length_squares + length_squared;
		}
	}
	//normalize weights
	if (sum_of_line_length_squares > 0.0)
		slope_all_lines = slope_all_lines / sum_of_line_length_squares;

	return slope_all_lines;
}


/*for each intersection line, compute the value of delta at either line beginning or end*/
void compute_deltas_at_ends_of_intersection_lines(LF_Workspace* workspace, int n_lines_total, LF_Options* options) {
	IntersectionLine* lines = workspace->lines;
	double V_i_temp, V_j_temp, delta_temp, *branch_data = workspace->branch_data, slope_all_lines = workspace->slope_all_lines;

	//handle first and last lines separately as they are special cases
	lines[0].point_end.delta = max_angle_for_given_current(lines[0].point_end.V_i, lines[0].point_end.V_j, branch_data);
	lines[n_lines_total - 1].point_end.delta = max_angle_for_given_current(lines[n_lines_total - 1].point_end.V_i,
		lines[n_lines_total - 1].point_end.V_j, branch_data);
	if (options->approximation == relaxed) {
		lines[0].point_begin.delta = lines[0].point_end.delta -
			slope_all_lines*(lines[0].point_end.V_i - lines[0].point_begin.V_i);
		lines[n_lines_total - 1].point_begin.delta = lines[n_lines_total - 1].point_end.delta -
			slope_all_lines*(lines[n_lines_total - 1].point_end.V_i - lines[n_lines_total - 1].point_begin.V_i);
	}
	else {
		lines[0].point_begin.delta = lines[0].point_end.delta;
		lines[n_lines_total - 1].point_begin.delta = lines[n_lines_total - 1].point_end.delta;
	}

	//deltas at end points of all other lines
	for (int i = 1; i < n_lines_total - 1; i++) {
		if (options->approximation == conservative)
			V_i_temp = point_on_conservative_intersection_line(&lines[i], branch_data, slope_all_lines, options);
		else
			V_i_temp = point_on_relaxed_intersection_line(&lines[i], workspace);

		//compute V_j and delta at the obtained point
		V_j_temp = branch_data[8] * V_i_temp + lines[i].offset;
		delta_temp = max_angle_for_given_current(V_i_temp, V_j_temp, branch_data);

		//compute values of delta at end of line segment
		lines[i].point_begin.delta = delta_temp - slope_all_lines*(V_i_temp - lines[i].point_begin.V_i);
		lines[i].point_end.delta = delta_temp - slope_all_lines*(V_i_temp - lines[i].point_end.V_i);
	}
}


/*compute value of V_i at the point at which the conservative intersection line touches the surface*/
double point_on_conservative_intersection_line(IntersectionLine* line, double *branch_data, 
	double slope_all_lines, LF_Options *options) {
	double data[3] = { line->offset , slope_all_lines, branch_data[8] }, V_i;

	//check which value of V_i will be used for getting value of delta at intersection line
	if (slope_all_lines<fmax(line->tangent_slope_begin, line->tangent_slope_end) &&
		slope_all_lines>fmin(line->tangent_slope_begin, line->tangent_slope_end)) {
		//find the point at which the tangent line has the same slope as intersection lines
		//choose initial point for NR
		V_i = line->point_begin.V_i + (slope_all_lines - line->tangent_slope_begin)*(line->point_end.V_i - line->point_begin.V_i) / 
			(line->tangent_slope_end - line->tangent_slope_begin);
		if (V_i<line->point_begin.V_i || V_i>line->point_end.V_i)
			V_i = (line->point_begin.V_i + line->point_end.V_i) / 2;
		//solve with NR
		V_i = solve_equation_by_Newton_Raphson(V_i, options, branch_data, update_in_NR_intersection_line, data);
		if (V_i == 0.0) //if NR did not converge, solve with bisection
			V_i = compute_value_by_bisection(line->point_begin.V_i, line->point_end.V_i, options->eps_tolerance, branch_data,
				bisection_check_for_intersection_line, &data);
		else if (V_i < line->point_begin.V_i)
			V_i = line->point_begin.V_i;
		else if (V_i > line->point_end.V_i)
			V_i = line->point_end.V_i;
	}
	else if (fabs(line->tangent_slope_begin - slope_all_lines) < fabs(line->tangent_slope_end - slope_all_lines))
		V_i = line->point_begin.V_i;
	else
		V_i = line->point_end.V_i;

	return V_i;
}


/*compute value of V_i at the point at which the initial relaxed intersection line touches the surface*/
double point_on_relaxed_intersection_line(IntersectionLine* line, LF_Workspace* workspace) {
	double delta_surface_end, delta_line_end, offset;
		
	offset = max_angle_for_given_current(line->point_begin.V_i, line->point_begin.V_j, workspace->branch_data) -
		line->point_begin.V_i*workspace->slope_all_lines;

	delta_line_end = line->point_end.V_i*workspace->slope_all_lines + offset;
	delta_surface_end = max_angle_for_given_current(line->point_end.V_i, line->point_end.V_j, workspace->branch_data);

	if (workspace->branch_data[9] * delta_line_end > workspace->branch_data[9] * delta_surface_end)
		return line->point_begin.V_i;
	else
		return line->point_end.V_i;
}


/*compute plane normals for all planes*/
void compute_plane_normals(LF_Workspace* workspace, int N_constraints, ApproximationType approximation) {
	IntersectionLine* lines = workspace->lines;
	int counter = workspace->constraint_counter, index_second_line;
	double coef_normalized, coef_ViVj = workspace->branch_data[8], *plane_normals = workspace->plane_normals,
		vector_1[3] = { 0 }, vector_2[3] = { 0 }, plane_normal[3] = { 0 };

	//to find the equation of the plane, we need two vectors and a point that all lie in this plane
	for (int i = 0; i < N_constraints; i++) {
		//determine index of second intersection line that will be used for parameters of current plane
		if (N_constraints == 1)
			index_second_line = 2;
		else
			index_second_line = i + 1;

		//record parameters of direction vector of intersection lines
		if (approximation == conservative &&
			((i == 0 && !workspace->feas_reg.bottom_right_corner_feasible) ||
				(i == N_constraints - 1 && !workspace->feas_reg.top_left_corner_feasible)))
			fill_vector_with_values(1, coef_ViVj, 0.0, vector_1);
		else
			fill_vector_with_values(1, coef_ViVj, workspace->slope_all_lines, vector_1);

		//get vector parallel to line segment connecting two points on adjacent intersection lines
		if (approximation == conservative)
			fill_vector_with_values(lines[index_second_line].point_begin.V_i - lines[i].point_begin.V_i,
				lines[index_second_line].point_begin.V_j - lines[i].point_begin.V_j,
				lines[index_second_line].point_begin.delta - lines[i].point_begin.delta, vector_2);
		else
			fill_vector_with_values(lines[index_second_line].point_end.V_i - lines[i].point_end.V_i,
				lines[index_second_line].point_end.V_j - lines[i].point_end.V_j,
				lines[index_second_line].point_end.delta - lines[i].point_end.delta, vector_2);

		//compute plane normal
		compute_vector_product(vector_1, vector_2, plane_normal);

		//compute normalization coefficient
		coef_normalized = plane_normal[2] * workspace->branch_data[9];

		//normalize plane normal such that coefficient for delta is 1 (for upper part) or -1 (for lower part)
		for (int j = 0; j < 3; j++)
			plane_normals[(i + counter) * 3 + j] = plane_normal[j] / coef_normalized;
	}
}


/*compute the offsets of planes (conservative or relaxed)*/
void compute_plane_offsets(LF_Workspace* workspace, int N_constraints, LF_Options* options) {
	IntersectionLine* lines = workspace->lines;
	double point_on_plane[3] = { 0 }, delta_shift, *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes;
	int counter = workspace->constraint_counter, index_last_line;
	ApproximationType approximation = options->approximation;

	index_last_line = max_number_of_intersection_lines(N_constraints) - 1;

	//to find the equation of the plane, we need two vectors and a point that all lie in this plane
	for (int i = 0; i < N_constraints; i++) {
		//compute parameters of a point on the plane
		if (approximation == relaxed)
			delta_shift = delta_shift_plane_relaxed(workspace, i, N_constraints, options);
		else
			delta_shift = 0;
		
		if (i == N_constraints - 1)
			fill_vector_with_values(lines[index_last_line].point_end.V_i, lines[index_last_line].point_end.V_j,
				lines[index_last_line].point_end.delta + delta_shift, point_on_plane);
		else
			fill_vector_with_values(lines[i].point_end.V_i, lines[i].point_end.V_j,
				lines[i].point_end.delta + delta_shift, point_on_plane);

		//compute value of b in Ax = b for this plane
		b_planes[i + counter] = 0.0;
		for (int j = 0; j < 3; j++)
			b_planes[i + counter] = b_planes[i + counter] + plane_normals[(i + counter) * 3 + j] * point_on_plane[j];
	}
}


/*compute the required shift in delta to ensure that a given plane represents a relaxed approximation to the surface*/
double delta_shift_plane_relaxed(LF_Workspace* workspace, int index, int N_constraints, LF_Options* options) {
	PlaneType plane_type;
	int line_indexes[2] = { index, 0 };
	double delta_shift_begin = 0, delta_shift_end = 0, delta_shift_max;

	//record indices of the intersection lines relevant for this plane
	if (N_constraints == 1)
		line_indexes[1] = 2;
	else
		line_indexes[1] = index + 1;

	//check if this plane is the outermost one and the corresponding box corner is infeasible
	if (index == 0 && !workspace->feas_reg.bottom_right_corner_feasible)
		plane_type = first;
	else if (index == N_constraints - 1 && !workspace->feas_reg.top_left_corner_feasible)
		plane_type = last;
	else
		plane_type = ordinary;

	//compute required shift in delta for the beginning of the plane
	delta_shift_begin = delta_shift_one_side_of_plane(workspace, line_indexes, options, LF_TRUE, plane_type);

	//compute required shift in delta for the end of the plane
	delta_shift_end = delta_shift_one_side_of_plane(workspace, line_indexes, options, LF_FALSE, plane_type);

	//compute the largest required shift in delta
	if (workspace->branch_data[9] == 1.0)
		delta_shift_max = fmax(delta_shift_begin, delta_shift_end);
	else
		delta_shift_max = fmin(delta_shift_begin, delta_shift_end);

	return delta_shift_max;
}


/*compute the shift in delta based on the edge located at either beginning or end of the box*/
double delta_shift_one_side_of_plane(LF_Workspace* workspace, int* line_indexes, LF_Options* options,
	BOOL bottom_left_corner, PlaneType plane_type) {
	LineSegmentEdge segment = { .convex = LF_TRUE };
	IntersectionLine* lines = workspace->lines;
	Point3D point1, point2, point_corner;
	double offset_corner, delta_surface_corner, *branch_data = workspace->branch_data, delta_shift_temp[3];

	//record values of relevant points on intersection lines
	if (bottom_left_corner) {
		point1 = lines[line_indexes[0]].point_begin;
		point2 = lines[line_indexes[1]].point_begin;
	}
	else {
		point1 = lines[line_indexes[0]].point_end;
		point2 = lines[line_indexes[1]].point_end;
	}

	//compute offset of line passing through corner 
	if (bottom_left_corner)
		offset_corner = workspace->V_limits.V_j_min - branch_data[8] * workspace->V_limits.V_i_min;
	else
		offset_corner = workspace->V_limits.V_j_max - branch_data[8] * workspace->V_limits.V_i_max;

	//compare with offsets of two intersection lines
	if (lines[line_indexes[0]].offset >= offset_corner || lines[line_indexes[1]].offset <= offset_corner) { //both lines intersect same edge of Vbox
		compute_parameters_of_plane_edge_segment(branch_data, &segment, &point1, &point2);
		return delta_shift_on_edge(workspace, &segment, options, plane_type);
	}

	if (plane_type == ordinary)
		plane_type = corner;

	//fill out parameters of the relevant corner point
	if (bottom_left_corner) {
		point_corner.V_i = workspace->V_limits.V_i_min;
		point_corner.V_j = workspace->V_limits.V_j_min;
		compute_delta_plane_in_Vbox_corner(workspace, line_indexes[0], &lines[line_indexes[0]].point_begin, &point_corner);
	}
	else {
		point_corner.V_i = workspace->V_limits.V_i_max;
		point_corner.V_j = workspace->V_limits.V_j_max;
		compute_delta_plane_in_Vbox_corner(workspace, line_indexes[0], &lines[line_indexes[0]].point_end, &point_corner);
	}

	//deal with shift required to move corner inside the feasible region
	delta_surface_corner = max_angle_for_given_current(point_corner.V_i, point_corner.V_j, branch_data);
	delta_shift_temp[0] = delta_surface_corner - point_corner.delta;

	//deal with shift for edge segment from first intersection line to corner
	compute_parameters_of_plane_edge_segment(branch_data, &segment, &point1, &point_corner);
	delta_shift_temp[1] = delta_shift_on_edge(workspace, &segment, options, plane_type);

	//deal with shift for edge segment from corner to second intersection line
	compute_parameters_of_plane_edge_segment(branch_data, &segment, &point2, &point_corner);
	delta_shift_temp[2] = delta_shift_on_edge(workspace, &segment, options, plane_type);

	//select the largest shift
	if (workspace->branch_data[9] == 1.0)
		return max_element_in_array(delta_shift_temp, 3);
	else
		return min_element_in_array(delta_shift_temp, 3);
}


/*compute parameters of the line segment that represents the intersection of a plane with the edge of the box*/
void compute_parameters_of_plane_edge_segment(double *branch_data, LineSegmentEdge* segment,
	Point3D* point1, Point3D* point2) {
	//record the values of fixed V, V_min and V_max
	if (point1->V_i == point2->V_i) {
		segment->tangent_line_type = fixed_V_i;
		segment->V_fixed = point1->V_i;
		segment->V_min = fmin(point1->V_j, point2->V_j);
		segment->V_max = fmax(point1->V_j, point2->V_j);
	}
	else {
		segment->tangent_line_type = fixed_V_j;
		segment->V_fixed = point1->V_j;
		segment->V_min = fmin(point1->V_i, point2->V_i);
		segment->V_max = fmax(point1->V_i, point2->V_i);
	}

	//record the values of delta at the ends of the interval
	if (point1->V_i + point1->V_j < point2->V_i + point2->V_j) {
		segment->delta_V_min = point1->delta;
		segment->delta_V_max = point2->delta;
	}
	else {
		segment->delta_V_min = point2->delta;
		segment->delta_V_max = point1->delta;
	}

	//compute slope of the line connecting end points of the interval
	if (segment->V_max != segment->V_min)
		segment->line_slope = (segment->delta_V_max - segment->delta_V_min) / (segment->V_max - segment->V_min);
	else
		segment->line_slope = 0.0;
}


/*compute value of delta on the plane in either (V_i_min, V_j_min) or (V_i_max, V_j_max) corner of Vbox*/
void compute_delta_plane_in_Vbox_corner(LF_Workspace* workspace, int index, Point3D* point_line, Point3D* point_corner) {
	double plane_offset, point_on_plane[3], *plane_normals = workspace->plane_normals;
	int ind = index + workspace->constraint_counter;

	//record point on the plane
	fill_vector_with_values(point_line->V_i, point_line->V_j, point_line->delta, point_on_plane);

	//compute value of b in Ax = b for this plane
	plane_offset = plane_normals[ind * 3] * point_on_plane[0] + plane_normals[ind * 3 + 1] * point_on_plane[1] +
		plane_normals[ind * 3 + 2] * point_on_plane[2];

	//compute value of delta on the plane in the corner
	point_corner->delta = (plane_offset - plane_normals[ind * 3] * point_corner->V_i -
		plane_normals[ind * 3 + 1] * point_corner->V_j) / plane_normals[ind * 3 + 2];
}


/*compute by how much the plane should be shifted to become relaxed based on a given edge of Vbox*/
double delta_shift_on_edge(LF_Workspace *workspace, LineSegmentEdge* segment, LF_Options* options, PlaneType plane_type) {
	double V_temp, delta_surface, delta_line, delta_shift = 0, *branch_data = workspace->branch_data, V_min, V_max;

	//compute limits on V to check against
	if (plane_type == corner) {
		V_min = segment->V_min;
		V_max = segment->V_max;
	}
	else if (segment->tangent_line_type == fixed_V_i) {
		V_min = fmax(workspace->V_limits.V_j_min, segment->V_fixed*branch_data[8] + workspace->feas_reg.offset_min);
		V_max = fmin(workspace->V_limits.V_j_max, segment->V_fixed*branch_data[8] + workspace->feas_reg.offset_max);
	}
	else {
		V_min = fmax(workspace->V_limits.V_i_min, (segment->V_fixed - workspace->feas_reg.offset_max) / branch_data[8]);
		V_max = fmin(workspace->V_limits.V_i_max, (segment->V_fixed - workspace->feas_reg.offset_min) / branch_data[8]);
	}

	//choose initial point for NR algorithm
	if (plane_type == ordinary || plane_type == corner)
		V_temp = (segment->V_min + segment->V_max) / 2;
	else if ((plane_type == first && segment->tangent_line_type == fixed_V_i) || (plane_type == last && segment->tangent_line_type == fixed_V_j))
		V_temp = segment->V_min + (segment->V_max - segment->V_min) / 3.0;
	else
		V_temp = segment->V_min + (segment->V_max - segment->V_min)*2.0 / 3.0;
	//determine the most distant point on on the curve from the line
	V_temp = solve_equation_by_Newton_Raphson(V_temp, options, branch_data, update_in_NR_line_with_fixed_V, segment);
	//check if the algorithm converged to a desired point and if not, use bisection, which is guaranteed to converge
	if (V_temp == 0.0)
		V_temp = compute_value_by_bisection(V_min, V_max, options->eps_tolerance, branch_data,
			bisection_check_for_line_with_fixed_V, segment);
	else if (V_temp < V_min)
		V_temp = V_min;
	else if (V_temp > V_max)
		V_temp = V_max;

	//compute the value of delta on the line connecting segment end points
	delta_line = segment->delta_V_min + (V_temp - segment->V_min)*segment->line_slope;

	//compute delta_on surface
	if (segment->tangent_line_type == fixed_V_i)
		delta_surface = max_angle_for_given_current(segment->V_fixed, V_temp, branch_data);
	else
		delta_surface = max_angle_for_given_current(V_temp, segment->V_fixed, branch_data);

	//check if the shift is in the correct direction
	if (branch_data[9] * delta_surface > branch_data[9] * delta_line)
		delta_shift = delta_surface - delta_line;
	return delta_shift;
}


/*update parameters of plane if a box corner is infeasible to make sure the approximation is conservative*/
void update_plane_for_infeasible_box_corner(LF_Workspace *workspace, int index_outer, int N_constraints, LF_Options *options) {
	double *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes, coef_normalized,
		vector_1[3], vector_2[3], plane_normal[3];
	int index_inner, ind_line, index, counter = workspace->constraint_counter;
	IntersectionLine line;
	Point3D point, point_line;

	//compute necessary indices
	if (index_outer == 0) {
		index_inner = 1;
		ind_line = 2;
	}
	else {
		index_inner = index_outer - 1;
		ind_line = index_inner;
	}
	index = (index_inner + counter);

	//obtain parameters of the point that has maximum violation of constraint
	obtain_most_distant_point(workspace, index_outer + counter, index_inner + counter, &point);

	if (point.V_i == 0.0) //approximation is conservative, no change required
		return;

	//obtain the point on the intersection of the surface, outer plane and line perpendicular to 'regular' intersection line
	if (N_constraints == 2 && !workspace->feas_reg.top_left_corner_feasible &&
		!workspace->feas_reg.bottom_right_corner_feasible) {
		point = workspace->lines[1].point_end;
		point.delta = max_angle_for_given_current(point.V_i, point.V_j, workspace->branch_data);
	}
	else
		obtain_point_on_updated_plane(workspace, index_outer, options, &point);

	//fill out parameters of a vector on the updated plane
	if (N_constraints == 2 && ((index_outer == 0 && !workspace->feas_reg.top_left_corner_feasible) ||
		(index_outer == 1 && !workspace->feas_reg.bottom_right_corner_feasible)))
		fill_vector_with_values(1, workspace->branch_data[8], 0, vector_1);
	else
		fill_vector_with_values(1, workspace->branch_data[8], workspace->slope_all_lines, vector_1);
	//fill oout parameters of a second point on the update plane that belongs to the intersection line
	if (N_constraints != 3)
		point_line = workspace->lines[ind_line].point_end;
	else {
		obtain_actual_intersection_line(workspace, 2 - index_outer + counter, index_inner + counter, &line);
		point_line = line.point_end;
	}
	fill_vector_with_values(point.V_i - point_line.V_i, point.V_j - point_line.V_j, point.delta - point_line.delta, vector_2);
	//recompute plane normal and offset
	compute_vector_product(vector_1, vector_2, plane_normal);
	fill_vector_with_values(point.V_i, point.V_j, point.delta, vector_1);
	coef_normalized = plane_normal[2] * workspace->branch_data[9];
	b_planes[index] = 0.0;
	for (int j = 0; j < 3; j++) {
		plane_normals[index * 3 + j] = plane_normal[j] / coef_normalized;
		b_planes[index] += plane_normals[index * 3 + j] * vector_1[j];
	}
}


/*compute parameters of a point that lies on the updated plane*/
void obtain_point_on_updated_plane(LF_Workspace *workspace, int index_outer, LF_Options *options, Point3D *point) {
	double *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes,
		data[5], V_limit, V_temp;
	int ind_line_outer, counter = workspace->constraint_counter;

	//compute necessary indices
	if (index_outer == 0) {
		ind_line_outer = 0;
		data[4] = 1;
	}
	else {
		ind_line_outer = index_outer + 1;
		data[4] = -1;
	}

	//prepare some data
	data[0] = -1.0 / workspace->branch_data[8];
	data[1] = point->V_j - data[0] * point->V_i;
	data[2] = (plane_normals[(index_outer + counter) * 3] + plane_normals[(index_outer + counter) * 3 + 1] * data[0]) /
		plane_normals[(index_outer + counter) * 3 + 2];
	data[3] = (plane_normals[(index_outer + counter) * 3 + 1] * data[1] - b_planes[index_outer + counter]) /
		plane_normals[(index_outer + counter) * 3 + 2];

	//compute the value of V_i at the interseciton of outermost intersection line and line perpendicular to it
	V_limit = (data[1] - workspace->lines[ind_line_outer].offset) / (workspace->branch_data[8] - data[0]);
	V_temp = solve_equation_by_Newton_Raphson(point->V_i, options, workspace->branch_data,
		update_in_NR_intersection_of_plane_and_surface, data);
	if (V_temp<fmin(V_limit, point->V_i) || V_temp>fmax(V_limit, point->V_i)) //if Newton Raphson failed, use bisection
		V_temp = compute_value_by_bisection(fmin(V_limit, point->V_i), fmax(V_limit, point->V_i),
			options->eps_tolerance, workspace->branch_data, bisection_intersection_of_plane_and_surface, data);
	point->V_i = V_temp;
	point->V_j = data[0] * point->V_i + data[1];
	point->delta = (b_planes[index_outer + counter] - plane_normals[(index_outer + counter) * 3] * point->V_i -
		plane_normals[(index_outer + counter) * 3 + 1] * point->V_j) / plane_normals[(index_outer + counter) * 3 + 2];
}


/*compute parameters of the point on the intersection line that violates the conservative approximation the most*/
void obtain_most_distant_point(LF_Workspace *workspace, int index_outer, int index_inner, Point3D *point) {
	double distance_old, distance_new;
	IntersectionLine line;
	Point3D point_temp;
	int N_points = (int)ceil((workspace->V_limits.V_i_max - workspace->V_limits.V_i_min) / 0.02);

	if (N_points < 5)
		N_points = 5; //check at least five points

	//determine parameters of the intersection line
	obtain_actual_intersection_line(workspace, index_outer, index_inner, &line);

	//check points along the line and record the one which is the most distant from the surface and not conservative
	distance_old = 0.0;
	for (int i = N_points; i > 0; i--) {
		point_temp.V_i = line.point_begin.V_i + (line.point_end.V_i - line.point_begin.V_i)*i / N_points;
		point_temp.V_j = line.point_begin.V_j + (line.point_end.V_j - line.point_begin.V_j)*i / N_points;
		point_temp.delta = line.point_begin.delta + (line.point_end.delta - line.point_begin.delta)*i / N_points;
		distance_new = (point_temp.delta - max_angle_for_given_current(point_temp.V_i, point_temp.V_j, workspace->branch_data))*
			workspace->branch_data[9];
		if (distance_new > 0 && distance_new > distance_old) {
			point->V_i = point_temp.V_i;
			point->V_j = point_temp.V_j;
			point->delta = point_temp.delta;
			distance_old = distance_new;
		}
	}
	if (distance_old == 0.0)
		point->V_i = 0.0;
}


/*compute parameters of the actual intersection line of two planes*/
void obtain_actual_intersection_line(LF_Workspace *workspace, int index_outer, int index_inner, IntersectionLine *line) {
	double *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes;
	Point3D point_1, point_2;

	//determine parameters of the intersection line
	//first point on the line
	point_1.V_i = 0.0;
	point_1.V_j = intersection_of_2D_lines(plane_normals[index_outer * 3 + 1], plane_normals[index_outer * 3 + 2],
		plane_normals[index_inner * 3 + 1], plane_normals[index_inner * 3 + 2], b_planes[index_outer], b_planes[index_inner]);
	point_1.delta = (b_planes[index_outer] - plane_normals[index_outer * 3 + 1] * point_1.V_j) / plane_normals[index_outer * 3 + 2];
	line->offset = point_1.V_j; //offset of this line in (V_i, V_j) plane
	//second point on the line
	if (fabs(point_1.V_j) > 1e-5)
		point_2.V_j = 0.0;
	else
		point_2.V_j = 1.0;
	point_2.V_i = intersection_of_2D_lines(plane_normals[index_outer * 3], plane_normals[index_outer * 3 + 2],
		plane_normals[index_inner * 3], plane_normals[index_inner * 3 + 2], b_planes[index_outer] -
		plane_normals[index_outer * 3 + 1] * point_2.V_j, b_planes[index_inner] - plane_normals[index_inner * 3 + 1] * point_2.V_j);
	point_2.delta = (b_planes[index_outer] - plane_normals[index_outer * 3] * point_2.V_i -
		plane_normals[index_outer * 3 + 1] * point_2.V_j) / plane_normals[index_outer * 3 + 2];

	//compute slope of this line in (Vi, Vj) plane
	line->slope = (point_2.V_j - point_1.V_j) / (point_2.V_i - point_1.V_i);

	compute_intersections_of_line_with_Vbox(line, line->slope, &workspace->V_limits);
	line->point_begin.delta = (b_planes[index_outer] - plane_normals[index_outer * 3] * line->point_begin.V_i -
		plane_normals[index_outer * 3 + 1] * line->point_begin.V_j) / plane_normals[index_outer * 3 + 2];
	line->point_end.delta = (b_planes[index_outer] - plane_normals[index_outer * 3] * line->point_end.V_i -
		plane_normals[index_outer * 3 + 1] * line->point_end.V_j) / plane_normals[index_outer * 3 + 2];
}


/*estimate maximum values of conservative approximation errors for each of constructed constraints*/
void estimate_conservative_approximation_errors(LF_Workspace* workspace, int N_constraints) {
	IntersectionLine middle_line, *lines = workspace->lines;
	double offset_corner_min, offset_corner_max, error_begin, error_end, error_corner, coef_ViVj = workspace->branch_data[8],
		*errors = workspace->errors;

	//compute offsets of the lines passing through (V_i_min, V_j_min) and (V_i_max, V_j_max) corners
	offset_corner_min = workspace->V_limits.V_j_min - workspace->branch_data[8] * workspace->V_limits.V_i_min;
	offset_corner_max = workspace->V_limits.V_j_max - workspace->branch_data[8] * workspace->V_limits.V_i_max;

	for (int i = 0; i < N_constraints; i++) {
		//compute the offset of the line that we will use for determining the error 
		if (N_constraints == 1)
			middle_line.offset = lines[1].offset;
		else
			middle_line.offset = (lines[i].offset + lines[i + 1].offset) / 2.0;

		compute_intersections_of_line_with_Vbox(&middle_line, coef_ViVj, &workspace->V_limits);

		//compute error at the beginning of the line
		error_begin = error_value_in_ViVj_point(workspace, middle_line.point_begin.V_i, middle_line.point_begin.V_j, i);

		//check if V_i_min V_j_min corner of the Vbox belongs to this plane
		if (offset_corner_min >= lines[i].offset && offset_corner_min <= lines[i + 1].offset)
			error_corner = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_min, workspace->V_limits.V_j_min, i);
		else
			error_corner = 0;

		error_begin = fmax(error_begin, error_corner);

		//compute error at the end of the line
		error_end = error_value_in_ViVj_point(workspace, middle_line.point_end.V_i, middle_line.point_end.V_j, i);

		//check if V_i_max V_j_max corner of the Vbox belongs to this plane
		if (offset_corner_max >= lines[i].offset && offset_corner_max <= lines[i + 1].offset)
			error_corner = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_max, workspace->V_limits.V_j_max, i);
		else
			error_corner = 0;

		error_end = fmax(error_end, error_corner);

		//record max error estimate for this constraint
		errors[i] = fmax(error_begin, error_end);
	}
	errors[N_constraints] = errors[N_constraints - 1]; //this is just a dummy term needed for convenience
}



/*estimate values of relaxed approximation errors for each of constructed constraints that will be used for refining line offsets*/
void estimate_relaxed_approximation_errors(LF_Workspace* workspace, int N_constraints) {
	IntersectionLine line_temp, *lines = workspace->lines;
	BOOL compute_error_at_corners;
	double line_slope, coef_ViVj = workspace->branch_data[8], *errors = workspace->errors, *b_planes = workspace->b_planes,
		abs_slope = fabs(workspace->slope_all_lines), *plane_normals = workspace->plane_normals;
	int ind_plane, ind_line, counter = workspace->constraint_counter, index;

	//check if we need to compute error in the corner for the case of N_constraints=1
	if (N_constraints == 1 &&
		((workspace->feas_reg.bottom_right_corner_feasible && workspace->feas_reg.top_left_corner_feasible) ||
			(!workspace->feas_reg.bottom_right_corner_feasible && !workspace->feas_reg.top_left_corner_feasible)))
		compute_error_at_corners = LF_TRUE;
	else 
		compute_error_at_corners = LF_FALSE;
	
	for (int i = 0; i <= N_constraints; i++) {
		if ((i == 0 && workspace->feas_reg.bottom_right_corner_feasible) || compute_error_at_corners)
			errors[i] = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_max, workspace->V_limits.V_j_min, 0);
		else if ((i == N_constraints && workspace->feas_reg.top_left_corner_feasible) || compute_error_at_corners)
			errors[i] = error_value_in_ViVj_point(workspace, workspace->V_limits.V_i_min, workspace->V_limits.V_j_max, i - 1);
		else {
			if (i > 0 && i < N_constraints) {
				ind_plane = i - 1;
				index = ind_plane + counter;
				line_slope = coef_ViVj;
				line_temp.offset = intersection_of_2D_lines(plane_normals[index * 3 + 1], plane_normals[index * 3 + 2], 
					plane_normals[(index + 1) * 3 + 1], plane_normals[(index + 1) * 3 + 2], b_planes[index], b_planes[index + 1]);
			}
			else {
				ind_plane = (N_constraints - 1)*(i>0);
				ind_line = i + (N_constraints == 1);
				line_slope = -plane_normals[(ind_plane + counter) * 3] / plane_normals[(ind_plane + counter) * 3 + 1];
				line_temp.offset = (workspace->b_planes[ind_plane + counter] - plane_normals[(ind_plane + counter) * 3 + 2] *
					lines[ind_line].point_end.delta) / plane_normals[(ind_plane + counter) * 3 + 1];
			}
			compute_intersections_of_line_with_Vbox(&line_temp, line_slope, &workspace->V_limits);

			if (fabs(lines[i].slope) > abs_slope)
				errors[i] = error_value_in_ViVj_point(workspace, line_temp.point_end.V_i, line_temp.point_end.V_j, ind_plane);
			else
				errors[i] = error_value_in_ViVj_point(workspace, line_temp.point_begin.V_i, line_temp.point_begin.V_j, ind_plane);
			errors[i] = fmax(errors[i], error_value_in_ViVj_point(workspace, (line_temp.point_begin.V_i+line_temp.point_end.V_i)/2,
				(line_temp.point_begin.V_j+line_temp.point_end.V_j)/2, ind_plane));
		}
	}
}


/*compute one coordinate of the point that lies on the intersection of two 2D lines*/
double intersection_of_2D_lines(double a11, double a12, double a21, double a22, double b1, double b2) {
	return (b1 * a22 - b2 * a12) / (a11*a22 - a21*a12);
}


/*compute value of the approximation error at given (V_i, V_j) point*/
double error_value_in_ViVj_point(LF_Workspace* workspace, double V_i, double V_j, int index_constraint) {
	double delta_line, error_max, offset, *branch_data=workspace->branch_data, *plane_normals = workspace->plane_normals;
	int index = workspace->constraint_counter + index_constraint;

	if (!workspace->all_points_reachable) {
		if (fabs(max_angle_for_given_current_slow(V_i, V_j, branch_data)) > fabs(workspace->delta_max)) {
			//find offset of the line || to intersection line on which this point lies
			offset = V_j - branch_data[8] * V_i;
			//then only part of the line should be used for determining the error
			V_i = update_end_of_intersection_line_segment(workspace, offset);
			V_j = branch_data[8] * V_i + offset;
		}
	}
	delta_line = (workspace->b_planes[index] - plane_normals[index * 3] * V_i -
		plane_normals[index * 3 + 1] * V_j) / plane_normals[index * 3 + 2];
	if (current_magnitude_squared(V_i, V_j, delta_line, branch_data) < 0.0)
		error_max = 100.0;
	else
		error_max = fabs(workspace->I_max_user - current_magnitude(V_i, V_j, delta_line, branch_data)) * 100.0 / workspace->I_max_user;

	return error_max;
}


/*compute minimum error among all constructed relaxed constraints*/
double min_error_relaxed_approximation(double *errors, int N_constraints) {
	double error_min = 1e6, error_temp;
	for (int i = 0; i < N_constraints; i++) {
		error_temp = fmax(errors[i], errors[i + 1]);
		if (error_temp < error_min)
			error_min = error_temp;
	}
	return error_min;
}


/*check if approximation for lower part can be obtained by reflecting approximation for upper part*/
BOOL is_approximation_symmetric(LF_Workspace* workspace) {
	if (fabs(workspace->branch_data[2] / workspace->branch_data[3]) < 0.001)
		return LF_TRUE;
	else
		return LF_FALSE;
}


/*compute parameters of approximation to lower part by reflecting the approximation to upper part*/
void reflect_approximation(LF_Workspace* workspace, int N_constraints_for_upper_part) {
	int counter_max = workspace->constraint_counter, counter_min = workspace->constraint_counter - N_constraints_for_upper_part;
	double delta_upper, delta_lower, *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes;
	
	if (N_constraints_for_upper_part > 0) {
		//compute offset difference
		delta_upper = max_angle_for_given_current_slow(workspace->V_limits.V_i_max, workspace->V_limits.V_i_max, workspace->branch_data);
		workspace->branch_data[9] = -1.0;
		delta_lower = max_angle_for_given_current_slow(workspace->V_limits.V_i_max, workspace->V_limits.V_i_max, workspace->branch_data);
		for (int i = 0; i < N_constraints_for_upper_part; i++) {
			plane_normals[(i + counter_max) * 3] = plane_normals[(i + counter_min) * 3];
			plane_normals[(i + counter_max) * 3 + 1] = plane_normals[(i + counter_min) * 3 + 1];
			plane_normals[(i + counter_max) * 3 + 2] = -plane_normals[(i + counter_min) * 3 + 2];
			b_planes[i + counter_max] = b_planes[i + counter_min] - (delta_lower + delta_upper);
		}
	}
	//update the counter value
	workspace->constraint_counter = counter_max + N_constraints_for_upper_part;
}


/*record results taking into account the type of line*/
void record_all_constraints(LF_Workspace* workspace, LF_Branch* branch, LF_Results* results) {
	int N_constructed_constraints = workspace->constraint_counter;
	double *plane_normals = workspace->plane_normals, *b_planes = workspace->b_planes, *A_matrix = results->A_matrix,
		*b_vector = results->b_vector, Kt_shift = branch->t_shift*LF_PI / 180.0, ratio = transformer_ratio(branch);

	//record results in the output structure
	for (int i = 0; i < N_constructed_constraints; i++) {
		//record constraint normal
		A_matrix[i * 3] = plane_normals[i * 3] / ratio;
		A_matrix[i * 3 + 1] = plane_normals[i * 3 + 1];
		A_matrix[i * 3 + 2] = plane_normals[i * 3 + 2];
		//record constraint offset
		if (A_matrix[i * 3 + 2]>0) //constraint approximates upper part of the surface
			b_vector[i] = b_planes[i] + Kt_shift;
		else //constraint approximates lower part of the surface
			b_vector[i] = b_planes[i] - Kt_shift;
	}

	update_constraints_for_fixed_V(results, branch, N_constructed_constraints);

	results->max_error_est = workspace->max_error_estimate;
	results->N_created_constraints = N_constructed_constraints;
	if (N_constructed_constraints > 0) {
		results->flag_result = success;
		strcpy(results->message, "Linear approximation successfully constructed.");
	}
	else {
		results->flag_result = non_binding;
		strcpy(results->message, "Nonlinear constraint is unreachable, no approximation needed.");
	}
}


/*update the values of constraint normals and offsets in case V_i or V_j is fixed*/
void update_constraints_for_fixed_V(LF_Results *results, LF_Branch *branch, int N_constructed_constraints) {
	double V_temp;
	//update values of coefficients in case V_i is fixed
	if (branch->V_i_min==branch->V_i_max) {
		V_temp = branch->V_i_min;
		for (int i = 0; i < N_constructed_constraints; i++) {
			results->b_vector[i] = results->b_vector[i] - results->A_matrix[i * 3] * V_temp;
			results->A_matrix[i * 3] = 0.0;
		}
	}

	//update values of coefficients in case V_i is fixed
	if (branch->V_j_min == branch->V_j_max) {
		V_temp = branch->V_j_min;
		for (int i = 0; i < N_constructed_constraints; i++) {
			results->b_vector[i] = results->b_vector[i] - results->A_matrix[i * 3 + 1] * V_temp;
			results->A_matrix[i * 3 + 1] = 0.0;
		}
	}
}


/*record relevant parameters of approximation for limits in both the beginning and end of the line*/
void record_relevant_constraints(LF_Workspace *workspace, LF_Branch *branch, 
	int N_constraints_end, LF_Results *results) {
	IntersectionLine line_lower = { { 0.0 } }, line_upper = { { 0.0 } };
	double 	Kt_shift = branch->t_shift*LF_PI / 180.0;
	int indexes_upper[4], indexes_lower[4], N_constructed_constraints=workspace->constraint_counter;

	compute_surface_intersection_lines(workspace, &line_lower, &line_upper);

	//check if intersection lines were successfully constructed
	if (line_lower.point_begin.V_i == 0.0 || line_upper.point_begin.V_i == 0.0 ||
		line_lower.point_end.V_i == 0.0 || line_upper.point_end.V_i == 0.0) {
		record_all_constraints(workspace, branch, results);
		return;
	}

	compute_boundary_indices_of_constraints(workspace, N_constraints_end, indexes_lower, indexes_upper);

	workspace->constraint_counter = 0;

	//record results for upper part of the surface
	record_desired_approximation_part(workspace, results, &line_upper, indexes_upper, Kt_shift);

	//record results for lower part of the surface
	record_desired_approximation_part(workspace, results, &line_lower, indexes_lower, -Kt_shift);

	//if no constraints were recorded, treat it is an error and record all constraints
	if (workspace->constraint_counter == 0) {
		workspace->constraint_counter = N_constructed_constraints;
		record_all_constraints(workspace, branch, results);
		return;
	}

	N_constructed_constraints = workspace->constraint_counter;
	update_constraints_for_fixed_V(results, branch, N_constructed_constraints);

	results->max_error_est = workspace->max_error_estimate;
	results->N_created_constraints = N_constructed_constraints;
	results->flag_result = success;
	strcpy(results->message, "Linear approximation successfully constructed.");
}


/*compute points at the beginning and end of lower and upper intersection lines of two surfaces*/
void compute_surface_intersection_lines(LF_Workspace *workspace, IntersectionLine *line_lower, IntersectionLine *line_upper) {
	double coefficients[4], coef_temp, *branch_data = workspace->branch_data;
	//precompute several coefficients
	coef_temp = branch_data[0] * branch_data[0] + branch_data[1] * branch_data[1] +
		2 * (branch_data[2] * branch_data[2] - branch_data[3] * branch_data[3]);
	coefficients[0] = branch_data[0] * branch_data[1] * coef_temp;
	coefficients[1] = coef_temp*(branch_data[2] * branch_data[2] - branch_data[3] * branch_data[3]);
	coefficients[2] = -2 * (branch_data[0] + branch_data[1])*branch_data[2] * branch_data[2] * branch_data[4];
	coefficients[3] = 4 * branch_data[2] * branch_data[2] * branch_data[4] * branch_data[4];

	//deal with points close to (V_i_min, V_j_min) corner
	compute_surface_intersection_points(workspace, line_lower, line_upper, coefficients, LF_TRUE);
	//deal with points close to (V_i_max, V_j_max) corner
	compute_surface_intersection_points(workspace, line_lower, line_upper, coefficients, LF_FALSE);
}


/*compute points at the intersection of two surfaces with the given edge of Vbox*/
void compute_surface_intersection_points(LF_Workspace *workspace, IntersectionLine *line_lower, IntersectionLine *line_upper,
	double *coefficients, BOOL line_beginning) {
	double V_fixed, V_values[2];
	Point3D points[2] = { 0 };
	int index_lower, index_upper, counter_points = 0;

	//fix V_i
	if (line_beginning)
		V_fixed = workspace->V_limits.V_i_min;
	else
		V_fixed = workspace->V_limits.V_i_max;
	compute_V_at_surfaces_intersection(coefficients, V_fixed*V_fixed, V_values);

	//record the intersection points if they lie within the box
	if (workspace->V_limits.V_j_min <= V_values[0] && V_values[0] <= workspace->V_limits.V_j_max) {
		points[counter_points].V_i = V_fixed;
		points[counter_points].V_j = V_values[0];
		counter_points++;
	}
	if (workspace->V_limits.V_j_min <= V_values[1] && V_values[1] <= workspace->V_limits.V_j_max) {
		points[counter_points].V_i = V_fixed;
		points[counter_points].V_j = V_values[1];
		counter_points++;
	}

	if (counter_points < 2) { //if not both points lie on the intersection with fixed V_i
		//fix V_j
		if (line_beginning)
			V_fixed = workspace->V_limits.V_j_min;
		else
			V_fixed = workspace->V_limits.V_j_max;
		compute_V_at_surfaces_intersection(coefficients, V_fixed*V_fixed, V_values);

		//record the intersection points if they lie within the box
		if (workspace->V_limits.V_i_min <= V_values[0] && V_values[0] <= workspace->V_limits.V_i_max) {
			points[counter_points].V_i = V_values[0];
			points[counter_points].V_j = V_fixed;
			counter_points++;
		}
		if (counter_points < 2 && workspace->V_limits.V_i_min <= V_values[1] && V_values[1] <= workspace->V_limits.V_i_max) {
			points[counter_points].V_i = V_values[1];
			points[counter_points].V_j = V_fixed;
			counter_points++;
		}
	}

	//check if two points were obtained
	if (counter_points < 2)
		return;

	//determine which point relates to the upper intersection line and which to the lower and record it
	compute_deltas_at_surfaces_intersection(workspace->branch_data, &points[0]);
	compute_deltas_at_surfaces_intersection(workspace->branch_data, &points[1]);

	if (points[0].delta == 5.0 || points[1].delta == 5.0) //if there was an error in obtaining the points
		return;

	if (points[0].delta < points[1].delta) {
		index_lower = 0;
		index_upper = 1;
	}
	else {
		index_lower = 1;
		index_upper = 0;
	}

	if (line_beginning) {
		line_lower->point_begin = points[index_lower];
		line_upper->point_begin = points[index_upper];
	}
	else {
		line_lower->point_end = points[index_lower];
		line_upper->point_end = points[index_upper];
	}
}


/*compute values of V at the intersection of two surfaces with the given edge of Vbox*/
void compute_V_at_surfaces_intersection(double *coefficients, double V_fixed, double *V_out) {
	double c1, c2, D;
	//set V values to zero initially
	V_out[0] = 0.0;
	V_out[1] = 0.0;

	//tolve biquadratic equation to get values of V
	c1 = coefficients[1] * V_fixed + coefficients[2];
	c2 = coefficients[0] * V_fixed*V_fixed + 2 * coefficients[2] * V_fixed + coefficients[3];
	D = c1 *c1 - coefficients[0] * c2;
	if (D < 0.0) //exit if there is no solution
		return;
	D = sqrt(D);
	if ((-c1 - D) / coefficients[0] >= 0.0)
		V_out[0] = sqrt((-c1 - D) / coefficients[0]);
	if ((-c1 + D) / coefficients[0] >= 0.0)
		V_out[1] = sqrt((-c1 + D) / coefficients[0]);
}


/*compute values of delta at the given point of the intersection of two surfaces*/
void compute_deltas_at_surfaces_intersection(double *branch_data, Point3D *point) {
	double deltas_point[4];
	//elements of deltas_point are: delta_upper limit at line begin, delta_upper limit at line end, delta_lower at line begin, delta_lower at line end

	//deltas for the surface corresponding to the limit at the beginning of the line
	branch_data[9] = 1.0;
	deltas_point[0] = max_angle_for_given_current(point->V_i, point->V_j, branch_data);
	branch_data[9] = -1.0;
	deltas_point[2] = max_angle_for_given_current(point->V_i, point->V_j, branch_data);

	//deltas for the surface corresponding to the limit at the end of the line
	change_branch_data_for_flow_side(branch_data);
	branch_data[9] = 1.0;
	deltas_point[1] = max_angle_for_given_current(point->V_i, point->V_j, branch_data);
	branch_data[9] = -1.0;
	deltas_point[3] = max_angle_for_given_current(point->V_i, point->V_j, branch_data);
	change_branch_data_for_flow_side(branch_data);

	if (fabs(deltas_point[0] - deltas_point[1]) < 5e-4 || fabs(deltas_point[2] - deltas_point[3]) < 5e-4)
		if (fabs(deltas_point[0] - deltas_point[1]) < fabs(deltas_point[2] - deltas_point[3]))
			point->delta = deltas_point[0];
		else
			point->delta = deltas_point[2];
	else
		point->delta = 5.0;
}


/*change parameters of branch data used for computing surface angle from the limit at one line end to the limit at another line end*/
void change_branch_data_for_flow_side(double *branch_data) {
	double temp = branch_data[0];
	branch_data[0] = branch_data[1];
	branch_data[1] = temp;
	branch_data[6] = -branch_data[6];
}


/*compute first and last indexes of constraints corresponding to lower and upper part of the approximation*/
void compute_boundary_indices_of_constraints(LF_Workspace *workspace, int N_constraints_end, int *indexes_lower, int *indexes_upper) {
	int N_constructed_constraints = workspace->constraint_counter, index_1, index_2;

	//check which planes should be retained from the beginning and end of line limit based on surface intersection
	if (workspace->branch_data[8] < 1) {
		index_1 = 0;
		index_2 = 2;
	}
	else {
		index_1 = 2;
		index_2 = 0;
	}

	indexes_upper[index_1] = 0;
	indexes_lower[index_1 + 1] = N_constraints_end;
	indexes_upper[index_2] = N_constraints_end;
	indexes_lower[index_2 + 1] = N_constructed_constraints;

	for (int i = 0; i < N_constraints_end; i++) {
		if (workspace->plane_normals[3 * i + 2] < 0) {
			indexes_upper[index_1 + 1] = i;
			indexes_lower[index_1] = i;
			break;
		}
	}

	for (int i = N_constraints_end; i < N_constructed_constraints; i++) {
		if (workspace->plane_normals[3 * i + 2] < 0) {
			indexes_upper[index_2 + 1] = i;
			indexes_lower[index_2] = i;
			break;
		}
	}
}


/*record upper or lower part of the approximation into the output structure*/
void record_desired_approximation_part(LF_Workspace *workspace, LF_Results *results, IntersectionLine *line, 
	int *indexes, double Kt_shift) {
	IntersectionLine line_temp;
	//record constraints before the intersection of surfaces
	for (int i = indexes[0]; i < indexes[1]; i++) {
		if (i == indexes[0] || i == indexes[1] - 1) {
			record_parameters_of_one_constraint(workspace, results, i, Kt_shift);
			continue;
		}
		obtain_actual_intersection_line(workspace, i, i + 1, &line_temp);
		if (constraint_should_be_recorded(&line_temp, line, 1.0))
			record_parameters_of_one_constraint(workspace, results, i, Kt_shift);
		else {
			record_parameters_of_one_constraint(workspace, results, i, Kt_shift);
			break;
		}
	}
	//record constraints after the intersection of surfaces
	for (int i = indexes[2]; i < indexes[3]; i++) {
		if (i == indexes[3] - 1) {
			record_parameters_of_one_constraint(workspace, results, i, Kt_shift);
			break;
		}
		obtain_actual_intersection_line(workspace, i, i + 1, &line_temp);
		if (constraint_should_be_recorded(&line_temp, line, -1.0))
			record_parameters_of_one_constraint(workspace, results, i, Kt_shift);
	}
}


/*record parameters of a given constraint into the desired place in the output array*/
void record_parameters_of_one_constraint(LF_Workspace *workspace, LF_Results *results, int index_in_workspace, double Kt_shift) {
	//record constraint normal
	results->A_matrix[workspace->constraint_counter * 3] = workspace->plane_normals[index_in_workspace * 3];
	results->A_matrix[workspace->constraint_counter * 3 + 1] = workspace->plane_normals[index_in_workspace * 3 + 1];
	results->A_matrix[workspace->constraint_counter * 3 + 2] = workspace->plane_normals[index_in_workspace * 3 + 2];
	//record constraint offset
	results->b_vector[workspace->constraint_counter] = workspace->b_planes[index_in_workspace] + Kt_shift;
	//update iteration counter
	workspace->constraint_counter++;
}


/*check if the constraint should be recorded based on the intersection line between this and next constraint and intersection line of two surfaces*/
BOOL constraint_should_be_recorded(IntersectionLine *plane_line, IntersectionLine *surface_line, double sign) {
	if (sign*plane_line->point_begin.V_j / plane_line->point_begin.V_i < sign*surface_line->point_begin.V_j / surface_line->point_begin.V_i ||
		sign*plane_line->point_end.V_j / plane_line->point_end.V_i < sign*surface_line->point_end.V_j / surface_line->point_end.V_i)
		return LF_TRUE;
	else
		return LF_FALSE;
}



/*LOW LEVEL MATH FUNCTIONS*/
/*bisection algorithm (general realization)*/
double compute_value_by_bisection(double x_min, double x_max, double eps_tolerance, double *branch_data,
	BOOL(*retain_first_half)(double, double*, void*), void* data) {
	BOOL split_interval = LF_TRUE;
	double x;

	while (split_interval) {
		if (x_max - x_min < eps_tolerance)
			split_interval = LF_FALSE;
		else {
			x = (x_min + x_max) / 2;
			if (retain_first_half(x, branch_data, data))
				x_max = x;
			else
				x_min = x;
		}
	}
	return (x_min + x_max) / 2;
}


/*check if the first half of the interval has to be retained while looking for the inflection point*/
BOOL bisection_check_for_inflection_point(double V, double *branch_data, void* data) {
	LineSegmentEdge* region = (LineSegmentEdge*)data;
	double delta_curve, delta_line, line_slope = region->line_slope, sign = branch_data[9];
	BOOL retain_first_half;

	//compute the value of delta on the curve
	if (region->tangent_line_type == fixed_V_i)
		delta_curve = max_angle_for_given_current(region->V_fixed, V, branch_data);
	else
		delta_curve = max_angle_for_given_current(V, region->V_fixed, branch_data);

	//compute the value of delta on the line
	delta_line = line_slope*V + region->line_offset;

	if (sign*delta_line>sign*delta_curve) { //if at this point initial approximation is not conservative
		if (sign*line_slope>sign*region->d_delta_V_min)
			retain_first_half = LF_FALSE;
		else
			retain_first_half = LF_TRUE;
	}
	else { //if at this point initial approximation is conservative
		if (sign*line_slope>sign*region->d_delta_V_min)
			retain_first_half = LF_TRUE;
		else
			retain_first_half = LF_FALSE;
	}
	return retain_first_half;
}


/*check if the first half of the interval has to be retained while looking for the surface point that has the same slope as given intersection line*/
BOOL bisection_check_for_intersection_line(double V_i, double *branch_data, void* data) {
	double* temp = (double*)data, sign = branch_data[9];
	BOOL retain_first_half;
	double V_j, d_delta, offset = temp[0], line_slope = temp[1];

	V_j = branch_data[8] * V_i + offset;
	d_delta = slope_of_tangent_line(V_i, V_j, branch_data, intersection_line);

	if (sign*d_delta < sign*line_slope)
		retain_first_half = LF_FALSE;
	else
		retain_first_half = LF_TRUE;

	return retain_first_half;
}


/*check if the first half of the interval has to be retained while looking for the surface point that has the same slope as given line with fixed V*/
BOOL bisection_check_for_line_with_fixed_V(double V, double *branch_data, void* data) {
	LineSegmentEdge* region = (LineSegmentEdge*)data;
	BOOL retain_first_half;
	double d_delta, V_i, V_j, sign = branch_data[9];

	if (region->tangent_line_type == fixed_V_i) {
		V_i = region->V_fixed;
		V_j = V;
	}
	else {
		V_i = V;
		V_j = region->V_fixed;
	}

	d_delta = slope_of_tangent_line(V_i, V_j, branch_data, region->tangent_line_type);

	if (!region->convex) {
		if (sign*d_delta < sign*region->line_slope)
			retain_first_half = LF_FALSE;
		else
			retain_first_half = LF_TRUE;
	}
	else {
		if (sign*d_delta > sign*region->line_slope)
			retain_first_half = LF_FALSE;
		else
			retain_first_half = LF_TRUE;
	}

	return retain_first_half;
}


/*check if the first half of the interval has to be retained while looking for the intersection of plane and surface*/
BOOL bisection_intersection_of_plane_and_surface(double V_i, double *branch_data, void* data_in) {
	double V_j, delta_surface, delta_plane, *data = (double*)data_in;

	V_j = data[0] * V_i + data[1];
	delta_surface = max_angle_for_given_current(V_i, V_j, branch_data);
	delta_plane = -V_i*data[2] - data[3];

	if (data[4] > 0)
		if (delta_surface*branch_data[9] > delta_plane*branch_data[9])
			return LF_TRUE;
		else
			return LF_FALSE;
	else
		if (delta_surface*branch_data[9] > delta_plane*branch_data[9])
			return LF_FALSE;
		else
			return LF_TRUE;
}


/*solve a one-variable nonlinear equation with a Newton-Raphson method*/
double solve_equation_by_Newton_Raphson(double x, LF_Options *options, double *branch_data,
	double(*update_in_Newton_Raphson)(double, double*, void*), void* data) {
	int iter = 0;
	double delta_x;

	while (LF_TRUE) {
		delta_x = update_in_Newton_Raphson(x, branch_data, data);
		x = x + delta_x;
		iter++;

		//check stopping criteria
		if (fabs(delta_x) < options->eps_tolerance)
			break;
		else if (delta_x == 0.0 || iter == options->iter_Max) {
			x = 0;
			break;
		}
	}
	return x;
}


/*compute variable update for Newton Raphson for line with fixed V*/
double update_in_NR_line_with_fixed_V(double V, double *branch_data, void* data_in) {
	LineSegmentEdge *segment = (LineSegmentEdge*)data_in;
	double V_i, V_j, D, D_root, w, w1, w1_root, dw, ddw, t, dt, ddt;

	if (segment->tangent_line_type == fixed_V_i) {
		V_i = segment->V_fixed;
		V_j = V;
		dt = -branch_data[0] * V_i / (V_j*V_j) + branch_data[1] / V_i + branch_data[4] / (V_i*V_j*V_j);
		ddt = 2 * branch_data[0] * V_i / (V_j*V_j*V_j) - 2 * branch_data[4] / (V_i*V_j*V_j*V_j);
	}
	else {
		V_i = V;
		V_j = segment->V_fixed;
		dt = branch_data[0] / V_j - branch_data[1] * V_j / (V_i*V_i) + branch_data[4] / (V_i*V_i*V_j);
		ddt = 2 * branch_data[1] * V_j / (V_i*V_i*V_i) - 2 * branch_data[4] / (V_i*V_i*V_i*V_j);
	}

	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	D = branch_data[5] - t*t;

	if (D > 0.0)
		D_root = sqrt(D);
	else
		return 0.0;

	w = branch_data[6] * t + branch_data[9] * branch_data[7] * D_root;
	dw = (branch_data[6] - branch_data[9] * branch_data[7] * t / D_root)*dt;
	w1 = 1.0 - w *w;
	if (w1 >= 0.0) {
		w1_root = sqrt(w1);
		ddw = branch_data[6] * ddt - branch_data[9] * branch_data[7] * ((dt*dt + t*ddt)*D + t*t*dt*dt) / (D*D_root);
		return (segment->line_slope*w1_root - dw)*w1 / (ddw*w1 + w*dw*dw);
	}
	else
		return 0.0;
}


/*compute variable update for Newton Raphson for interseciton line*/
double update_in_NR_intersection_line(double V_i, double *branch_data, void* data_in) {
	double V_j, D, D_root, w, w1, w1_root, dw, ddw, t, dt, ddt, *data = (double*)data_in;

	V_j = data[2] * V_i + data[0];
	w = 2.0 * data[2] * V_i + data[0];
	dt = branch_data[0] * data[0] / (V_j*V_j) - branch_data[1] * data[0] /
		(V_i*V_i) + branch_data[4] * w / (V_i*V_i*V_j*V_j);
	ddt = -2.0* data[2] * data[0] * branch_data[0] / (V_j*V_j*V_j) + 2 * data[0] * branch_data[1] /
		(V_i*V_i*V_i) + 2 * branch_data[4] * (data[2] * V_i*V_j - w*w) / pow(V_i*V_j, 3);

	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	D = branch_data[5] - t*t;

	if (D > 0.0)
		D_root = sqrt(D);
	else
		return 0.0;

	w = branch_data[6] * t + branch_data[9] * branch_data[7] * D_root;
	dw = (branch_data[6] - branch_data[9] * branch_data[7] * t / D_root)*dt;
	w1 = 1.0 - w *w;
	if (w1 >= 0.0) {
		w1_root = sqrt(w1);
		ddw = branch_data[6] * ddt - branch_data[9] * branch_data[7] * ((dt*dt + t*ddt)*D + t*t*dt*dt) / (D*D_root);
		return (data[1] * w1_root - dw)*w1 / (ddw*w1 + w*dw*dw);
	}
	else
		return 0.0;
}


/*compute variable update for Newton Raphson for interseciton line*/
double update_in_NR_intersection_of_plane_and_surface(double V_i, double *branch_data, void* data_in) {
	double V_j, D, w, w1, dw, t, dt, f, df, *data = (double*)data_in;

	V_j = data[0] * V_i + data[1];
	dt = branch_data[0] * data[1] / (V_j*V_j) - branch_data[1] * data[1] /
		(V_i*V_i) + branch_data[4] * (2.0 * data[0] * V_i + data[1]) / (V_i*V_i*V_j*V_j);
	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	D = branch_data[5] - t*t;

	if (D > 0.0)
		D = sqrt(D);
	else
		return 0.0;

	w = branch_data[6] * t + branch_data[9] * branch_data[7] * D;
	dw = (branch_data[6] - branch_data[9] * branch_data[7] * t / D)*dt;
	w1 = 1.0 - w *w;
	if (w1 >= 0.0) {
		f = max_angle_for_given_current(V_i, V_j, branch_data) + V_i*data[2] + data[3];
		df = dw / sqrt(w1) + data[2];
		return -f / df;
	}
	else
		return 0.0;
}


/*compute the slope of the tangent line at the given point for a line direction*/
double slope_of_tangent_line(double V_i, double V_j, double *branch_data, TangentLineType line_type) {
	double D, w, w1, t, dt, offset, d_delta = 0.0;

	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	if (line_type == fixed_V_i)
		dt = -branch_data[0] * V_i / (V_j*V_j) + branch_data[1] / V_i + branch_data[4] / (V_i*V_j*V_j);
	else if (line_type == fixed_V_j)
		dt = branch_data[0] / V_j - branch_data[1] * V_j / (V_i*V_i) + branch_data[4] / (V_i*V_i*V_j);
	else {
		offset = V_j - branch_data[8] * V_i;
		dt = branch_data[0] * offset / (V_j*V_j) - branch_data[1] * offset /
			(V_i*V_i) + branch_data[4] * (2.0 * branch_data[8] * V_i + offset) / (V_i*V_i*V_j*V_j);
	}

	D = branch_data[5] - t* t;

	if (D > 0.0)
		D = sqrt(D);
	else
		return 0.0;

	w = (branch_data[6] * t + branch_data[9] * branch_data[7] * D);
	w1 = 1.0 - w *w;
	if (w1 >= 0.0)
		d_delta = (branch_data[6] - branch_data[9] * branch_data[7] * t / D)*dt / sqrt(w1);
	return d_delta;
}



/*compute the updated value of the end of intersection line segment that still has delta <= delta_max_user*/
double update_end_of_intersection_line_segment(LF_Workspace* workspace, double offset) {
	double V_i_out, a33, c1, c2, c3, temp_1, delta_max = workspace->delta_max, *branch_data = workspace->branch_data;

	a33 = 2.0 * (branch_data[2] * sin(delta_max) - branch_data[3] * cos(delta_max));
	c1 = 2.0 * branch_data[0] + a33*branch_data[8];
	c2 = 2.0 * offset * (branch_data[1] * branch_data[8] + a33);
	c3 = branch_data[1] * offset * offset - branch_data[4];
	temp_1 = c2 * c2 - 4.0 * c1*c3;
	if (temp_1 < 0)
		temp_1 = 0;
	else
		temp_1 = sqrt(temp_1);
	V_i_out = (-c2 + temp_1) / (2.0 * c1); //only one solution suits us

	return V_i_out;
}


/*compute magnitude of current at given point*/
double current_magnitude(double V_i, double V_j, double delta, double* branch_data) {
	//this function computes the magnitude of a current flowing through the line
	double I_value;
	I_value = (V_i*V_i * branch_data[0] + V_j*V_j * branch_data[1] - 2.0 * V_i*V_j*(branch_data[3] * cos(delta) - branch_data[2] * sin(delta)));
	if (I_value < 0.0)
		I_value = 0;
	else
		I_value = sqrt(I_value);
	return I_value;
}


/*compute square of current magnitude at given point*/
double current_magnitude_squared(double V_i, double V_j, double delta, double* branch_data) {
	return V_i*V_i * branch_data[0] + V_j*V_j * branch_data[1] - 2.0 * V_i*V_j*(branch_data[3] * cos(delta) - branch_data[2] * sin(delta));
}


/*Compute delta at which I=I_max_user (for upper or lower side of boundary surface). Can only handle abs(delta)<pi/2*/
double max_angle_for_given_current(double V_i, double V_j, double* branch_data) {
	double D, t, sin_delta;

	//set some intermediate variables for ease of computation
	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	D = branch_data[5] - t*t;

	//check the determinant and compute the sine of delta
	if (D < 0.0) //based on how we select points, this can happen only due to small numerical error => we assume that D=0
		sin_delta = branch_data[6] * t;
	else
		sin_delta = branch_data[6] * t + branch_data[9] * branch_data[7] * sqrt(D);
	return asin(sin_delta);
}


/*Compute delta at which I=I_max_user (for upper or lower side of boundary surface). Can any value of delta*/
double max_angle_for_given_current_slow(double V_i, double V_j, double* branch_data) {
	double D, t, sin_delta, delta_critical, delta_critical1, w, w1;

	//set some intermediate variables for ease of computation
	t = branch_data[0] * V_i / V_j + branch_data[1] * V_j / V_i - branch_data[4] / (V_i*V_j);
	D = branch_data[5] - t*t;

	//check the determinant and compute the sine of delta
	if (D < 0.0) //based on how we select points, this can happen only due to small numerical error => we assume that D=0
		delta_critical = asin(branch_data[6] * t);
	else {
		//compute the value of the sine of delta
		sin_delta = branch_data[6] * t + branch_data[9] * branch_data[7] * sqrt(D);

		//compute the angle
		delta_critical = asin(sin_delta);
		delta_critical1 = branch_data[9] * LF_PI - delta_critical;

		//check if the equation is satisfied
		w = fabs(branch_data[0] * V_i*V_i + branch_data[1] * V_j*V_j + 2.0 * V_i*V_j*(branch_data[2] * sin(delta_critical) - 
			branch_data[3] * cos(delta_critical)) - branch_data[4]);
		w1 = fabs(branch_data[0] * V_i*V_i + branch_data[1] * V_j*V_j + 2.0 * V_i*V_j*(branch_data[2] * sin(delta_critical1) - 
			branch_data[3] * cos(delta_critical1)) - branch_data[4]);
		if (w > w1)
			delta_critical = delta_critical1;
	}
	return delta_critical;
}


/*record value of transformer tap ratio (1 if no transformer)*/
double transformer_ratio(LF_Branch* branch) {
	if (branch->t_ratio == 0.0)
		return 1.0;
	else
		return branch->t_ratio;
}


/*compute vector product of two vectors*/
void compute_vector_product(double* input_1, double* input_2, double* output) {
	output[0] = input_1[1] * input_2[2] - input_1[2] * input_2[1];
	output[1] = input_1[2] * input_2[0] - input_1[0] * input_2[2];
	output[2] = input_1[0] * input_2[1] - input_1[1] * input_2[0];
}


/*fill 3-element vector with given values*/
void fill_vector_with_values(double value_1, double value_2, double value_3, double* output_vector) {
	output_vector[0] = value_1;
	output_vector[1] = value_2;
	output_vector[2] = value_3;
}


/*compute number of intersection lines for a given number of constraints*/
int max_number_of_intersection_lines(int N_constraints) {
	if (N_constraints == 1)
		return 3; //the algorithm is tricked into producing an extra intersection line in order to compute the slope in (V_i, delta) plane
	else
		return N_constraints + 1;
}


/*find value of maximum element in 1D array*/
double max_element_in_array(double* Array, int Length) {
	double max_value;
	max_value = Array[0];
	for (int i = 1; i < Length; i++)
		if (Array[i] > max_value)
			max_value = Array[i];
	return max_value;
}


/*find value of minimum element in 1D array*/
double min_element_in_array(double* Array, int Length) {
	double min_value;
	min_value = Array[0];
	for (int i = 1; i < Length; i++)
		if (Array[i] < min_value)
			min_value = Array[i];
	return min_value;
}




/*getter functions*/
double* LF_get_A_matrix(LF_Results* results) {
	return results->A_matrix;
}

double* LF_get_b_vector(LF_Results* results) {
	return results->b_vector;
}

LF_ResultFlag LF_get_flag(LF_Results* results) {
	return results->flag_result;
}

int LF_get_number_constraints(LF_Results* results) {
	return results->N_created_constraints;
}

double LF_get_error(LF_Results* results) {
	return results->max_error_est;
}

char* LF_get_message(LF_Results* results) {
	return results->message;
}

int LF_get_mode(LF_Options* options) {
	return options->computation_mode;
}

int LF_get_iter_max(LF_Options* options) {
	return options->iter_Max;
}

int LF_get_N_constraints_max(LF_Options* options) {
	return options->N_constraints_max;
}

int LF_get_N_adjustments(LF_Options* options) {
	return options->N_adjustments;
}

int LF_get_transformer_model_type(LF_Options* options) {
	return options->tr_model_type;
}

double LF_get_delta_max_user(LF_Options* options) {
	return options->delta_max_user;
}

double LF_get_eps_tolerance(LF_Options* options) {
	return options->eps_tolerance;
}

double LF_get_error_max(LF_Options* options) {
	return options->error_max;
}

double LF_get_max_error_change(LF_Options* options) {
	return options->max_error_change;
}

double LF_get_ratio_threshold(LF_Options* options) {
	return options->ratio_threshold;
}

int LF_get_approximation_type(LF_Options* options) {
	return (int)options->approximation;
}



/*setter functions*/
void LF_set_mode(int computation_mode, LF_Options* options) {
	options->computation_mode = computation_mode;
}

void LF_set_iter_max(int iter_max, LF_Options* options) {
	options->iter_Max = iter_max;
}

void LF_set_N_constraints_max(int N_constraints_max, LF_Options* options) {
	options->N_constraints_max = N_constraints_max;
}

void LF_set_N_adjustments(int N_adjustments, LF_Options* options) {
	options->N_adjustments = N_adjustments;
}

void LF_set_transformer_model_type(int transformer_model_type, LF_Options* options) {
	options->tr_model_type = transformer_model_type;
}

void LF_set_delta_max_user(double delta_max_user, LF_Options* options) {
	options->delta_max_user = delta_max_user*LF_PI / 180;
}

void LF_set_eps_tolerance(double eps_tolerance, LF_Options* options) {
	options->eps_tolerance = eps_tolerance;
}

void LF_set_error_max(double error_max, LF_Options* options) {
	options->error_max = error_max;
}

void LF_set_max_error_change(double max_error_change, LF_Options* options) {
	options->max_error_change = max_error_change;
}

void LF_set_ratio_threshold(double ratio_threshold, LF_Options* options) {
	options->ratio_threshold = ratio_threshold;
}

void LF_set_approximation_type(int approximation, LF_Options* options) {
	if (approximation == 0)
		options->approximation = conservative;
	else
		options->approximation = relaxed;
}

void LF_set_branch_parameters(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
	double g, double b, double b_sh, double t_ratio, double t_shift, double I_max, LF_Branch* branch) {
	branch->V_i_min = V_i_min;
	branch->V_i_max = V_i_max;
	branch->V_j_min = V_j_min;
	branch->V_j_max = V_j_max;
	branch->g = g;
	branch->b = b;
	branch->b_sh = b_sh;
	branch->t_ratio = t_ratio;
	branch->t_shift = t_shift;
	branch->I_max = I_max;
}
