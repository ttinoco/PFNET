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

#ifndef __LF_HEADER__
#define __LF_HEADER__

#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#define LF_PI 3.14159265358979323846


/*CUSTOM TYPES*/
typedef struct LF_Branch { //contains branch parameters
	double V_i_min;
	double V_i_max;
	double V_j_min;
	double V_j_max;
	double g;
	double b;
	double b_sh;
	double t_ratio;
	double t_shift;
	double I_max;
} LF_Branch;

typedef struct LF_Results LF_Results;

typedef struct LF_Options LF_Options;

typedef struct VoltageLimits VoltageLimits;

typedef struct FeasibleRegion FeasibleRegion;

typedef struct Point3D Point3D;

typedef struct IntersectionLine IntersectionLine;

typedef struct LF_Workspace LF_Workspace;

typedef struct LineSegmentEdge LineSegmentEdge;

typedef enum LF_ResultFlag { non_binding, infeasible, success, error_branch_data, 
error_options, zero_limit, error_other } LF_ResultFlag;

typedef enum BranchType { gen_gen, gen_load, load_load } BranchType;

typedef enum TangentLineType { fixed_V_i, fixed_V_j, intersection_line } TangentLineType;

typedef enum ApproximationType { conservative, relaxed } ApproximationType;

typedef enum PlaneType {first, ordinary, corner, last} PlaneType;

typedef char BOOL;
enum { LF_FALSE = 0, LF_TRUE = 1 };


/*Outer interface function. Check input branch branch_data, initialize workspace and output structures,
and construct linear approximation*/
LF_Results* LF_construct(LF_Branch* branch, int flow_side, LF_Options* options);

/*create default algorithm options*/
LF_Options* LF_get_default_options(void);

/*initialize output structure*/
LF_Results* LF_initialize_results(LF_Options* options, int flow_side);

/*free memory taken by output structure*/
void LF_free_results(LF_Results* results);

/*check correctness of algorithm options*/
BOOL LF_check_options(LF_Options* options, LF_Results* results, int flow_side);

/*check correctness of line branch_data*/
BOOL LF_check_branch_data(LF_Branch* branch, LF_Results* results);

/*allocate memory for workspace*/
LF_Workspace* LF_initizalize_workspace(LF_Options* options, int flow_side);

/*free memory taken by workspace*/
void LF_free_workspace(LF_Workspace* workspace);

/*construct linear approximation using provided workspace and results structures*/
void LF_linearize_one_line(LF_Branch* branch, int flow_side, LF_Options* options, LF_Workspace* workspace, LF_Results* results);

/*construct linear approximation for one end of the line*/
void linearize_one_line_end(int flow_side, LF_Options* options, LF_Workspace* workspace, LF_Results* results);

/*compute values derived from branch parameters that will be extensively used throughout the algorithm*/
void compute_branch_parameters(LF_Branch* branch, int flow_side, LF_Options* options, LF_Workspace* workspace);

/*check if intersection of ViVj box with feasible region of line flow constraint is empty*/
BOOL is_feasible_region_empty(LF_Workspace* workspace);

/*extract the information on the feasible region of line flow constraint*/
void compute_parameters_of_feasible_region(LF_Workspace* workspace, LF_Options* options);

/*initialize parameters of feasible region to their default values*/
void set_default_flags_for_feasible_region(FeasibleRegion* feas_reg);

/*compute values of current at box corners and delta=+-delta_max_user-Kt_shift*/
void compute_I_max_at_Vbox_corners(LF_Workspace* workspace);

/*tighten voltage limits to make sure that at least three corners of the box can have I>=I_max_user for delta<=delta_max_user*/
void tighten_V_limits_for_reachability(LF_Workspace* workspace);

/*tighten limit on V_i or V_j to make sure that at this limit I=I_max_user for delta=delta_max_user*/
double tightened_V_limit(LF_Workspace* workspace, double V_value, int flag_V_fixed);

/*check and record whether top-left and bottom-right box corners are feasible*/
void establish_feasibility_of_Vbox_corners(LF_Workspace* workspace);

/*compute min amd max offsets of a line that is || to the ridge and contains at least one feasible point*/
void compute_offsets_of_outermost_lines(LF_Workspace* workspace);

/*tighten limits on V for line between gen and load to make sure all V points within limits are feasible*/
void tighten_V_limits_for_feasibility(LF_Workspace* workspace);

/*check if the feasible region is convex for line between gen and load and if it is non-convex, retain only the part that causes nonconvexity*/
void establish_convexity_of_feasible_region(LF_Workspace* workspace, LF_Options* options);

/*compute specific parameters of feasible region for line between load and generator*/
void compute_parameters_of_edge_line_segment(double *branch_data, VoltageLimits* Vbox, LineSegmentEdge* segment);

/*Construct linear approximation in a nested loop.*/
void compute_parameters_of_approximation(LF_Workspace* workspace, LF_Options* options);

/*coumpute the relative width of the constraint in the middle*/
double compute_width_of_middle_constraint(LF_Workspace* workspace, int N_constructed_constraints);

/*compute the values of offsets of intersection lines given that their number has changed*/
void compute_offsets_of_new_intersection_lines(LF_Workspace* workspace, int N_constructed_constraints, int N_new_constraints);

/*insert new intersection line based on the offsets of already constructed intersection lines*/
void insert_intersection_line(IntersectionLine* lines, double dif_offset, int N_constructed_constraints);

/*remove intersection line based on the offsets of already constructed intersection lines*/
void remove_intersection_line(IntersectionLine* lines, double dif_offset, int N_constructed_constraints);

/*refine the values of offsets of existing intersection lines in order to reduce the approximation error*/
void refine_offsets_of_existing_interseciton_lines(LF_Workspace* workspace, 
	int N_constructed_constraints, ApproximationType approximation);

/*compute the offsets of intersection lines such that the difference between two adjacent offsets is the same for all lines*/
void distribute_offsets_evenly(LF_Workspace* workspace, int N_constructed_constraints, double dif_alpha);

/*determine parameters of linear constraints for a given number of constraints*/
void construct_linear_constraints(LF_Workspace* workspace, LF_Options* options, int N_constraints);

/*construct conservative linear approximation for line between gen and load*/
void construct_approximation_to_nonconvex_region(LF_Workspace* workspace, LF_Options* options);

/*estimate approximation error for line between load and generator in case of non-convex feasible region*/
void estimate_approximation_error_for_nonconvex_region(LF_Workspace* workspace, LF_Options* options);

/*compute all parameters related to intersection lines given their offsets*/
void construct_intersection_lines(LF_Workspace* workspace, LF_Options* options, int N_constraints);

/*compute parameters of points at which a line intersects with Vbox*/
void compute_intersections_of_line_with_Vbox(IntersectionLine* line, double line_slope, VoltageLimits* Vbox);

/*update values of end points of the intersection line to make sure that all points between them can have I>=I_max for delta<=delta_max*/
void shrink_intersection_line_segment(double delta_begin, double delta_end, LF_Workspace* workspace, IntersectionLine* line);

/*compute the slope that all intersection lines will have in (V_i, delta) plane*/
double slope_of_intersection_lines(IntersectionLine* lines, int n_lines_total);

/*for each intersection line, compute the value of delta at either line beginning or end*/
void compute_deltas_at_ends_of_intersection_lines(LF_Workspace* workspace, int n_lines_total, LF_Options* options);

/*compute value of V_i at the point at which the conservative intersection line touches the surface*/
double point_on_conservative_intersection_line(IntersectionLine* line, double *branch_data,
	double slope_all_lines, LF_Options *options);

/*compute value of V_i at the point at which the initial relaxed intersection line touches the surface*/
double point_on_relaxed_intersection_line(IntersectionLine* line, LF_Workspace* workspace);

/*compute plane normals for all planes*/
void compute_plane_normals(LF_Workspace* workspace, int N_constraints, ApproximationType approximation);

/*compute the offsets of planes (conservative or relaxed)*/
void compute_plane_offsets(LF_Workspace* workspace, int N_constraints, LF_Options* options);

/*compute the required shift in delta to ensure that a given plane represents a relaxed approximation to the surface*/
double delta_shift_plane_relaxed(LF_Workspace* workspace, int index, int N_constraints, LF_Options* options);

/*compute the shift in delta based on the edge located at either beginning or end of the box*/
double delta_shift_one_side_of_plane(LF_Workspace* workspace, int* line_indexes, LF_Options* options,
	BOOL bottom_left_corner, PlaneType plane_type);

/*compute parameters of the line segment that represents the intersection of a plane with the edge of the box*/
void compute_parameters_of_plane_edge_segment(double *branch_data, LineSegmentEdge* segment,
	Point3D* point1, Point3D* point2);

/*compute value of delta on the plane in either (V_i_min, V_j_min) or (V_i_max, V_j_max) corner of Vbox*/
void compute_delta_plane_in_Vbox_corner(LF_Workspace* workspace, int index, Point3D* point_line, Point3D* point_corner);

/*compute by how much the plane should be shifted to become relaxed based on a given edge of Vbox*/
double delta_shift_on_edge(LF_Workspace* workspace, LineSegmentEdge* segment, LF_Options* options, PlaneType plane_type);

/*update parameters of plane if a box corner is infeasible to make sure the approximation is conservative*/
void update_plane_for_infeasible_box_corner(LF_Workspace *workspace, int index_outer, int N_constraints, LF_Options *options);

/*compute parameters of a point that lies on the updated plane*/
void obtain_point_on_updated_plane(LF_Workspace *workspace, int index_outer, LF_Options *options, Point3D *point);

/*compute parameters of the point on the intersection line that violates the conservative approximation the most*/
void obtain_most_distant_point(LF_Workspace *workspace, int index_outer, int index_inner, Point3D *point);

/*compute parameters of the actual intersection line of two planes*/
void obtain_actual_intersection_line(LF_Workspace *workspace, int index_outer, int index_inner, IntersectionLine *line);

/*estimate maximum values of conservative approximation errors for each of constructed constraints*/
void estimate_conservative_approximation_errors(LF_Workspace* workspace, int N_constraints);

/*estimate values of relaxed approximation errors for each of constructed constraints*/
void estimate_relaxed_approximation_errors(LF_Workspace* workspace, int N_constraints);

/*compute one coordinate of the point that lies on the intersection of two 2D lines*/
double intersection_of_2D_lines(double a11, double a12, double a21, double a22, double b1, double b2);

/*compute value of the approximation error at given (V_i, V_j) point*/
double error_value_in_ViVj_point(LF_Workspace* workspace, double V_i, double V_j, int index_constraint);

/*compute minimum error among all constructed relaxed constraints*/
double min_error_relaxed_approximation(double *errors, int N_constraints);

/*check if approximation for lower part can be obtained by reflecting approximation for upper part*/
BOOL is_approximation_symmetric(LF_Workspace* workspace);

/*compute parameters of approximation to lower part by reflecting the approximation to upper part*/
void reflect_approximation(LF_Workspace* workspace, int N_constraints_for_upper_part);

/*record results taking into account the type of line*/
void record_all_constraints(LF_Workspace* workspace, LF_Branch* branch, LF_Results* results);

/*update the values of constraint normals and offsets in case V_i or V_j is fixed*/
void update_constraints_for_fixed_V(LF_Results *results, LF_Branch *branch, int N_constructed_constraints);

/*record relevant parameters of approximation for limits in both the beginning and end of the line*/
void record_relevant_constraints(LF_Workspace *workspace, LF_Branch *branch, int N_constraints_end, LF_Results *results);

/*compute points at the beginning and end of lower and upper intersection lines of two surfaces*/
void compute_surface_intersection_lines(LF_Workspace *workspace, IntersectionLine *line_lower, IntersectionLine *line_upper);

/*compute points at the intersection of two surfaces with the given edge of Vbox*/
void compute_surface_intersection_points(LF_Workspace *workspace, IntersectionLine *line_lower, IntersectionLine *line_upper,
	double *coefficients, BOOL line_beginning);

/*compute values of V at the intersection of two surfaces with the given edge of Vbox*/
void compute_V_at_surfaces_intersection(double *coefficients, double V_fixed, double *V_out);

/*compute values of delta at the given point of the intersection of two surfaces*/
void compute_deltas_at_surfaces_intersection(double *branch_data, Point3D *point);

/*change parameters of branch data used for computing surface angle from the limit at one line end to the limit at another line end*/
void change_branch_data_for_flow_side(double *branch_data);

/*compute first and last indexes of constraints corresponding to lower and upper part of the approximation*/
void compute_boundary_indices_of_constraints(LF_Workspace *workspace, int N_constraints_end, int *indexes_lower, int *indexes_upper);

/*record upper or lower part of the approximation into the output structure*/
void record_desired_approximation_part(LF_Workspace *workspace, LF_Results *results, IntersectionLine *line, 
	int *indexes, double Kt_shift);

/*record parameters of a given constraint into the desired place in the output array*/
void record_parameters_of_one_constraint(LF_Workspace *workspace, LF_Results *results, 
	int index_in_workspace, double Kt_shift);

/*check if the constraint should be recorded based on the intersection line between this and next constraint and intersection line of two surfaces*/
BOOL constraint_should_be_recorded(IntersectionLine *plane_line, IntersectionLine *surface_line, double sign);


/*LOW LEVEL MATH FUNCTIONS*/
/*bisection algorithm (general realization)*/
double compute_value_by_bisection(double x_min, double x_max, double eps_tolerance, double *branch_data,
	BOOL(*retain_first_half)(double, double*, void*), void* data);

/*check if the first half of the interval has to be retained while looking for the inflection point*/
BOOL bisection_check_for_inflection_point(double V, double *branch_data, void* data);

/*check if the first half of the interval has to be retained while looking for the surface point that has the same slope as given intersection line*/
BOOL bisection_check_for_intersection_line(double V_i, double *branch_data, void* data);

/*check if the first half of the interval has to be retained while looking for the surface point that has the same slope as given line with fixed V*/
BOOL bisection_check_for_line_with_fixed_V(double V, double *branch_data, void* data);

/*check if the first half of the interval has to be retained while looking for the intersection of plane and surface*/
BOOL bisection_intersection_of_plane_and_surface(double V_i, double *branch_data, void* data_in);

/*solve a one-variable nonlinear equation with a Newton-Raphson method*/
double solve_equation_by_Newton_Raphson(double x, LF_Options *options, double *branch_data,
	double(*update_in_Newton_Raphson)(double, double*, void*), void* data);

/*compute variable update for Newton Raphson for line with fixed V*/
double update_in_NR_line_with_fixed_V(double V, double *branch_data, void* data_in);

/*compute variable update for Newton Raphson for interseciton line*/
double update_in_NR_intersection_line(double V_i, double *branch_data, void* data_in);

/*compute variable update for Newton Raphson for interseciton line*/
double update_in_NR_intersection_of_plane_and_surface(double V_i, double *branch_data, void* data_in);

/*compute the slope of the tangent line at the given point for a line direction*/
double slope_of_tangent_line(double V_i, double V_j, double *branch_data, TangentLineType line_type);

/*compute the updated value of the end of intersection line segment that still has delta <= delta_max_user*/
double update_end_of_intersection_line_segment(LF_Workspace* workspace, double offset);

/*compute magnitude of current at given point*/
double current_magnitude(double V_i, double V_j, double delta, double* branch_data);

/*compute square of current magnitude at given point*/
double current_magnitude_squared(double V_i, double V_j, double delta, double* branch_data);

/*Compute delta at which I=I_max_user (for upper or lower side of boundary surface). Can only handle abs(delta)<pi/2*/
double max_angle_for_given_current(double V_i, double V_j, double* branch_data);

/*Compute delta at which I=I_max_user (for upper or lower side of boundary surface). Can any value of delta*/
double max_angle_for_given_current_slow(double V_i, double V_j, double* branch_data);

/*record value of transformer tap ratio (1 if no transformer)*/
double transformer_ratio(LF_Branch* branch);

/*compute vector product of two vectors*/
void compute_vector_product(double* input_1, double* input_2, double* output);

/*fill 3-element vector with given values*/
void fill_vector_with_values(double value_1, double value_2, double value_3, double* output_vector);

/*compute number of intersection lines for a given number of constraints*/
int max_number_of_intersection_lines(int N_constraints);

/*find value of maximum element in 1D array*/
double max_element_in_array(double* Array, int Length);

/*find value of minimum element in 1D array*/
double min_element_in_array(double* Array, int Length);


/*getter functions*/
double* LF_get_A_matrix(LF_Results* results);

double* LF_get_b_vector(LF_Results* results);

LF_ResultFlag LF_get_flag(LF_Results* results);

int LF_get_number_constraints(LF_Results* results);

double LF_get_error(LF_Results* results);

char* LF_get_message(LF_Results* results);

int LF_get_mode(LF_Options* options);

int LF_get_iter_max(LF_Options* options);

int LF_get_N_constraints_max(LF_Options* options);

int LF_get_N_adjustments(LF_Options* options);

int LF_get_transformer_model_type(LF_Options* options);

double LF_get_delta_max_user(LF_Options* options);

double LF_get_eps_tolerance(LF_Options* options);

double LF_get_error_max(LF_Options* options);

double LF_get_max_error_change(LF_Options* options);

double LF_get_ratio_threshold(LF_Options* options);

int LF_get_approximation_type(LF_Options* options);

/*setter functions*/
void LF_set_mode(int computation_mode, LF_Options* options);

void LF_set_iter_max(int iter_max, LF_Options* options);

void LF_set_N_constraints_max(int N_constraints_max, LF_Options* options);

void LF_set_N_adjustments(int N_adjustments, LF_Options* options);

void LF_set_transformer_model_type(int transformer_model_type, LF_Options* options);

void LF_set_delta_max_user(double delta_max_user, LF_Options* options);

void LF_set_eps_tolerance(double eps_tolerance, LF_Options* options);

void LF_set_error_max(double error_max, LF_Options* options);

void LF_set_max_error_change(double max_error_change, LF_Options* options);

void LF_set_ratio_threshold(double ratio_threshold, LF_Options* options);

void LF_set_approximation_type(int approximation, LF_Options* options);

void LF_set_branch_parameters(double V_i_min, double V_i_max, double V_j_min, double V_j_max,
	double g, double b, double b_sh, double t_ratio, double t_shift, double I_max, LF_Branch* branch);

#endif

