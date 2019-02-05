Unreleased
----------
* Constant current constraint for CSC HVDC.

Version 1.3.4
-------------
* Updated parser init/free routines to make parser params work.
* Modularized heuristics.
* Added constrained function.
* Added problem num_extra_vars bug fix for constrained functions.
* Added "apply" method to heuristic.
* Made heuristic accept array of constraint pointers instead of linked list of constraints.
* Added contingency name and corresponding setter/getter.
* Updated json float sig digits from 11 to 16.
* Modular flow routines for branches.
* Added array attribute getters for network component structs.
* Moved constraint alloc/clear to "parent" struct (advanced custom constraints can extend).
* Made constraint init/free optional.
* Removed NBOUND constraint.
* Changed PVPQ switching heuristic to HEUR_PVPQ_SWITCHING for consistency with constraint.
* Renamed LBOUND constraint to BOUND.
* Moved function alloc/clear to "parent" struct (advanced custom functions can extend).
* Made function init/free optional.
* Changed count/analyze/eval/etc loop to be bus-based.
* Fixed bug with PVPQ switching constraint when no Q vars.
* Voltage dependent loads: comp_cp, comp_cq, comp_ci, comp_cj, comp_cg, and comp_cb attributes, is_vdep routine, LOAD_PROP_VDEP property, updated parsers.
* Added CONSTR_LOAD_VDEP constraint.
* DC buses, DC branches, and VSC converters for HVDC.
* Added shunt types and net counters.
* Added bus property BUS_PROP_VSET_REG and "is_v_set_regulated" method.
* CSC converters for HVDC.
* FACTS devices.
* Extended count/analyze/eval framework (and added safeguards to existing constr/func/heur) to loop through ac and then dc buses. 
* Added regulating object interface.
* Added constraint for VSC DC voltage control (CONSTR_VSC_DC_VSET).
* Added constraint for VSC DC power control (CONSTR_VSC_DC_PSET).
* Added constraint for HVDC power balance (CONSTR_HVDCPF).
* Added constraint for VSC equations (CONSTR_VSC_EQ).
* Updated PVPQ switching constraints to use reg object interface.
* Updated PVPQ switching heuristic to use reg object interface.
* Added function for encouraging VSC DC power control (FUNC_VSC_DC_PSET).
* Changed CONSTR_REG_GEN to CONSTR_REG_VSET and made it use the reg object interface and changed name to "voltage set point regulation".
* Changed bus sens_v_reg_by_gen to sens_v_set_reg.
* Made PVPQ switching heuristic consider reg buses that are slack (to be consistent with constr_REG_VSET).
* Added switching constraint for power factor regulation (CONSTR_REG_PF_SWITCH).
* Added switching heuristic for power factor regulation (HEUR_REG_PF_SWITCH).
* Added smooth constraint for power factor regulation (CONSTR_REG_PF).
* Added switching constraints for FACTS active/reactive power control (CONSTR_FACTS_PSET_SWITCH and CONSTR_FACTS_QSET_SWITCH).
* Added functions for FACTS active/reactive power control (FUNC_FACTS_PSET and FUNC_FACTS_QSET).
* Added constraint for FACTS equations (CONSTR_FACTS_EQ).
* Added load "in_service" field, getter and setter. Not used anywhere.
* Added routines for updating PQ components of load according to provided weights.
* Added very basic CSC constraints and functions (CONSTR_CSC_EQ, CONSTR_CSC_DC_PSET, CONSTR_CSC_DC_VSET, FUNC_CSC_DC_PSET).
* Added JSON support (read/write) for new components (vdep loads, csc, vsc, facts, dc bus, dc branch).
* Added switched shunt control mode (discrete, continuous) and rounding capability.
* Added network routing for rounding susceptance of discrete switched shunts and count.
* Added branch routine for using power flow count/analyze/eval subroutines to construct Jacobian of (P_km, Q_km) or (P_mk, Q_mk).
* Removed graphviz interface and dependency.
* Removed number of actions from network properties.
* Added redundant buses.
* Added output_level to nework component summary output.

Version 1.3.3
-------------
* Improved gen Q participation to look exactly at which Qs are vars to add correct number of constraints.
* Added support for changing function and constraint parameters (CONSTR_set_parameter, FUNC_set_parameter).
* Added "variable regularization" function or FUNC_REG_VAR, which has parameters w and x0 and computes (x-x0)^Tdiag(w)(x-x0).
* Symmetric connectors/removers for all bus-connected components (connecting A to B also connects B to A).
* Exposed all obj.set_bus and bus.add/remove_obj routines in Python, and made "obj.bus = None" disconnect obj from bus.
* Renamed add_batteries/add_var_generators to add_battery_from_parameters/add_var_generators_from_parameters.
* Added net routines for adding and removing generators from the network.
* Added net routines for adding and removing loads from the network.
* Added net routines for adding and removing shunts from the network.
* Added net routines for adding and removing branches from the network.
* Added net routines for adding and removing buses from the network.
* Added net routines for adding and removing batteries from the network.
* Added net routines for adding and removing var generators from the network.
* Added net routine for extracting subnetwork containing a specific set of buses.
* Extended problem.show() to show number of vars and constraints of each type.
* Fix bug involving problem/constraint/function's network going out of scope in Python.
* Added load reactive power limits.
* Integrated line_flow library source and header for linearized AC branch flow limits.
* Added Q_par (reactive power participation) field to generator.
* Changed REG_GEN constraint to treat generators separately so that constraint can be used without participations.
* Changed PAR_GEN_Q constraint to PVPQ_SWITCHING, which enforces flexible participations based on Q_par and performs all required modifications for PV-PQ switching heuristics.
* Updated PVPQ switching heuristics to utilize PVPQ_SWITCHING constraint.
* Removed net adjust_generators.
* Improved handling of branch and gen outages: Outages can be enabled by setting the outage flag. Contingencies no longer disconnect components.
* Added bus.is_star(), branch.is_part_of_3_winding_transformer(), and net.get_num_star_buses().
* Added routines for getting number of gens/branches on outage.
* Added utility routines to check function gradients and Hessians in pfnet.tests.utils.
* Added utility routines to check constraint Jacobians and Hessians in pfnet.tests.utils.
* Created network state tag to track net changes and make constr/funcs robust to net changes like outage flag updates.
* Updated bus.get_index_P and bus.get_index_Q to make bus.get_index_P work for both AC and DC PF as well as for multiple periods.
* Added bus.get_index_t for getting unique indices for each bus and time.
* Added bus area and zone.
* Started julia wrapper.
* Moved python and julia wrappers to separate repos.

Version 1.3.2
-------------
* Moved to setuptools.
* Changed configure.ac to look for raw_parser using RAW_PARSER env variable.
* Removed C library tarball from python/docs/_static.
* Made Python wrapper self-contained, i.e., it contains and installs the C library from a tarball in lib directory.
* Updated travis instructions to directly install python wrapper.
* Removed "build" script for readthedocs since it was no longer needed.
* Removed setup.cfg since it was no longer needed.
* Added property mask to projection operator.
* Added projection getters for network and extra variables of constraints.
* Made generator Qmin, Qmax writable in Python.
* Added routines for extracting constraint row info strings with format (constr_name:obj_type:obj_index:constr_info:time).
* Made AC_FLOW_LIM constraint store row info.
* Made LBOUND constraint store row info.
* Updated LOAD_PF (load constant power factor) constraint to maintain sign differences of current P and Q.
* Made load "set_target_power_factor" routine accept negative power factors.
* Fixed memory leaks coming from "VEC_new_from_array" in Python wrapper.
* Added routines for extracting info strings about entries of var values vector (obj_type:obj_index:quantity:time).
* Added "v_base" field to bus structure/object to store base voltage in kv and updated all parsers to store this info.
* Made Network object pickle-able.
* Added automatic enforcement of lower-triangularity of Hessian of objective functions.
* Made Contingency object pickle-able.
* Added test utilities in new module pfnet.tests.
* Exposed component flags bit masks in Python for network comparisons.
* Added network "get_copy" and "copy_from_net" routines and enhanced network comparison test utility.
* Made "get_index"-type routines of network components return -1 for NULL pointers to avoid silent errors.
* Changed Contigency class methods to use word "generator" instead of "gen" abbreviation.
* Added network getters for branch current and apparent power magnitudes.
* Added cmake windows build script invokation in setup.py.
* Changed autotools and cmake rules to incorporate raw_parser using conditional joint compilation of sources as opposed to linkage with external library.
* Changed slack limits in branch flow constraints from [0,thermal_rating] to [-thermal_rating,thermal_rating].
* Improved contingency to disable slack bus if all its generators are on outage.
* Added name attributes to all bus-connected components and removed vargens name hash.
* Added branch "get_rating" method that takes 'A', 'B', or 'C' as argument.
* Eliminated CustomParser and dummy python parser example.
* Fixed bug with treatment of outage branches in PV-PQ switching heuristics.
* Added arrays branch_outages and generator_outages to contingency object to store indices of outage components.
* Added unittests to check that examples run without errors.
* Extended support for storing and retrieving constraint sensitivity information.
* Updated Python wrapper documentation to show how to install with pip or download/run tests.
* Added routines to network to get bus-connected components from names and bus numbers uses internal hash tables.
* Updated examples, documentation (macros, intersphinx).
* Distributed pfnet python wrapper through pypi.

Version 1.3.1
-------------
* Support for init values for constraint extra variables.
* Improved post processing of structures of constraint Hessians (now constr.c ensures lower triangular and fills H_combined ij, fixing bug in AC_FLOW_LIM).
* Mat parser detection of branches out of service.
* H_combined is now completely handled by "base" constraint in constr.c. Custom constraints no longer need to allocate this matrix.
* -Wall -Werror had no effect in Makefile.am and were moved to configure.ac. Now they work (requires autoconf-archive).
* Support for adding nonlinear constraints in Python and documentation.
* JSON network representation and parser with read/write capability.
* Base parser defaults to number of periods associated with data file and format-specific parsers are responsible for "propagating data on time".
* Added version string to info dictionary of pfnet python module.

Version 1.3.0
-------------
* Load Q variables.
* Support for variable generator curtailments.
* User-friendly way to add batteries and variable generators.
* Improved naming consistency in Python network class.
* Support for all battery and load variables in ACPF and network properties.
* Battery dynamics and boundary conditions.
* Updated documentation and examples.
* Full support for constraint auxiliary variables (lin eq, nonlin eq, line ineq).
* Elimination of obscure variables "voltage magnitude deviation".
* Elimination of obscure variables "voltage magnitude violation", "tap ratio deviation", "susceptance deviation".
* Improved setup.py without argparse that relies on existing build_ext commands for custom builds.
* Upper and lower bounds for constraint extra variables.
* Three types of voltage magnitude limits (normal, regulation, emergency).
* Load power factor, target power factor, and constraint for constant power factor.
* Removed Problem "set_network" method and required that Problem constructor takes network as argument for consistency with functions and constraints.
* Added default arguments to routines "add_batteries" and "add_var_generators".
* Fixed indentation bug in problem.show().
* (Conservative) linearized AC thermal limits via external and optional "line_flow" library.
* Changed key "raw parser" to "raw_parser" in "info" dictionary of pfnet python wrapper, and added "line_flow".
* Added constraint/function/network error checks in problem analyze and eval routines.
* Added branch phase and ratio python setters.
* Fixed sign error with second derivative of current mag with respect to phase shift in AC_FLOW_LIM constraints.

Version 1.2.9
-------------
* Changed function Hessian nnz counter (Hcounter to Hphi_nnz).
* Problem add_function takes Function object/struct and not function name.
* Problem add_constraint takes Constraint object/struct and not function name.
* Function class constructor takes function name for convenience.
* Constraint class constructor takes constraint name for convenience.
* CustomFunction class for creating functions written in python that work with the C library.
* CustomConstraint class for creating constraints written in python that work with the C library.
* Changed constraint (file) names PF and BOUND to ACPF and NBOUND, respectively.
* Separated netowrk and parser.
* (Generic) Parser class constructor takes file extension or sample filename for convenience.
* Format-specific parsers for mat,raw,art files also available (ParserMAT,ParserRAW,ParserART).
* CustomParser class for creating parsers written in python that work with the C library.
* Raw parser availability is detected when library is placed in pfnet/lib folder (no environment variable needed anymore).

Version 1.2.8
-------------
* Python-based parsers.
* Info dictionary in Python wrapper that indictes availability of graphviz, python-based parser, raw parser.
* Autotools build system (automatic detection of raw_parser, graphviz, python-based parser).
* Travis continuous integration.
* Readthedocs-based documentation.
* Fixed mpc2mat and parser_MAT gen cost ordering.

Version 1.2.7
-------------
* Added constraint that enforces AC branch flow limits using current magnitudes (ignores branches with 0 ratingA).
* Added support for extra variables in constraints (Jbar, Gbar matrices) and in problem.
* Changed constraint nnz counters (Acounter,Jcounter,Gcounter) to (A_nnz,J_nnz,G_nnz).
* Changed constraint row counters (Aconstr_index,Jcounstr_index,Gconstr_index) to (A_row,J_row,G_row).
* Made constraint that enforces DC branch flow limits ignore branches with 0 ratingA.

Version 1.2.6
-------------
* Branch bus name changes (from/to to k/m).
* Branch AC flow getters.
* Improved error handling in Problem Python class (has_error, clear_error, error checks in combine_H).
* Separated python wrapper pyx into multiple files.
* Bug fix: voltage magnitude limits in MAT parser.
* Sphinx C docs.

Version 1.2.5
-------------
* Strings instead of constants in python wrapper for object types, flag types, properties, object quantities, function and constraint types.

Version 1.2.4
-------------
* Multi-period support.
* Function value independent of variable flags.
* Name changes in Python wrapper (trying to eliminate abbreviations, e.g. gens, bats, etc).
* Improved memory management and bookkeeping when adding variable generators.
* Routine for getting number of variables of Bus.

Version 1.2.3
-------------
* Linear power flow constraints (LINPF).
* Improved Makefile (Linux and Mac).
* Documentation build rules with PFNET_DOCS.
* Removed DEBUG conditional compilation and added output levels to parsers.

Version 1.2.2
-------------
* Battery objects.
* Net consumption function.
* Bus price attribute.
* Load support to FIX constraint.

Version 1.2.1
-------------
* Python 3 and Jupyter compatibility.
* Variable load active powers.

Version 1.2.0
-------------
* Variable generators (e.g., wind and solar)
* l <= Gx <= u constraints.
* LBOUND constraint type.
* Contingencies.

Version 1.1
-----------
* Artere parser.

Version 1.0
-----------
* Initial version.
