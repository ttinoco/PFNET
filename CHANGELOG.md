Targets
-------
* Mechanism for changing constraint parameters.
* MAT or ART data writer.
* CIM parser written in python.
* Rich object comparisons for all network components.
* Instructions for adding new functions and constraints in C.
* Area/zone projections.
* Variable generator curtailment cost.

Unreleased
----------

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
* Extended support for storing and retrieving constraint sensitivity information (done).
* Update Python wrapper documentation to show how to install with pip or download/run tests (done).
* Added routines to network to get bus-connected components from names and bus numbers uses internal hash tables (done).
* Update examples, documentation (macros, intersphinx) and create release (todo).
* Distribute pfnet python wrapper through pypi (todo).

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
