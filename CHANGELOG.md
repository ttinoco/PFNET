Targets
-------
* Mechanism for changing constraint parameters.
* MAT or ART data writer.
* CIM parser written in python.
* Rich object comparisons for all network components.
* Instructions for adding new functions and constraints.
* Area/zone projections.
* Battery initial and final energy levels.

Unreleased
----------
* Python-based parsers.
* Info dictionary in Python wrapper that indictes availability of graphviz, python-based parser, raw parser.
* Autotools build system (automatic detection of raw_parser, graphviz, python-based parser).
* Travis continuous integration.
* Readthedocs-based documentation.

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
