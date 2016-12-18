Targets
-------
* Mechanism for changing constraint parameters.
* MAT or ART data writer.
* CIM parser in python.
* Rich object comparisons for all network components.
* Need AC branch flow limits.
* Instructions for adding new functions and constraints.
* Use Autotools.

Unreleased
----------
* Area/zone projections.
* Branch bus name changes.
* Branch AC flow getters.
* Battery initial and final energy levels.
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
