Targets
-------
* Mechanism for changing constraint parameters.
* MAT or ART data writer.
* Rich object comparisons for all network components.
* Need AC branch flow limits.
* Instructions for adding new functions and constraints.
* Make an Autotools build. (long-term)

Unreleased
----------
* Multi-period support.
* Option flag for getting scalar outputs when number of periods is one.
* Function value independent of variable flags.
* Bus number of variables.
* Area/zone projections.
* Branch bus name changes.
* Branch AC flow getters.
* Name changes in Python wrapper (trying to eliminate abbreviations, e.g. gens, bats, etc).
* Improved memory management and bookeeing when adding variable generators.

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
