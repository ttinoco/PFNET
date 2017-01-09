# Automake
# AM_CPPFLAGS = -I$(
AUTOMAKE_OPTIONS = subdir-objects

# libraries
lib_LTLIBRARIES = libpfnet.la

# sources (note all headers are in $PFNET/include/pfnet)

graph_src = graph/graph.c graph.h

math_src = math/matrix.c matrix.h \
	 	math/vector.c vector.h

net_src = net/bat.c bat.h \
		net/branch.c branch.h  \
		net/bus.c bus.h  \
		net/contingency.c contingency.h  \
		net/gen.c gen.h  \
		net/load.c load.h \
		net/net.c net.h  \
		net/shunt.c shunt.h  \
		net/vargen.c vargen.h

parser_src = parser/parser_ART.c parser_ART.h \
		parser/parser_CSV.c parser_CSV.h \
		parser/parser_MAT.c parser_MAT.h \
		parser/parser_RAW.c parser_RAW.h

problem_src = problem/constr.c constr.h \
		problem/func.c func.h \
		problem/heur.c heur.h \
		problem/heur_PVPQ.c heur_PVPQ.h \
		problem/problem.c problem.h

problem_constr_src = problem/constr/constr_BOUND.c constr_BOUND.h \
		problem/constr/constr_DC_FLOW_LIM.c constr_DC_FLOW_LIM.h \
		problem/constr/constr_DCPF.c constr_DCPF.h \
		problem/constr/constr_FIX.c constr_FIX.h \
		problem/constr/constr_GEN_RAMP.c constr_GEN_RAMP.h \
		problem/constr/constr_LBOUND.c constr_LBOUND.h \
		problem/constr/constr_LINPF.c constr_LINPF.h \
		problem/constr/constr_PAR_GEN_P.c constr_PAR_GEN_P.h \
		problem/constr/constr_PAR_GEN_Q.c constr_PAR_GEN_Q.h \
		problem/constr/constr_PF.c constr_PF.h \
		problem/constr/constr_REG_GEN.c constr_REG_GEN.h \
		problem/constr/constr_REG_SHUNT.c constr_REG_SHUNT.h \
		problem/constr/constr_REG_TRAN.c constr_REG_TRAN.h

problem_func_src = problem/func/func_GEN_COST.c func_GEN_COST.h \
		problem/func/func_LOAD_UTIL.c func_LOAD_UTIL.h \
		problem/func/func_NETCON_COST.c func_NETCON_COST.h \
		problem/func/func_REG_PHASE.c func_REG_PHASE.h \
		problem/func/func_REG_PQ.c func_REG_PQ.h \
		problem/func/func_REG_RATIO.c func_REG_RATIO.h \
		problem/func/func_REG_SUSC.c func_REG_SUSC.h \
		problem/func/func_REG_VANG.c func_REG_VANG.h \
		problem/func/func_REG_VMAG.c func_REG_VMAG.h \
		problem/func/func_SLIM_VMAG.c func_SLIM_VMAG.h \
		problem/func/func_SP_CONTROLS.c func_SP_CONTROLS.h

utils_src = utils/utils.c utils.h

libpfnet_la_SOURCES = $(graph_src) $(math_src) $(net_src) $(parser_src) \
		      $(problem_src) $(problem_constr_src) $(problem_func_src) $(utils_src) \
		      array.h constants.h flag_types.h list.h obj_types.h \
		      types.h uthash.h

libpfnet_la_CFLAGS = -I../include
# libpfnet_la_LDFLAGS = -shared  # want both static and shared
libpfnet_la_LIBADD = -L../lib -lm -lraw_parser

