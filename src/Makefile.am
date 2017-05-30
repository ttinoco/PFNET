# Automake options
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

# local paths for includes
inc_path = $(top_srcdir)/include/pfnet

nobase_nodist_include_HEADERS = $(inc_path)/pfnet_config.h

# libraries
lib_LTLIBRARIES = libpfnet.la

# sources (note all headers are in $PFNET/include/pfnet)
graph_src = 	graph/graph.c 

graph_hdr = 	$(inc_path)/graph.h

math_src = 	math/matrix.c \
		math/vector.c

math_hdr = 	$(inc_path)/matrix.h \
		$(inc_path)/vector.h

net_src = 	net/bat.c \
		net/branch.c \
		net/bus.c \
		net/contingency.c \
		net/gen.c \
		net/load.c \
		net/net.c \
		net/shunt.c \
		net/vargen.c

net_hdr = 	$(inc_path)/bat.h \
		$(inc_path)/branch.h \
		$(inc_path)/bus.h \
		$(inc_path)/contingency.h \
		$(inc_path)/gen.h \
		$(inc_path)/load.h \
		$(inc_path)/net.h \
		$(inc_path)/shunt.h \
		$(inc_path)/vargen.h

parser_src = 	parser/parser.c \
		parser/parser_ART.c \
		parser/parser_CSV.c \
		parser/parser_MAT.c \
		parser/parser_RAW.c \
		parser/parser_JSON.c

parser_hdr = 	$(inc_path)/parser.h \
		$(inc_path)/parser_ART.h \
		$(inc_path)/parser_CSV.h \
		$(inc_path)/parser_MAT.h \
		$(inc_path)/parser_RAW.h \
		$(inc_path)/parser_JSON.h

problem_src = 	problem/constr.c \
		problem/func.c \
		problem/heur.c \
		problem/heur_PVPQ.c \
		problem/problem.c

problem_hdr = 	$(inc_path)/constr.h \
	  	$(inc_path)/func.h \
		$(inc_path)/heur.h \
		$(inc_path)/heur_PVPQ.h \
		$(inc_path)/problem.h

problem_constr_src = 	problem/constr/constr_NBOUND.c  \
			problem/constr/constr_DC_FLOW_LIM.c \
			problem/constr/constr_AC_FLOW_LIM.c \
			problem/constr/constr_AC_LIN_FLOW_LIM.c \
			problem/constr/constr_DCPF.c \
			problem/constr/constr_ACPF.c \
			problem/constr/constr_LINPF.c \
			problem/constr/constr_FIX.c \
			problem/constr/constr_GEN_RAMP.c \
			problem/constr/constr_BAT_DYN.c \
			problem/constr/constr_LBOUND.c \
			problem/constr/constr_PAR_GEN_P.c \
			problem/constr/constr_PAR_GEN_Q.c \
			problem/constr/constr_REG_GEN.c \
			problem/constr/constr_REG_SHUNT.c \
			problem/constr/constr_REG_TRAN.c \
			problem/constr/constr_LOAD_PF.c

problem_constr_hdr = 	$(inc_path)/constr_NBOUND.h \
			$(inc_path)/constr_DC_FLOW_LIM.h \
			$(inc_path)/constr_AC_FLOW_LIM.h \
			$(inc_path)/constr_AC_LIN_FLOW_LIM.h \
			$(inc_path)/constr_DCPF.h \
			$(inc_path)/constr_ACPF.h \
			$(inc_path)/constr_LINPF.h \
			$(inc_path)/constr_FIX.h \
			$(inc_path)/constr_GEN_RAMP.h \
			$(inc_path)/constr_BAT_DYN.h \
			$(inc_path)/constr_LBOUND.h \
			$(inc_path)/constr_PAR_GEN_P.h \
			$(inc_path)/constr_PAR_GEN_Q.h \
			$(inc_path)/constr_REG_GEN.h \
			$(inc_path)/constr_REG_SHUNT.h \
			$(inc_path)/constr_REG_TRAN.h \
			$(inc_path)/constr_LOAD_PF.h

problem_func_src = 	problem/func/func_GEN_COST.c \
			problem/func/func_LOAD_UTIL.c \
			problem/func/func_NETCON_COST.c \
			problem/func/func_REG_PHASE.c  \
			problem/func/func_REG_PQ.c \
			problem/func/func_REG_RATIO.c  \
			problem/func/func_REG_SUSC.c  \
			problem/func/func_REG_VANG.c  \
			problem/func/func_REG_VMAG.c  \
			problem/func/func_SLIM_VMAG.c  \
			problem/func/func_SP_CONTROLS.c

problem_func_hdr = 	$(inc_path)/func_GEN_COST.h \
			$(inc_path)/func_LOAD_UTIL.h \
			$(inc_path)/func_NETCON_COST.h \
			$(inc_path)/func_REG_PHASE.h \
			$(inc_path)/func_REG_PQ.h \
			$(inc_path)/func_REG_RATIO.h \
			$(inc_path)/func_REG_SUSC.h \
			$(inc_path)/func_REG_VANG.h \
			$(inc_path)/func_REG_VMAG.h \
			$(inc_path)/func_SLIM_VMAG.h \
			$(inc_path)/func_SP_CONTROLS.h

utils_src = 	utils/utils.c 

utils_hdr = 	$(inc_path)/utils.h

other_hdr = 	$(inc_path)/array.h $(inc_path)/constants.h $(inc_path)/flag_types.h \
		$(inc_path)/list.h $(inc_path)/obj_types.h $(inc_path)/types.h \
		$(inc_path)/uthash.h $(inc_path)/pfnet.h

libpfnet_la_SOURCES = 	$(graph_src) $(math_src) $(net_src) $(parser_src) \
		      	$(problem_src) $(problem_constr_src) $(problem_func_src) $(utils_src)

# Have to move back a directory $PFNET/include/pfnet/*.h
libpfnet_la_CFLAGS = -I$(inc_path)/..
libpfnet_la_LDFLAGS = -shared # want both static and shared
libpfnet_la_LIBADD = -lm

pkginclude_HEADERS = 	$(graph_hdr) $(math_hdr) $(net_hdr) $(parser_hdr) \
			$(problem_hdr) $(problem_constr_hdr) $(problem_func_hdr) \
			$(utils_hdr) $(other_hdr)
