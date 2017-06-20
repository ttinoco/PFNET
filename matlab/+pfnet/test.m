pfnet.load_library

net = pfnet.Network();

net.load('../data/ieee14.mat')

% Constants
%%%%%%%%%%%
assert(exist('CURRENT','var')==1)
assert(exist('BUS_VAR_VANG','var')==1)

% Show
%%%%%%
net.show_components()
net.show_properties()

% Numbers
%%%%%%%%%
assert(net.num_buses == 14);
assert(net.num_gens == 5);
assert(net.num_branches == 20);
assert(net.num_shunts == 1);
assert(net.num_loads == 11);

% Buses
%%%%%%%
bus = net.get_bus(13);
assert(bus.index == 13);
bus.is_slack();
bus.is_regulated_by_gen();
bus.has_flags(FLAG_VARS,BUS_VAR_VMAG);
bus.get_total_gen_P();
bus.get_total_gen_Q();
bus.get_total_load_P();
bus.get_total_load_Q();
bus.show();
assert(bus.index_v_mag==0);
assert(bus.index_v_ang==0);
bus.number;
bus.degree;
bus.v_mag;
bus.v_ang;
bus.v_set;
bus.v_max_reg;
bus.v_min_reg;
bus.v_max_norm;
bus.v_min_norm;
bus.v_max_emer;
bus.v_min_emer;
assert(size(bus.gens,1)>0);
assert(size(bus.gens,2)==1);
assert(bus.gens{1}.bus.index==13);

% Branches
%%%%%%%%%%
br = net.get_branch(3);
assert(br.index == 3);

% Gens
%%%%%%
gen = net.get_gen(4);
assert(gen.index == 4);
assert(gen.bus.number > 0);

% Shunts
%%%%%%%%
sh = net.get_shunt(0);
assert(sh.index == 0);

% Load
%%%%%%
load = net.get_load(9);
assert(load.index == 9);
assert(load.bus.number > 0)

% Variables
%%%%%%%%%%%

assert(net.num_vars == 0);

net.set_flags(OBJ_BUS,....
	      FLAG_VARS,...
	      BUS_PROP_ANY,...
	      BUS_VAR_VMAG);
assert(net.num_vars > 0);
assert(net.num_vars == net.num_buses);

x = net.get_var_values();
assert(all(size(x) == [1 net.num_vars]));

% Graph
%%%%%%%

g = pfnet.Graph(net);
g.set_layout()
%g.write('ps2','graph.ps2');

% Function
%%%%%%%%%%

f = pfnet.Function(FUNC_TYPE_REG_VMAG,2.5,net);
assert(f.type==FUNC_TYPE_REG_VMAG);
assert(f.type~=FUNC_TYPE_REG_VANG);
assert(f.weight==2.5);
assert(f.phi==0.);
f.analyze();
f.eval(x);
assert(f.phi>0.);
assert(all(size(f.gphi)==size(x)));
