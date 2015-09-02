pfnet.load_library

net = pfnet.Network();

net.load('../data/ieee14.mat')

% Show
net.show_components()
net.show_properties()

% Numbers
assert(net.num_buses == 14);
assert(net.num_gens == 5);
assert(net.num_branches == 20);
assert(net.num_shunts == 1);
assert(net.num_loads == 11);

% Buses
bus = net.get_bus(5);
assert(bus.index == 5);

% Branches
br = net.get_branch(3);
assert(br.index == 3);

% Gens
gen = net.get_gen(4);
assert(gen.index == 4);

% Shunts
sh = net.get_shunt(0);
assert(sh.index == 0);

% Load
load = net.get_load(9);
assert(load.index == 9);



