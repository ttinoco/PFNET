addpath(strcat(getenv('PFNET'),'/matlab'));

pfnet.load_library

net = pfnet.Network();
net.num_buses

net.load('../../data/ieee14.mat');
net.show_components();

bus = net.get_bus(10);

bus.index == 10

other_bus = net.get_bus_by_number(bus.number);

other_bus.index == bus.index

buses = net.buses;
reg_buses = [];
for i=1:net.num_buses
    bus = buses(i);
    if bus.is_regulated_by_gen()
       reg_buses = [bus reg_buses];
    end
end

s = size(reg_buses); 
[s(2) net.get_num_buses_reg_by_gen()]

branch = net.get_branch(5);

branch.index == 5;

branches = net.branches;
lines = [];
for i=1:net.num_branches
    branch = branches(i);
    if branch.is_line()
       lines = [branch lines];
    end
end

s = size(lines); 
[s(2) net.get_num_lines()]
    
gen = net.get_gen(2);

gen.index == 2

gens = net.generators;
slack_gens = [];
for i=1:net.num_gens
    gen = gens(i);
    if gen.is_slack()
       slack_gens = [gen slack_gens];
    end
end

s = size(slack_gens); 
[s(2) net.get_num_slack_gens()]

P = 0;
for i=1:net.num_gens
    P = P + gens(i).P*net.base_power;
end
P

shunt = net.get_shunt(0);

shunt.index == 0;

gbuses = net.get_gen_buses();

s = size(gbuses);


