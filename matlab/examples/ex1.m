addpath(strcat(getenv('PFNET'),'/matlab'));

pfnet.load_library

net = pfnet.Network();
net.load('../../data/ieee14.mat');

deg = 0;
for i=1:net.num_buses
    bus = net.get_bus(i-1);
    deg = deg + bus.degree/net.num_buses;
end
disp(deg);

