load_pfnet

net = Network();

net.load('../data/ieee14.mat')

net.show_components()

net.num_buses

bus = net.get_bus(5);

bus.index
