addpath(strcat(getenv('PFNET'),'/matlab'));

pfnet.load_library

net = pfnet.Network();
net.load('../../data/ieee14.mat');

g = pfnet.Graph(net);

g.set_layout()

g.write('ps2','graph.ps2');
