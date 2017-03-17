import pfnet
net = pfnet.Network()
net.load('../data/ieee14.mat')

f = pfnet.functions.DummyGenCost(1.,net)

net.set_flags('generator',
              'variable',
              'any',
              'active power')

x = net.get_var_values()

f.analyze()
