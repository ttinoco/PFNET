import pfnet as pf
from numpy import hstack
from numpy.linalg import norm
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve

def NRsolve(net):

    net.clear_flags()

    # bus voltage angles
    net.set_flags('bus',
                  'variable',
                  'not slack',
                  'voltage angle')
    
    # bus voltage magnitudes
    net.set_flags('bus',
                  'variable',
                  'not regulated by generator',
                  'voltage magnitude')
    
    # slack gens active powers
    net.set_flags('generator',
                  'variable',
                  'slack',
                  'active power')
    
    # regulator gens reactive powers
    net.set_flags('generator',
                  'variable',
                  'regulator',
                  'reactive power')

    p = pf.Problem()
    p.set_network(net)
    p.add_constraint('AC power balance')  
    p.add_constraint('generator active power participation')  
    p.add_constraint('generator reactive power participation')
    p.analyze()
    
    x = p.get_init_point()
    p.eval(x)

    residual = lambda x: hstack((p.A*x-p.b,p.f))

    while norm(residual(x)) > 1e-4:
        x = x + spsolve(bmat([[p.A],[p.J]],format='csr'),-residual(x))
        p.eval(x)

    net.set_var_values(x)
    net.update_properties()

