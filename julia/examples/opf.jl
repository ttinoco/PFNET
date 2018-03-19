using pfnet

# Case
case = ARGS[1]

# Parser
parser = Parser(case)

# Network
net = parse_case(parser, case)
show_components(net)

# Clear flags
clear_flags(net)

# Voltage magnitudes
set_flags(net,
          "bus",
          ["variable", "bounded"], 
          "any",
          "voltage magnitude")

# Voltage angles
set_flags(net,
          "bus",
          "variable",
          "not slack",
          "voltage angle")

# Generator powers
set_flags(net,
          "generator",
          ["variable","bounded"],
          "not on outage",
          ["active power","reactive power"])

assert(num_vars(net) == (2*num_buses(net)-
                         num_slack_buses(net) +
                         2*num_generators_not_on_outage(net))*num_periods(net))

assert(num_bounded(net) == (2*num_generators_not_on_outage(net)+
                            num_buses(net))*num_periods(net))

# Objective function
gen_cost = pfnet.Function("generation cost", 1., net)

# Constaints
acpf = pfnet.Constraint("AC power balance", net)
bounds = pfnet.Constraint("variable bounds", net)
th_limits = pfnet.Constraint("linearized AC branch flow limits", net)

# Problem
problem = pfnet.Problem(net)
add_function(problem, gen_cost)
add_constraint(problem, acpf)
add_constraint(problem, bounds)
add_constraint(problem, th_limits)
analyze(problem)
show_problem(problem)
