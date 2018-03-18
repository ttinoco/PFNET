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

@printf("%d %d\n",
        num_vars(net),
        (2*num_buses(net)-
         num_slack_buses(net) +
         2*num_generators_not_on_outage(net))*num_periods(net))

@printf("%d %d\n",
        num_bounded(net),
        (2*num_generators_not_on_outage(net)+
         num_buses(net))*num_periods(net))


