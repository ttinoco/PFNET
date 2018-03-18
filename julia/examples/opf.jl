using pfnet

using pfnet

case = ARGS[1]
parser = Parser(case)

net = parse_case(parser, case)
show_components(net)

clear_flags(net)



@printf("%d %d\n",
        num_vars(net),
        (2*num_buses(net)-
         num_slack_buses(net) +
         2*num_generators_not_on_outage(net))*num_periods(net))

@printf("%d %d\n",
        num_bounded(net),
        (2*num_generators_not_on_outage(net)+
         num_buses(net))*num_periods(net))


