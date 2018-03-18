
using pfnet

case = ARGS[1]

parser = Parser(case)

net = parse_case(parser, case)

show_components(net)
