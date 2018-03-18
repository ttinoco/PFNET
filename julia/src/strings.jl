
# Objects
str2obj = Dict("all" => 0,
               "bus" => 1,
               "generator" => 2,
               "branch" => 3,
               "shunt" => 4,
               "load" => 5,
               "variable generator" => 6,
               "battery" => 7,
               "unknown" => 8)

# Flags
str2flag = Dict("variable" => 0x01,
                "fixed" => 0x02,
                "bounded" => 0x04,
                "sparse" => 0x08)

# Quantities
str2q_bus = Dict("all" => 0xFF,
                 "voltage magnitude" => 0x01,
                 "voltage angle" => 0x02)

str2q_branch = Dict("all" => 0xFF,
                    "tap ratio" => 0x01,
                    "phase shift" => 0x02)

str2q_gen = Dict("all" => 0xFF,
                 "active power" => 0x01,
                 "reactive power" => 0x02)

str2q_shunt = Dict("all" => 0xFF,
                   "susceptance" => 0x01)

str2q_load = Dict("all" => 0xFF,
                  "active power" => 0x01,
                  "reactive power" => 0x02)

str2q_vargen = Dict("all" => 0xFF,
                    "active power" => 0x01,
                    "reactive power" => 0x02)

str2q_bat = Dict("all" => 0xFF,
                 "charging power" => 0x01,
                 "energy level" => 0x02)

str2q = Dict("all" => Dict("all" => 0),
             "bus" => str2q_bus,
             "branch" => str2q_branch,
             "generator" => str2q_gen,
             "shunt" => str2q_shunt,
             "load" => str2q_load,
             "variable generator" => str2q_vargen,
             "battery" => str2q_bat)

# Properties
str2prop_bus = Dict("any" => 0x00,
                    "slack" => 0x01,
                    "regulated by generator" => 0x02,
                    "regulated by transformer" => 0x04,
                    "regulated by shunt" => 0x08,
                    "not regulated by generator" => 0x10,
                    "not slack" => 0x20)

str2prop_branch = Dict("any" => 0x00,
                       "tap changer" => 0x01,
                       "tap changer - v" => 0x02,
                       "tap changer - Q" => 0x04,
                       "phase shifter" => 0x08,
                       "not on outage" => 0x10)

str2prop_gen = Dict("any" => 0x00,
                    "slack" => 0x01,
                    "regulator" => 0x02,
                    "not regulator" => 0x04,
                    "not slack" => 0x08,
                    "not on outage" => 0x10,
                    "adjustable active power" => 0x20)

str2prop_shunt =  Dict("any" => 0x00,
                       "switching - v" => 0x01)

str2prop_load = Dict("any" => 0x00,
                     "adjustable active power" => 0x01)

str2prop_vargen = Dict("any" => 0x00)

str2prop_bat = Dict("any" => 0x00)

str2prop = Dict("all" => Dict("all" => 0),
                "bus" => str2prop_bus,
                "branch" => str2prop_branch,
                "generator" => str2prop_gen,
                "shunt" => str2prop_shunt,
                "load" => str2prop_load,
                "variable generator" => str2prop_vargen,
                "battery" => str2prop_bat)


                 
                 
               


