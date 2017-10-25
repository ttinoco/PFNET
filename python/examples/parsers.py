#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

# Parsers - Overview

import os
import sys
sys.path.append('.')
import pfnet

parser = pfnet.Parser(sys.argv[1])

network = parser.parse(sys.argv[1])

pfnet.ParserJSON().write(network, "new_network.json")

os.remove("new_network.json")
