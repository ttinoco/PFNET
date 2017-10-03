#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import glob

DATA_DIR = './data/'
MAT_CASES = glob.glob(DATA_DIR+'*.mat')
RAW_CASES = glob.glob(DATA_DIR+'*.raw')
CASES = glob.glob(DATA_DIR+'*')
