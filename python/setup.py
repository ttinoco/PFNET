#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import os
import sys
import numpy as np
from Cython.Build import cythonize
from distutils.core import setup, Extension
    
setup(name='PFNET',
      version='1.2.5',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      url='https://github.com/ttinoco/PFNET/python',
      packages=['pfnet',
                'pfnet.parser'],
      ext_modules=cythonize([Extension("pfnet.cpfnet", 
                                       [os.path.join("pfnet",'cpfnet.pyx')],
                                       library_dirs=[],
                                       libraries=['pfnet'],
                                       extra_compile_args=[],
                                       include_dirs=["../include",np.get_include()])]))
