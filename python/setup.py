1#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import numpy as np
from Cython.Build import cythonize
from distutils.core import setup, Extension

setup(name='PFNET',
      version='1.2.9',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      url='https://github.com/ttinoco/PFNET/python',
      packages=['pfnet',
                'pfnet.parsers',
                'pfnet.functions',
                'pfnet.constraints'],
      ext_modules=cythonize([Extension(name="pfnet.cpfnet", 
                                       sources=["./pfnet/cpfnet.pyx"],
                                       include_dirs=[np.get_include()])]))

