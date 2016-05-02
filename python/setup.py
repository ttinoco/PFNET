#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
from Cython.Build import cythonize
from distutils.core import setup, Extension
import numpy as np
import os

libraries = ['pfnet']
# raw parser
if '--no_raw_parser' in sys.argv:
    sys.argv.remove('--no_raw_parser')
else:
    libraries.append('raw_parser')
    
# graphviz
if '--no_graphviz' in sys.argv:
    sys.argv.remove('--no_graphviz')
else:
    libraries.append('gvc')


module_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfnet')
pfnet_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # Go one folder up
library_path = os.path.join(pfnet_path, 'lib')
include_path = os.path.join(pfnet_path, 'include')

print("Module Path: " + module_path)
print("Include Path: " + library_path)
print("Library Path: " + include_path)


setup(name='PFNET',
      version='1.3',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      packages=['pfnet'],
      ext_modules=cythonize([Extension("pfnet.cpfnet", 
                                       [os.path.join(module_path, 'cpfnet.pyx')],
                                       libraries=libraries,
                                       library_dirs=[library_path],
                                       include_dirs=[include_path, np.get_include()])]))
      
