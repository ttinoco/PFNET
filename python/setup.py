#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import os
import numpy as np
from subprocess import call
from Cython.Build import cythonize
from setuptools import setup, Extension, find_packages

return_code = call(["./build_lib.sh"])

DIR = os.getcwd()+'/lib/pfnet/build'

setup(name='PFNET',
      version='1.3.1',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      url='https://github.com/ttinoco/PFNET/python',
      packages=find_packages(),
      ext_modules=cythonize([Extension(name="pfnet.cpfnet",
                                       sources=["./pfnet/cpfnet.pyx"],
                                       libraries=['pfnet'],
                                       include_dirs=[np.get_include(),DIR+'/include'],
                                       library_dirs=[DIR+'/lib'],
                                       extra_link_args=['-Wl,-rpath,'+DIR+'/lib'])]))
                                       
