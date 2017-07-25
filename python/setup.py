#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import numpy as np
from Cython.Build import cythonize
from setuptools import setup, Extension

package_data = {}
extra_link_args = []

# need to check if building distributable wheel and make sure to package libpfnet.*
if 'bdist_wheel' in sys.argv:
    package_data={'pfnet': ["libpfnet.*"]}
    # if on OSX add loader_path to rpath, so libpfnet.* is located
    if sys.platform.lower() == 'darwin':
        extra_link_args.append("-Wl,-rpath,@loader_path/")
    # if on Linux add origin to rpath, so libpfnet.so is located
    elif 'linux' in sys.platform.lower():
	extra_link_args.append("-Wl,-rpath=$ORIGIN")

setup(name='PFNET',
      version='1.3.0',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      url='https://github.com/ttinoco/PFNET/python',
      include_package_data=True,
      package_data=package_data,
      packages=['pfnet',
                'pfnet.parsers',
                'pfnet.functions',
                'pfnet.constraints'],
      ext_modules=cythonize([Extension(name="pfnet.cpfnet", 
                                       sources=["./pfnet/cpfnet.pyx"],
                                       include_dirs=[np.get_include()],
                                       extra_link_args=extra_link_args
                                       )]))
