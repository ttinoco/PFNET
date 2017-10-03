#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015, Tomas Tinoco De Rubira.       #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import sys
import numpy as np
from subprocess import call
from Cython.Build import cythonize
from setuptools import setup, Extension

# C library build
if 'darwin' in sys.platform.lower() or 'linux' in sys.platform.lower():
    return_code = call(["./build_lib.sh"])
else:
    return_code = call(["build_lib.bat"])
if return_code != 0:
    raise ValueError('Unable to build C library')

# Extra link args
if 'darwin' in sys.platform.lower():
    extra_link_args=['-Wl,-rpath,@loader_path/']
elif 'linux' in sys.platform.lower():
    extra_link_args=['-Wl,-rpath=$ORIGIN']
else:
    extra_link_args=['']

setup(name='PFNET',
      zip_safe=False,
      version='1.3.1rc9',
      description='Power Flow Network Library',
      url='https://github.com/ttinoco/PFNET/python',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      include_package_data=True,
      license='BSD 2-Clause License',
      packages=['pfnet',
                'pfnet.parsers',
                'pfnet.functions',
                'pfnet.constraints',
                'pfnet.tests'],
      install_requires=['cython>=0.20.1',
                        'numpy>=1.11.2',
                        'scipy>=0.18.1',
                        'nose'],
      package_data={'pfnet':['libpfnet*']},
      classifiers=['Development Status :: 5 - Production/Stable',
                   'License :: OSI Approved :: BSD License',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.5'],
      ext_modules=cythonize([Extension(name="pfnet.cpfnet",
                                       sources=["./pfnet/cpfnet.pyx"],
                                       libraries=['pfnet'],
                                       include_dirs=[np.get_include(),'./lib/pfnet/build/include'],
                                       library_dirs=['./lib/pfnet/build/lib'],
                                       extra_link_args=extra_link_args)]))
                                       
