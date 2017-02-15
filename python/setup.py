#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy
import os
import sys

# simple function to process command line arguments
def get_args(arg_list, arg):
    # get indices
    vals = [x.partition("=")[2] for x in arg_list if x.startswith(arg)]
    # clear out parameters
    arg_list = [x for x in arg_list if not x.startswith(arg)]
    # return information
    return (vals, arg_list)

# collect graphviz option
no_gvc, sys.argv = get_args(sys.argv, '--no_graphviz')

# collect raw-parser option
no_raw_parser, sys.argv = get_args(sys.argv, '--no_raw_parser')

# collect pfnet static library
pfnet_args, sys.argv = get_args(sys.argv, '--pfnet_lib')

# collect include directories
include_dirs, sys.argv = get_args(sys.argv, '--include_dir')

# collect library directories
library_dirs, sys.argv = get_args(sys.argv, '--library_dir')

# add numpy directories
include_dirs.append(numpy.get_include())

libraries = []
extra_compile_args = []
extra_objects = []

# graphviz
if no_gvc:
    extra_compile_args.append("-DNO_GRAPHVIZ")
else:
    libraries.append('gvc')

# raw-parser
if no_raw_parser:
    extra_compile_args.append("-DNO_RAW_PARSER")
else:
    libraries.append('raw_parser')

if pfnet_args:
    # static link
    pfnet_lib = pfnet_args[-1]
    extra_objects.append(pfnet_lib)
else:
    # dynamic link
    libraries.append("pfnet")

ext = Extension("pfnet.cpfnet",
                [os.path.join("pfnet", 'cpfnet.pyx')],
                libraries=libraries,
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                extra_objects=extra_objects,
                extra_compile_args=extra_compile_args)

setup(name='PFNET',
      version='1.2.5',
      license='BSD 2-clause license',
      description='Power Flow Network Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      url='https://github.com/ttinoco/PFNET',
      packages=['pfnet'],
      install_requires=['numpy', 'scipy'],
      ext_modules=cythonize([ext]))
