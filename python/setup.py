#***************************************************#
# This file is part of PFNET.                       #
#                                                   #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.  #
#                                                   #
# PFNET is released under the BSD 2-clause license. #
#***************************************************#

import numpy
import os
import sys
import argparse
import numpy as np
from Cython.Build import cythonize
#from distutils.core import setup, Extension
from setuptools import setup, Extension

parser = argparse.ArgumentParser()
parser.add_argument('--libdirs', dest='libdirs', action='store',nargs='*',default=[])
parser.add_argument('--incdirs', dest='incdirs', action='store',nargs='*',default=[])
parser.add_argument('--libpfnet', dest='libpfnet', action='store',nargs='*',default=[])
args,unknown = parser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

if args.libdirs:
    library_dirs=args.libdirs
    extra_link_args=["-Wl,-rpath,%s" %s for s in args.libdirs]
else:
    library_dirs=[]
    extra_link_args=["-Wl,-rpath,/usr/local/lib"]
if args.incdirs:
    include_dirs=[np.get_include()]+args.incdirs
else:
    include_dirs=[np.get_include(),"../include"]
if args.libpfnet:
    extra_objects=args.libpfnet
    libraries=[]
else:
    extra_objects=[]
    libraries=["pfnet"]

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
#      install_requires=['numpy', 'scipy'],
      ext_modules=cythonize([Extension("pfnet.cpfnet",
                                       [os.path.join("pfnet","cpfnet.pyx")],
                                       include_dirs=include_dirs,
                                       library_dirs=library_dirs,
                                       libraries=libraries,
                                       extra_compile_args=[],
                                       extra_objects=extra_objects,
                                       extra_link_args=extra_link_args)]))
