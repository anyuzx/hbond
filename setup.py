from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
	cmdclass = {'build_ext':build_ext},
	ext_modules = [Extension("hbond",
							  sources = ["hbond.pyx","c_hbond.cpp"],
							  language="c++",
							  include_dirs = [numpy.get_include()])],
)