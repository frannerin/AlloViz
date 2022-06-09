import os, sys, numpy
from distutils.sysconfig import get_python_lib
from setuptools import setup, Extension


on_rtd = os.environ.get('READTHEDOCS') == 'True'


if not on_rtd:

  # Necessary to create the version.py file for mdentropy to be able to use the original repo directly as a submodule
  # from .Packages.mdentropy.basesetup import write_version_py
  sys.path.append(".")
  from src.Packages.mdentropy.basesetup import write_version_py
  sys.path.append("./src/Packages/mdentropy")
  from src.Packages.mdentropy.setup import VERSION, ISRELEASED
  write_version_py(VERSION, ISRELEASED, 'src/Packages/mdentropy/mdentropy/version.py')





  libinteract = \
        Extension("libinteract.innerloops",
                  ["src/Packages/pyinteraph2/libinteract/innerloops.pyx",
                   "src/Packages/pyinteraph2/libinteract/clibinteract.c"], \
                  include_dirs = [numpy.get_include()])






  mdtrajdir = get_python_lib() + "/mdtraj/core/lib"
  mdtraj_capi = {'include_dir': mdtrajdir, 'lib_dir': mdtrajdir}

  libdistance = \
      Extension('msmbuilder.libdistance',
                language='c++',
                sources=['src/Packages/msmbuilder/msmbuilder/libdistance/libdistance.pyx'],
                # msvc needs to be told "libtheobald", gcc wants just "theobald"
                libraries=['theobald'],#'%stheobald' % ('lib' if compiler.msvc else '')],
                include_dirs=["src/Packages/msmbuilder/msmbuilder/libdistance/src",
                              mdtraj_capi['include_dir'], numpy.get_include()],
                library_dirs=[mdtraj_capi['lib_dir']],
                )




  ext_modules = [libinteract, libdistance]




else:
	ext_modules = []





# Where the magic happens:
setup(
    ext_modules = ext_modules
    )
