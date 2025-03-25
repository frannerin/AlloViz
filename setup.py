import os, sys, numpy
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
                  include_dirs = [numpy.get_include(), "src/Packages/pyinteraph2/libinteract"])




  from distutils.sysconfig import get_python_lib
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




  # import platform
  # import distutils.ccompiler

  # # this code checks for OS. If OS is OSx then it checks for GCC as default compiler
  # #if GCC is the default compiler adds -fopenmp to linker and compiler args.
  # if 'darwin' in platform.system().lower():
  # 	if 'gcc' in  distutils.ccompiler.get_default_compiler():
  # 		use_openmp = True
  # 	else:
  # 		use_openmp = False
  # else:
  # 	use_openmp = True

  # extra_compile_args = ['-Wno-unreachable-code']
  # extra_link_args = []

  # if use_openmp:
  # 	extra_compile_args += ['-fopenmp']
  # 	extra_link_args = ['-fopenmp']





  enspara_extensions = [
    Extension(
        "enspara.info_theory.libinfo",
        ["src/Packages/enspara/enspara/info_theory/libinfo.pyx"],
        #extra_compile_args=extra_compile_args,
        #extra_link_args=extra_link_args,
        include_dirs=[numpy.get_include()],
    ), Extension(
        "enspara.geometry.libdist",
        ["src/Packages/enspara/enspara/geometry/libdist.pyx"],
        #extra_compile_args=extra_compile_args,
        #extra_link_args=extra_link_args,
        include_dirs=[numpy.get_include()],
    ), Extension(
        "enspara.msm.libmsm",
        ["src/Packages/enspara/enspara/msm/libmsm.pyx"],
        #extra_compile_args=extra_compile_args,
        #extra_link_args=extra_link_args,
        include_dirs=[numpy.get_include()],
    )]


  #from Cython.Build import cythonize
  ext_modules = [libinteract, libdistance] + enspara_extensions #cythonize()




else:

  
  # from src.AlloViz.AlloViz import info

  # df = info.df
  # df.index = pandas.MultiIndex.from_tuples(list(df.index), names=["Residue information extracted from trajectories",
  #                                                               "Package",
  #                                                               "Correlation measurement",
  #                                                               "Atom/angle tracked"])

  # with open("docs/source/table.html", "w") as f:
  #   f.writelines(
  #     df.to_html(buf=, header=False).replace(' valign="top"', '')
  #     )




	ext_modules = []





# Where the magic happens:
setup(
    ext_modules = ext_modules
    )
