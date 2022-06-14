#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import sys
from shutil import rmtree
from glob import glob

from setuptools import find_packages, setup, Command, Extension

# Package meta-data.
NAME = 'AlloViz'
DESCRIPTION = 'A Python package to interactively computate, analyze and visualize protein allosteric communication (residue interaction) networks and delta-networks.'
URL = 'https://github.com/frannerin/AlloViz'
EMAIL = 'francho.nerin01@estudiant.upf.edu'
AUTHOR = 'Francho Nerín-Fonz'
REQUIRES_PYTHON = '>=3.7.0'
VERSION = '1.0.0'

# What packages are required for this module to be executed?
REQUIRED = [
	"multiprocess", # multiprocess and not multiprocessing to pickle object methods
	"mdanalysis>=2.0.0",
	"numpy>=1.21.0",
	"pandas>=1.3.5", "pyarrow>=6.0.0",
	"beautifulsoup4>=4.0.0", "certifi", "requests>=2.0.0", # for downloading GPCRmd files
	"mdtraj>=1.9", "vmd-python>=3", "matplotlib" # pkgs dependencies
]

# What packages are optional?
EXTRAS = {
	'interactive and visualization': ['jupyterlab', "nglview>=3.0.0"]
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()


# Necessary to create the version.py file for mdentropy to be able to use the original repo directly as a submodule
from Packages.mdentropy.basesetup import write_version_py
VERSION = "0.4.0dev0"
ISRELEASED = False
write_version_py(VERSION, ISRELEASED, 'Packages/mdentropy/mdentropy/version.py')





import numpy
#libinteract/innerloops
libinteract = \
      Extension("libinteract.innerloops",
                ["Packages/pyinteraph2/libinteract/innerloops.pyx",
                 "Packages/pyinteraph2/libinteract/clibinteract.c"], \
                include_dirs = [numpy.get_include()])



# Apparently the only needed Extension of msmbuilder to make mdentropy work; this block is mostly copied from msmbuilder's setup.py
#from basesetup import CompilerDetection
#if '--disable-openmp' in sys.argv:
#    sys.argv.remove('--disable-openmp')
#    DISABLE_OPENMP = True
#else:
#    DISABLE_OPENMP = False
#compiler = CompilerDetection(DISABLE_OPENMP)

from distutils.sysconfig import get_python_lib
mdtrajdir = get_python_lib() + "/mdtraj/core/lib"
mdtraj_capi = {'include_dir': mdtrajdir, 'lib_dir': mdtrajdir}

libdistance = \
    Extension('msmbuilder.libdistance',
              language='c++',
              sources=['Packages/msmbuilder/msmbuilder/libdistance/libdistance.pyx'],
              # msvc needs to be told "libtheobald", gcc wants just "theobald"
              libraries=['theobald'],#'%stheobald' % ('lib' if compiler.msvc else '')],
              include_dirs=["Packages/msmbuilder/msmbuilder/libdistance/src",
                            mdtraj_capi['include_dir'], numpy.get_include()],
              library_dirs=[mdtraj_capi['lib_dir']],
              )


# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    package_dir={"AlloViz": "AlloViz", "AlloViz.Packages": "Packages",
                 "pyinteraph": "Packages/pyinteraph2/pyinteraph", "libinteract": "Packages/pyinteraph2/libinteract",
                 "msmbuilder": "Packages/msmbuilder/msmbuilder"}, # let's test this: , "AlloViz.Packages": "Packages"
    # , "libdistance": "Packages/msmbuilder/msmbuilder/libdistance", 
    packages=["AlloViz", "AlloViz.Packages", "pyinteraph", "libinteract", "msmbuilder"], # "libinteract", "libdistance"
    #find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"])
    #data_files={"AlloViz/Forks": glob("AlloViz/Forks/*", recursive=True)},
    #data_files=[ ("Forks", glob("AlloViz/Forks/*", recursive=True)) ],
    #data_files=[ ("Forks", list(os.walk("AlloViz/Forks"))) ],
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    ext_modules = [libinteract, libdistance],
    #install_requires=[f"pytraj @ file://{os.getcwd()}/Packages/pytraj"],
    #install_requires = ["pytraj @ git+https://github.com/Amber-MD/pytraj/archive/refs/tags/v.2.0.5.zip"],
    #install_requires = ["pytraj @ git+https://github.com/Amber-MD/pytraj --depth 1"], #@2.0.5
    #install_requires=REQUIRED,
    #extras_require=EXTRAS,
    include_package_data=True,
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
"""License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Framework :: IPython
Framework :: Jupyter
Framework :: Jupyter :: JupyterLab
Intended Audience :: Science/Research
Intended Audience :: Healthcare Industry
Topic :: Scientific/Engineering :: Bio-Informatics""".split("\n")
    ]#,
# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev
    # $ setup.py publish support.
#    cmdclass={
#        'upload': UploadCommand,
#    },
)
