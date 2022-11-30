# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'AlloViz'
copyright = '2022, Francho Nerín-Fonz'
author = 'Francho Nerín-Fonz'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'numpydoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
	#'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    "nbsphinx",
    "nbsphinx_link",
    "matplotlib.sphinxext.plot_directive",
    #"autoapi.extension"
]

autoapi_dirs = ["../../src/AlloViz/"]
autoapi_ignore = ['*migrations*', '^.*', '*-checkpoint.py*']
autoapi_add_toctree_entry = False
autoapi_keep_files = True

nbsphinx_execute = 'never'
nbsphinx_kernel_name = 'python3'

intersphinx_mapping = {
'https://docs.mdanalysis.org/stable/': None,
'https://biopython.org/docs/latest/api/': None,
'http://nglviewer.org/nglview/latest/': None,
'https://networkx.org/documentation/stable/': None,
'https://pandas.pydata.org/docs/': None,
'https://matplotlib.org/stable/': None,
}

# Make sure the target is unique
autosectionlabel_prefix_document = True

# #numpydoc WITH AUTOMODULE
# autosummary_generate = False
# #autosummary_imported_members = True
# numpydoc_show_class_members = False
# #numpydoc_class_members_toctree = False

#numpydoc WITHOUT AUTOMODULE
#autosummary_generate = ["API/public_api", "API/complete_api"]
#autosummary_imported_members = True
#numpydoc_show_class_members = True
#numpydoc_class_members_toctree = False
numpydoc_attributes_as_param_list = True # this doesn't seem to be used
#autosummary_imported_members = True
templates_path = ["_templates"]
autosummary_mock_imports = [
	"nglview", 
	"vmd-python", "getcontacts",
	"sklearn", 
	"AlloViz.Packages.mdentropy.mdentropy.version",
	"mdtraj",
	"libinteract",
	"prody", 
	"numba",
	"pytraj",
	"cython", "h5py", "python-louvain", "community", "colorama", #dynetan deps
	"pyprind", "panedr", "natsort", "click", "PyQt5",
	"enspara.info_theory.libinfo", "enspara.geometry.libdist", "enspara.msm.libmsm", "tables",
]

# autodoc
#autodoc_member_order = "bysource"
# autosummary doesn't use this
autodoc_default_options = {
	'member_order': "bysource",
	'private-members': True,
	'undoc-members': True,
	'show-inheritance': True,
	#'inherited-members': 'pandas.DataFrame',
}

# matplotlib plot directive # copied from pandas docs conf.py
plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_pre_code = """import numpy as np
import pandas as pd"""

# Napoleon settings
# napoleon_google_docstring = False
# napoleon_include_special_with_doc = False
# napoleon_use_ivar = False # default
# napoleon_use_param = True # default
# napoleon_use_keyword = True # default
# look into napoleon_type_aliases to put `Protein` as Delta parameters

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

html_theme_options = {
	"show_prev_next": False,

    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/frannerin/AlloViz",
            "icon": "fab fa-github",
            "type": "fontawesome",
        }],

    "logo": {
        "text": "AlloViz",},
}

html_sidebars = {
  "tutorial": [],
  "table": [],
}



sys.path.insert(0, os.path.abspath('../..'))
from AlloViz.AlloViz import info

df = info.df

# header = "Network construction methods\n============================\n\n.. raw:: html\n\n"
firstlines = """
Network construction methods
============================

.. raw:: html

"""

with open("table.rst", "w") as f:
	f.write(firstlines)

	for line in df.to_html(header=False).replace(' valign="top"', '').split("\n"):
		f.write(f"\t{line}\n")

# tabulated_table = ""
# for line in df.to_html(header=False).replace(' valign="top"', '').split("\n"):
# 	tabulated_table += f"\t{line}\n"


# import fileinput, sys

# replace_line = False
# for line in fileinput.input("../../README.rst", inplace=True):
# 	if line.startswith('Cite'):
# 		replace_line = False

# 	if replace_line:
# 		continue

# 	if line.startswith('Available information sources for network generation'):
# 		line = header + tabulated_table + "\n|\n\n"
# 		replace_line = True

# 	sys.stdout.write(line)


sys.path.pop(0)




with open("../../README.rst", 'r') as f:
  lines = f.read()

# lines[-3:] = [
# 	".. _options: table.html\n",
# 	".. _tutorial: tutorial.html\n",
# 	".. _documentation: index.html\n",
# ]



with open('README.rst', 'w') as f:
  f.write(lines.replace(
  		"`options <https://alloviz.readthedocs.io/en/latest/table.html>`__",
  		":ref:`options <table:Network construction methods>`"
  	).replace(
  		"`tutorial <https://alloviz.readthedocs.io/en/latest/tutorial.html>`__",
  		":ref:`tutorial <tutorial:Tutorial>`"
  	))