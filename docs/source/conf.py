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
	'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    "nbsphinx",
    "nbsphinx_link",
]

nbsphinx_execute = 'never'
nbsphinx_kernel_name = 'python3'

intersphinx_mapping = {
'https://docs.mdanalysis.org/stable/': None,
}

# Make sure the target is unique
autosectionlabel_prefix_document = True

# autodoc
autodoc_member_order = "bysource"

# Napoleon settings
napoleon_google_docstring = False
napoleon_include_special_with_doc = False
napoleon_use_ivar = False # default
napoleon_use_param = True # default
napoleon_use_keyword = True # default
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
#html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']



sys.path.insert(0, os.path.abspath('../..'))
from AlloViz.AlloViz import info

df = info.df

header = "Available information sources for network generation\n----------------------------------------------------\n\n.. raw:: html\n\n"

tabulated_table = ""
for line in df.to_html(header=False).replace(' valign="top"', '').split("\n"):
	tabulated_table += f"\t{line}\n"


import fileinput, sys

replace_line = False
for line in fileinput.input("../../README.rst", inplace=True):
	if line.startswith('Cite'):
		replace_line = False

	if replace_line:
		continue

	if line.startswith('Available information sources for network generation'):
		line = header + tabulated_table + "\n|\n\n"
		replace_line = True

	sys.stdout.write(line)


sys.path.pop(0)