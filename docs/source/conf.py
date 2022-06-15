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
sys.path.insert(0, os.path.abspath('..'))#/..'), os.path.abspath('..'), os.path.abspath('../../src'), os.path.abspath('../../src/AlloViz'))


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
	#'readthedocs_ext.readthedocs'
    #'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    #'sphinx.ext.napoleon',
    #'sphinx.ext.napoleon',
    'numpydoc',
    'sphinx.ext.autosectionlabel',
    "myst_nb",
]

# Make sure the target is unique
autosectionlabel_prefix_document = True

# autodoc
autodoc_member_order = "bysource"
numpydoc_show_class_members = False

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = True
napoleon_custom_sections = None



# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']



sys.path.insert(0, os.path.abspath('../..'))
from AlloViz.AlloViz import info

df = info.df

# with open("../../README.rst", "a") as f:
# 	f.write(".. raw:: html\n\n")
# 	for line in df.to_html(header=False).replace(' valign="top"', '').split("\n"):
# 		f.write(f"\t{line}\n")

# with open("../../table.html", "w") as f:
# 	f.writelines(
# 		df.reset_index().to_markdown(tablefmt="grid", index=False)
# 		)


header = "Available information sources for network generation\n----------------------------------------------------\n\n.. raw:: html\n\n"

tabulated_table = ""
for line in df.to_html(header=False).replace(' valign="top"', '').split("\n"):
	tabulated_table += f"\t{line}\n"

import fileinput, sys
# for line in fileinput.input("../../README.rst", inplace=True):
#     if line.startswith('Available information sources for network generation'):
#     	line = header + tabulated_table
#     sys.stdout.write(line)

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


# with open("../../README.rst", "r") as in_file:
#     buf = in_file.readlines()

# with open("../../README.rst", "w") as out_file:
#     for line in buf:
#         if line.startswith('Available information sources for network generation'):
#         	line = header + tabulated_table
#         out_file.write(line)