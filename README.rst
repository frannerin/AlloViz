.. raw:: html

   <!-- [![PyPI - Python Version](https://img.shields.io/pypi/pyversions/correlationplus)](https://pypi.org/project/correlationplus/)
   [![PyPI](https://img.shields.io/pypi/v/correlationplus)](https://pypi.org/project/correlationplus/)
   [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/correlationplus/README.html)
   [![Open Source License: GPL v3](https://img.shields.io/badge/License-LGPLv3-blue.svg)](https://opensource.org/licenses/LGPL-3.0)
   [![Doc](https://readthedocs.org/projects/correlationplus/badge/?version=latest)](http://correlationplus.readthedocs.org/en/latest/#)
   [![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/structuraldynamicslab/correlationplus/latest)](https://hub.docker.com/repository/docker/structuraldynamicslab/correlationplus)
   ![Conda](https://img.shields.io/conda/pn/bioconda/correlationplus)
   [![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/tekpinar/correlationplus/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/tekpinar/correlationplus) -->

..
	.. role::  raw-html(raw)
	    :format: html
	:raw-html:`&#128301;`

AlloViz
=======

|Docs build|

.. |Maintenance yes| image:: https://readthedocs.org/projects/alloviz/badge/?version=latest
   :target: https://alloviz.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

A Python package to interactively compute, analyze and visualize protein
allosteric communication (residue interaction) networks and
delta-networks.

AlloViz binds together some newly written modules with 8 Python packages
that provide different ways of calculating residue interactions:
`GetContacts <https://github.com/getcontacts/getcontacts>`__,
`correlationplus <https://github.com/tekpinar/correlationplus>`__,
`dynetan <https://github.com/melomcr/dynetan>`__,
`PyInteraph2 <https://github.com/ELELAB/pyinteraph2>`__,
`pytraj <https://github.com/Amber-MD/pytraj>`__,
`MD-TASK <https://github.com/RUBi-ZA/MD-TASK>`__,
`gRINN <https://bitbucket.org/onursercinoglu/grinn>`__ (needs
`namd <https://www.ks.uiuc.edu/Research/namd/>`__) and
`MDEntropy <https://github.com/msmbuilder/mdentropy>`__.

.. raw:: html

   <!-- [g_correlation](https://www.mpinat.mpg.de/grubmueller/g_correlation), [GSAtools](https://github.com/AllosterIt/GSAtools) -->

For the same topology and molecular dynamics (MD) trajectory, the
network can be constructed based on residue contact frequency,
correlation of atom movement or dihedrals, or interaction energies,
depending on the package selected. Moreover, for example for movement
correlation, the movement tracked can be that of the whole residue, its
center of mass, its alpha-C or its beta-C; and it can be calculated as
the Pearson’s correlation coefficient, Mutual Information (MI) or Linear
MI (LMI). See
`below <#available-information-sources-for-network-generation>`__.

The resulting network can be analyzed with edge centrality metrics
algorithms provided by the Python package
`NetworkX <https://github.com/networkx/networkx>`__, and they can be
visualized in an interactive Python Notebook (i.e.,
`Jupyter <https://jupyter.org/>`__) using
`nglview <https://github.com/nglviewer/nglview>`__.

Install
-------

It is recommended to use a virtual environment
(`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__). This
repository includes submodules that need to be appropriately cloned
alongside the main repository using the ``--recursive`` flag. At
present, virtual environment dependencies can only be correctly
installed with conda.

.. code:: bash

   git clone --recursive --shallow-submodules -j 9 https://github.com/frannerin/AlloViz
   conda create AlloViz --file AlloViz/conda_explicit.txt
   conda activate AlloViz

Then go to the package folder (``cd AlloViz``) and install the package,
preferably with ``pip install .``.

   > If environment creation with `conda_explicit.txt` fails, the non-explicit requirements/dependencies file `conda_minimal.txt` can be used, providing the conda channel `conda-forge` (`-c conda-forge`).

Quickstart
----------

The package is designed for use in interactive notebooks. The main class
to be used is ``Protein``; ``Delta`` allows to bind two Proteins for
delta-network analysis.

for GPCRmd files
~~~~~~~~~~~~~~~~

A Protein is simply defined with the GPCRmd dynamics ID number, with
which files are automatically retrieved from the database. For example:

.. code:: python

   import AlloViz
   activeMuOR = AlloViz.Protein(GPCR=169)

And then network computations, analyses and visualizations are performed
with associated class methods. For example:

.. code:: python

   activeMuOR.calculate(pkg = "pytraj_CA")
   activeMuOR.analyze(metrics = "btw")
   activeMuOR.view("pytraj_CA", "btw")

Available information sources for network generation
----------------------------------------------------

.. raw:: html

	<table border="1" class="dataframe">
	  <thead>
	    <tr>
	      <th>Residue information extracted from trajectories</th>
	      <th>Package</th>
	      <th>Correlation measurement</th>
	      <th>Atom/angle tracked</th>
	      <th></th>
	    </tr>
	  </thead>
	  <tbody>
	    <tr>
	      <th rowspan="9">Movement correlation</th>
	      <th>MD-TASK</th>
	      <th>Pearson's</th>
	      <th>Carbon α</th>
	      <td>MDTASK</td>
	    </tr>
	    <tr>
	      <th rowspan="2">pytraj</th>
	      <th rowspan="2">Pearson's</th>
	      <th>Carbon α</th>
	      <td>pytraj_CA</td>
	    </tr>
	    <tr>
	      <th>Carbon β</th>
	      <td>pytraj_CB</td>
	    </tr>
	    <tr>
	      <th rowspan="2">dynetan</th>
	      <th rowspan="2">Mutual Information (MI)</th>
	      <th>Whole residue</th>
	      <td>dynetan</td>
	    </tr>
	    <tr>
	      <th>Residue COM</th>
	      <td>dynetan_COM</td>
	    </tr>
	    <tr>
	      <th rowspan="4">correlationplus</th>
	      <th rowspan="2">Pearson's</th>
	      <th>Carbon α</th>
	      <td>correlationplus_CA_Pear</td>
	    </tr>
	    <tr>
	      <th>Residue COM</th>
	      <td>correlationplus_COM_Pear</td>
	    </tr>
	    <tr>
	      <th rowspan="2">Linear MI (LMI)</th>
	      <th>Carbon α</th>
	      <td>correlationplus_CA_LMI</td>
	    </tr>
	    <tr>
	      <th>Residue COM</th>
	      <td>correlationplus_COM_LMI</td>
	    </tr>
	    <tr>
	      <th rowspan="10">Dihedral correlation</th>
	      <th rowspan="4">correlationplus</th>
	      <th rowspan="4">Pearson's</th>
	      <th>Phi</th>
	      <td>correlationplus_Phi</td>
	    </tr>
	    <tr>
	      <th>Psi</th>
	      <td>correlationplus_Psi</td>
	    </tr>
	    <tr>
	      <th>Omega</th>
	      <td>correlationplus_Omega</td>
	    </tr>
	    <tr>
	      <th>Backbone dihedrals (Phi, psi and omega)</th>
	      <td>correlationplus_Dihs</td>
	    </tr>
	    <tr>
	      <th rowspan="4">AlloViz</th>
	      <th rowspan="4">MI</th>
	      <th>Phi</th>
	      <td>AlloViz_Phi</td>
	    </tr>
	    <tr>
	      <th>Psi</th>
	      <td>AlloViz_Psi</td>
	    </tr>
	    <tr>
	      <th>Omega</th>
	      <td>AlloViz_Omega</td>
	    </tr>
	    <tr>
	      <th>Backbone dihedrals (Phi, psi and omega)</th>
	      <td>AlloViz_Dihs</td>
	    </tr>
	    <tr>
	      <th rowspan="2">MDEntropy</th>
	      <th rowspan="2">MI</th>
	      <th>Backbone dihedrals (Phi, psi and omega)</th>
	      <td>MDEntropy_Dihs</td>
	    </tr>
	    <tr>
	      <th>Alpha angle</th>
	      <td>MDEntropy_AlphaAngle</td>
	    </tr>
	    <tr>
	      <th rowspan="3">Contact frequency</th>
	      <th>MDEntropy</th>
	      <th>MI</th>
	      <th>Whole residue</th>
	      <td>MDEntropy_Contacts</td>
	    </tr>
	    <tr>
	      <th>GetContacts</th>
	      <th>-</th>
	      <th>Whole residue</th>
	      <td>GetContacts</td>
	    </tr>
	    <tr>
	      <th>PyInteraph2</th>
	      <th>-</th>
	      <th>Whole residue</th>
	      <td>PyInteraph2_Contacts</td>
	    </tr>
	    <tr>
	      <th>Interaction energy</th>
	      <th>PyInteraph2</th>
	      <th>-</th>
	      <th>Whole residue</th>
	      <td>PyInteraph2_Energy</td>
	    </tr>
	  </tbody>
	</table>

|

Cite
-------

License
---------
