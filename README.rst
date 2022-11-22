AlloViz
=======

.. image:: https://readthedocs.org/projects/alloviz/badge/?version=latest
    :target: https://alloviz.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Code style

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
`namd <https://www.ks.uiuc.edu/Research/namd/>`__),
`MDEntropy <https://github.com/msmbuilder/mdentropy>`__ and 
`CARDS <https://github.com/sukritsingh/cardsReader>`__.

For the same topology and molecular dynamics (MD) trajectory, the
network can be constructed based on residue contacts,
correlation of atom movement or dihedrals, or interaction energies,
depending on the package selected. Moreover, for example for movement
correlation, the movement tracked can be that of the whole residue, its
center of mass, its alpha-C or its beta-C; and it can be calculated as
the Pearson’s correlation coefficient, Mutual Information (MI) or Linear
MI (LMI). See all the `options <https://alloviz.readthedocs.io/en/latest/table.html>`__.

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
   conda create -n AlloViz --file AlloViz/conda_explicit.txt
   conda activate AlloViz

Then go to the package folder (``cd AlloViz``) and install the package,
preferably with ``pip install .``.

   If environment creation with `conda_explicit.txt` fails, the non-explicit requirements/dependencies file `conda_minimal.txt` can be used, providing the conda channel `conda-forge` (`-c conda-forge`).

Tutorial
--------

Check the `tutorial <https://alloviz.readthedocs.io/en/latest/tutorial.html>`__.

Cite
-------

License
---------


