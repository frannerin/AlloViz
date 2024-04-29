AlloViz
=======

.. image:: https://readthedocs.org/projects/alloviz/badge/?version=latest
    :target: https://alloviz.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://github.com/frannerin/AlloViz/actions/workflows/test_conda_newenv.yml/badge.svg?branch=main
   :target: https://github.com/frannerin/AlloViz/actions/workflows/test_conda_newenv.yml
   :alt: Conda installation

.. image:: https://github.com/frannerin/AlloViz/actions/workflows/test_pip_tcl_ubuntu_newenv.yml/badge.svg?branch=main
   :target: https://github.com/frannerin/AlloViz/actions/workflows/test_pip_tcl_ubuntu_newenv.yml
   :alt: pip installation
   

AlloViz is a Python package to interactively compute, analyze and visualize protein
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
`MDEntropy <https://github.com/msmbuilder/mdentropy>`__ and 
`CARDS <https://github.com/sukritsingh/cardsReader>`__.

..
    `gRINN <https://bitbucket.org/onursercinoglu/grinn>`__ (needs
    `namd <https://www.ks.uiuc.edu/Research/namd/>`__),

For the same topology and molecular dynamics (MD) trajectory, the
network can be constructed based on residue contacts,
correlation of atom movement or dihedrals, or interaction energies,
depending on the package selected. Moreover, for example for movement
correlation, the movement tracked can be that of the whole residue, its
center of mass, its alpha-C or beta-C; and it can be calculated as
the Pearson’s correlation coefficient, Mutual Information (MI) or Linear
MI (LMI). See all the `options <https://alloviz.readthedocs.io/en/latest/table.html>`__.

The resulting network can be analyzed with edge centrality metrics
algorithms provided by the Python package
`NetworkX <https://github.com/networkx/networkx>`__, and they can be
visualized in an interactive Python Notebook (i.e.,
`Jupyter <https://jupyter.org/>`__) using
`nglview <https://github.com/nglviewer/nglview>`__.

AlloViz can also be use through a `GUI <https://alloviz.readthedocs.io/en/latest/tutorials/gui.html>`__.

Installation
-------------------
1. Clone the repository


The repository must be cloned along with all the submodules using the ``--recursive`` flag.
Additional flags are recommended for speed:

.. code:: bash

   git clone --recursive --shallow-submodules -j 9 https://github.com/frannerin/AlloViz

---------------

2. Create the virtual environment


It is recommended to create a **virtual environment** with `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__
or similars using the ``conda-forge`` channel (a fast dependency solver is recommended for speed:  
`libamba solver for Miniconda <https://conda.github.io/conda-libmamba-solver/getting-started/>`__
or the `Mamba <https://mamba.readthedocs.io/en/latest/>`__ version of Conda):

.. code:: bash

   conda create -n alloviz -c conda-forge --solver libmamba --file AlloViz/conda_environment.txt
   conda activate alloviz

If you wish to create a virtual environment without conda, see below for the alternative procedure.


Finally, AlloViz is installed into the environment with: 

.. code:: bash

    pip install ./AlloViz


**⚠ Note for MacOS M1/M2 (ARM) users.** 
Porting of conda software to the ARM architecture is hit-and-miss, and
numerous dependencies are missing.
You may want to install x64 packages instead, as follow:

.. code:: bash

   CONDA_SUBDIR=osx-64 conda create -n alloviz -c conda-forge --solver libmamba --file AlloViz/conda_environment.txt

Then activate the environment and install AlloViz using ``pip install ./AlloViz``



2.1 Create the virtual environment - alternative procedure



Although not recommended, the virtual environment can also be created with **pip**:

.. code:: bash

   python -m venv alloviz/env
   source alloviz/env/bin/activate
   pip install -r AlloViz/pip_requirements.txt
   pip install ./AlloViz

..


   Python <3.10 is recommended (i.e., 3.9.16). ``pytraj`` and the construction of delta-networks won't be available in a pip environment,
   as `AmberTools <http://ambermd.org/AmberTools.php>`__ and `pymol-open-source <https://github.com/schrodinger/pymol-open-source/>`__ 
   are needed (respectively) for that, and they aren't distributed through PyPi. Other additional dependencies might also need to be installed by hand.


Quickstart
--------

Check the `tutorial notebooks <https://alloviz.readthedocs.io/en/latest/tutorials.html>`__ or the
`quickstart <https://alloviz.readthedocs.io/en/latest/tutorials/quickstart.html>`__.

Cite
-------

License
---------


