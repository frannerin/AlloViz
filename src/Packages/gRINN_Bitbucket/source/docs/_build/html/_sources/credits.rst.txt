Credits
=======

The core functionality of gRINN is provided by NAMD and GROMACS MD simulation software. As such, gRINN can be considered a "wrapper" around these, but with a purpose.

gRINN was written in `Python programming language <https://www.python.org/>`_ (version 2.7). In addition to the Python core library, several open source Python packages are used to provide the following functionality:

`PyQt5 <https://sourceforge.net/projects/pyqt/>`_ (a python wrapper around the Qt desktop application development environment, v5.6.0) is used for all GUI elements.

`matplotlib <https://matplotlib.org/>`_ (v2.0.2) and `seaborn <https://seaborn.pydata.org>`_ (v0.8.1) are used to display two-dimensional line and scatter plots as well as heatmaps included in the "View Results" interface.

`ProDy <https://prody.csb.pitt.edu>`_ (v1.9.3) is used for PDB and DCD trajectory manipulations, atom selections and all other general geometric tasks related to protein structures. 

`Mdtraj <https://mdtraj.org>`_ (v1.9.0) is used to convert GROMACS trajectory file formats to DCD for further processing using ProDy.

`PyMol by Schrödinger <https://www.pymol.org>`_ (open source version v1.9.0.0) is used as the embedded molecule viewer in the "View Results" interface of gRINN.

pexpect (v4.3.1) is used to interact with gmx executable from within Python environment.

`Numpy <https://www.numpy.org>`_ (v1.13.3) is used for all operations related to matrices, which occur throughout the computation workflow of gRINN.

`pandas <https://pandas.pydata.org>`_ (v0.20.3) is used to store, process and save data throughout the computational workflow of gRINN.

`networkx <https://networkx.github.io>`_ (v2.0) is used to construct Protein Energy Networks and calculation of local network metrics and shortest paths.

gRINN's logo, Shamrock lucky Icon was designed by `www.iconka.com <http://www.iconka.com>`_