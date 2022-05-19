Usage
=====

.. _installation:

Installation
------------

It is recommended to use a virtual environment ([Miniconda](https://docs.conda.io/en/latest/miniconda.html)). This repository includes submodules that need to be appropriately cloned alongside the main repository using the `--recursive` flag. At present, some virtual environment dependencies are exclusively available in Anaconda.

.. code-block:: console
   
   git clone --recursive --shallow-submodules -j 9 https://github.com/frannerin/AlloViz
   conda create AlloViz --file AlloViz/conda_explicit.txt
   conda activate AlloViz



Quickstart
----------------

The package is designed for use in interactive notebooks. The main class to be used is ``State``; ``Pair`` allows to bind two States for delta-network analysis.

## for GPCRmd files

A State is simply defined with the GPCRmd dynamics ID number, with which files are automatically retrieved from the database. For example:

>>> import AlloViz
>>> activeMuOR = AlloViz.State(GPCRmdID=169)

And then network computations, analyses and visualizations are performed with associated class methods. For example:

>>> activeMuOR.calculate(pkg = "corrplusCOM")
>>> activeMuOR.analyze(metrics = "btw", filterby="whole")
>>> activeMuOR.view("corrplusCOM", "btw_avg", filterby="whole")

