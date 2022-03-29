<!-- [![PyPI - Python Version](https://img.shields.io/pypi/pyversions/correlationplus)](https://pypi.org/project/correlationplus/)
[![PyPI](https://img.shields.io/pypi/v/correlationplus)](https://pypi.org/project/correlationplus/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/correlationplus/README.html)
[![Open Source License: GPL v3](https://img.shields.io/badge/License-LGPLv3-blue.svg)](https://opensource.org/licenses/LGPL-3.0)
[![Doc](https://readthedocs.org/projects/correlationplus/badge/?version=latest)](http://correlationplus.readthedocs.org/en/latest/#)
[![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/structuraldynamicslab/correlationplus/latest)](https://hub.docker.com/repository/docker/structuraldynamicslab/correlationplus)
![Conda](https://img.shields.io/conda/pn/bioconda/correlationplus)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/tekpinar/correlationplus/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/tekpinar/correlationplus) -->

# AlloViz 🔭

A Python package to interactively compute, analyze and visualize protein allosteric communication (residue interaction) networks and delta-networks.

AlloViz binds together some newly written modules with 8 Python packages that provide different ways of calculating residue interactions: [getcontacts](https://github.com/getcontacts/getcontacts), [correlationplus](https://github.com/tekpinar/correlationplus), [dynetan](https://github.com/melomcr/dynetan), [PyInteraph2](https://github.com/ELELAB/pyinteraph2), [pytraj](https://github.com/Amber-MD/pytraj), [MD-TASK](https://github.com/RUBi-ZA/MD-TASK), [gRINN](https://bitbucket.org/onursercinoglu/grinn) and [g_correlation](https://www.mpinat.mpg.de/grubmueller/g_correlation).

For the same topology and molecular dynamics (MD) trajectory, the network can be constructed based on residue contacts, movement correlation or interaction energies, depending on the package selected. Moreover, for movement correlation, the movement tracked can be that of the whole residue, its center of mass, its alpha-C or its beta-C; and it can be calculated as the Pearson's correlation coefficient, Mutual Information (MI) or Linear MI (LMI). See [below](#available-information-sources-for-network-generation).

The network can be analyzed with edge centrality metrics algorithms provided by the Python package [networkx](https://github.com/networkx/networkx), and they can be visualized in a Notebook using [nglview](https://github.com/nglviewer/nglview).

## Installation

It is recommended to use a virtual environment. This repository includes submodules that need to be appropriately cloned alongside the main repository using the `--recursive` flag. At present, virtual environment dependencies can only be installed with conda due to vmd-python not being available in PyPi.

```bash
git clone --recursive --shallow-submodules -j 6 https://github.com/frannerin/AlloViz
cd AlloViz
conda create AlloViz -c conda-forge --file environment.txt
```

Then activate the environment with `conda activate AlloViz` and install the package preferably with `pip install .`. Alternatively, use `python setup.py install`.

## Quickstart

The package is designed for use in interactive notebooks (i.e., [Jupyter](https://jupyter.org/)). The main class to be used is `State`; `Pair` allows to bind two States for delta-network analysis.

### for GPCRmd files

A State is simply defined with the GPCRmd dynamics ID number, with which files are automatically retrieved. For example:

```python
import AlloViz
activeMuOR = AlloViz.State(GPCRmdID=169)
```

And then network computations, analyses and visualizations are performed with associated class methods. For example:

```python
activeMuOR.calculate(pkg = "corrplusCOM")
activeMuOR.analyze(metrics = "btw", filterby="whole")
activeMuOR.view("corrplusCOM", "btw_avg", filterby="whole")
```

## Available information sources for network generation

<!-- https://www.tablesgenerator.com/html_tables -->

<table style="undefined;table-layout: fixed; width: 1070px">
<colgroup>
<col style="width: 334px">
<col style="width: 104px">
<col style="width: 183px">
<col style="width: 178px">
<col style="width: 70px">
</colgroup>
<thead>
  <tr>
    <th>Residue information extracted from trajectories</th>
    <th>Package</th>
    <th>Correlation measurement</th>
    <th>Subset of atoms tracked</th>
    <th>Name in AlloViz</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="2">Contact frequency</td>
    <td>getcontacts</td>
    <td>-</td>
    <td>Whole residue</td>
    <td>Getcontacts</td>
  </tr>
  <tr>
    <td rowspan="2">PyInteraph2</td>
    <td>-</td>
    <td>Whole residue</td>
    <td>PyInteraph</td>
  </tr>
  <tr>
    <td rowspan="3">Interaction energies</td>
    <td>-</td>
    <td>Whole residue</td>
    <td>PyInteraphEne</td>
  </tr>
  <tr>
    <td rowspan="2">gRINN</td>
    <td>-</td>
    <td>Whole residue</td>
    <td>GRINN</td>
  </tr>
  <tr>
    <td>Pearson's</td>
    <td>Whole residue</td>
    <td>GRINNcorr</td>
  </tr>
  <tr>
    <td rowspan="13">Movement correlation</td>
    <td rowspan="2">dynetan</td>
    <td rowspan="2">Mutual Information (MI)</td>
    <td>Whole residue</td>
    <td>Dynetan</td>
  </tr>
  <tr>
    <td>Residue COM</td>
    <td>DynetanCOM</td>
  </tr>
  <tr>
    <td rowspan="4">g_correlation</td>
    <td rowspan="2">MI</td>
    <td>alpha-C</td>
    <td>G_corrCAMI</td>
  </tr>
  <tr>
    <td>Residue COM</td>
    <td>G_corrCOMMI</td>
  </tr>
  <tr>
    <td rowspan="2">Linear MI (LMI)</td>
    <td>alpha-C</td>
    <td>G_corrCALMI</td>
  </tr>
  <tr>
    <td>Residue COM</td>
    <td>G_corrCOMLMI</td>
  </tr>
  <tr>
    <td rowspan="2">pytraj</td>
    <td rowspan="2">Pearson's</td>
    <td>alpha-C</td>
    <td>PytrajCA</td>
  </tr>
  <tr>
    <td>beta-C</td>
    <td>PytrajCB</td>
  </tr>
  <tr>
    <td>MD-TASK</td>
    <td>Pearson's</td>
    <td>alpha-C</td>
    <td>MDTASK</td>
  </tr>
  <tr>
    <td rowspan="5">correlationplus</td>
    <td rowspan="2">Pearson's</td>
    <td>alpha-C</td>
    <td>Corrplus</td>
  </tr>
  <tr>
    <td>Residue COM</td>
    <td>CorrplusCOM</td>
  </tr>
  <tr>
    <td rowspan="2">LMI</td>
    <td>alpha-C</td>
    <td>CorrplusLMI</td>
  </tr>
  <tr>
    <td>Residue COM</td>
    <td>CorrplusCOMLMI</td>
  </tr>
  <tr>
    <td rowspan="2">Dihedral correlation</td>
    <td>Pearson's</td>
    <td>Individual backbone dihedrals (Phi, psi and omega) and their combination</td>
    <td>CorrplusDihs (Corrplus[Psi, Phi, Omega])</td>
  </tr>
  <tr>
    <td>AlloViz</td>
    <td>MI</td>
    <td>Individual backbone dihedrals (Phi, psi and omega) and their combination</td>
    <td>AlloVizDihs (AlloViz[Psi, Phi, Omega])</td>
  </tr>
</tbody>
</table>

## Cite
--

## Licensing
:upside_down_face:

