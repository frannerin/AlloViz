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

AlloViz binds together some newly written modules with 8 Python packages that provide different ways of calculating residue interactions: [getcontacts](https://github.com/getcontacts/getcontacts), [correlationplus](https://github.com/tekpinar/correlationplus), [dynetan](https://github.com/melomcr/dynetan), [PyInteraph2](https://github.com/ELELAB/pyinteraph2), [pytraj](https://github.com/Amber-MD/pytraj), [MD-TASK](https://github.com/RUBi-ZA/MD-TASK), [gRINN](https://bitbucket.org/onursercinoglu/grinn) (needs [namd](https://www.ks.uiuc.edu/Research/namd/)) and [MDEntropy](https://github.com/msmbuilder/mdentropy).

<!-- [g_correlation](https://www.mpinat.mpg.de/grubmueller/g_correlation), [GSAtools](https://github.com/AllosterIt/GSAtools) -->

For the same topology and molecular dynamics (MD) trajectory, the network can be constructed based on residue contact frequency, correlation of atom movement or dihedrals, or interaction energies, depending on the package selected. Moreover, for example for movement correlation, the movement tracked can be that of the whole residue, its center of mass, its alpha-C or its beta-C; and it can be calculated as the Pearson's correlation coefficient, Mutual Information (MI) or Linear MI (LMI). See [below](#available-information-sources-for-network-generation).

The resulting network can be analyzed with edge centrality metrics algorithms provided by the Python package [networkx](https://github.com/networkx/networkx), and they can be visualized in an interactive Python Notebook (i.e., [Jupyter](https://jupyter.org/)) using [nglview](https://github.com/nglviewer/nglview).

## Installation

It is recommended to use a virtual environment ([Miniconda](https://docs.conda.io/en/latest/miniconda.html)). This repository includes submodules that need to be appropriately cloned alongside the main repository using the `--recursive` flag. At present, virtual environment dependencies can only be correctly installed with conda.

```bash
git clone --recursive --shallow-submodules -j 9 https://github.com/frannerin/AlloViz
conda create AlloViz --file AlloViz/conda_explicit.txt
conda activate AlloViz
```

Then go to the package folder (`cd AlloViz`) and install the package, preferably with `pip install .`.
<!--`pip install .`. Alternatively, use-->

## Quickstart

The package is designed for use in interactive notebooks. The main class to be used is `Protein`; `Pair` allows to bind two States for delta-network analysis.

### for GPCRmd files

A Protein is simply defined with the GPCRmd dynamics ID number, with which files are automatically retrieved from the database. For example:

```python
import AlloViz
activeMuOR = AlloViz.Protein(GPCR=169)
```

And then network computations, analyses and visualizations are performed with associated class methods. For example:

```python
activeMuOR.calculate(pkg = "corrplusCOM")
activeMuOR.analyze(metrics = "btw", filterby="whole")
activeMuOR.view("corrplusCOM", "btw_avg", filterby="whole")
```

## Available information sources for network generation

<!-- https://www.tablesgenerator.com/html_tables
https://github.com/msmbuilder/msmbuilder/blob/515fd5c27836c797692d600216b5eb224dfc1c5d/msmbuilder/featurizer/featurizer.py#L802
 -->

<table>
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
    <td rowspan="9">Movement correlation</td>
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
    <td rowspan="4">Dihedral correlation</td>
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
  <tr>
    <td rowspan="3">MDEntropy</td>
    <td rowspan="3">MI</td>
    <td>Combination of the backbone dihedrals (Phi, psi and omega)</td>
    <td>MDEntropyDihs</td>
  </tr>
  <tr>
    <td><a href="https://github.com/msmbuilder/msmbuilder/blob/515fd5c27836c797692d600216b5eb224dfc1c5d/msmbuilder/featurizer/featurizer.py#L802" target="_blank" rel="noopener noreferrer">Alpha angle</a> (dihedral between i-1, i, i+1 and i+2's alpha-Cs)</td>
    <td>MDEntropyAlphaAngle</td>
  </tr>
  <tr>
    <td rowspan="3">Contact frequency<br></td>
    <td>Whole residue</td>
    <td>MDEntropyContacts</td>
  </tr>
  <tr>
    <td>getcontacts</td>
    <td>-</td>
    <td>Whole residue</td>
    <td>Getcontacts</td>
  </tr>
  <tr>
    <td rowspan="2">PyInteraph2</td>
    <td rowspan="2">-</td>
    <td>Whole residue</td>
    <td>PyInteraph</td>
  </tr>
  <tr>
    <td rowspan="3">Interaction energies</td>
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
</tbody>
</table>

## Cite
--

## Licensing
:upside_down_face:

