# gRINN: get Residue Interaction eNergies and Networks

gRINN is a software for residue interaction enery-based analysis of protein MD simulation trajectories.

This repository contains the source code. If you're only interested in using gRINN in your research, 
visit http://grinn.readthedocs.io to obtain the compiled binary for your operating system (Linux x64 or Mac OSX)

## Version

1.1.0.hf1

## Authors

Onur Serçinoğlu (onursercin@gmail.com)

Pemra Ozbek (pemra.ozbek@marmara.edu.tr)

## License

Please read the license.rst file for usage terms and conditions. 

## Dependencies

gRINN depends on either NAMD or GROMACS. Please see http://grinn.readthedocs.io/en/latest/download.html for more details.

## Availability

gRINN is available for Linux x64 and MacOSX operating systems as a standalone executable. 

## Documentation

Documentation for gRINN is located at http://grinn.readthedocs.io

## Tutorial

Best way to learn about the features of gRINN and how to use it is to follow the tutorial at 
http://grinn.readthedocs.io/en/latest/tutorial.html

## Using the source code

Note that you need either GROMACS or NAMD executable for gRINN to function.

You can just download the compiled binary of gRINN from http://grinn.readthedocs.io and follow along the Tutorial there.

If you're interested in tweaking the source code and using gRINN via your local python installation, 
beware that gRINN is dependent on several python packages (still preparing a list of dependencies).
If you know your way around python packages, then you can figure out which ones are necessary via trial and error.

Starting gRINN is straightforward. Simply type `python source/grinn.py` to start the Main GUI (after you clone/fork this repository, of course).

Accessing the calc mode from the (command-line interface) CLI: `python source/grinn.py -calc <args>`

Accessing the results GUI from the CLI: `python source/grinn.py -results <args>`

## Compiling your own gRINN

If you really want to do it, gRINN can be compiled using Pyinstaller. 
Build specification files can be found in source/build_spec folder. Using *_onefile.spec files is recommended (others might not be up-to-date)

**Please note that this is a perilous adventure!** Definitely not recommended due to a lot of library clashes between Pyinstaller, PyMol, Qt, OpenGL and matplotlib.
Getting the pyinstaller pick up all the necessary libraries can be frustrating and I can't give you one-for-all recipe for your own custom Python installation. 


## History/Change Log

### v1.1.0.hf1 (2018/06/21)
^^^^^^^^^^^^^^^^^^^^^^^

This hf (hot-fix) version fixes two bugs which rendered gRINN unusable in some cases and an addition to sample input files:

Bug fixes:

* A major bug in gRINN which leads to a failure in processing TPR files without chain IDs is corrected. gRINN will assign a default chain ID of "P" to residues which have no chain IDs assigned in input TPR.
* IEM annotation is now shown only for smallest proteins (with sizes of at most 20 amino acids).

Additions:

* charm27.ff files (used by GROMACS sample trajectory data) are included in the distribution (considering that this force-field may not be included in GROMACS installation of some users).


### v1.1.0 (2018/04/06)
^^^^^^^^^^^^^^^^^^^

This version introduces a major internal code rehaul, leaving major features of gRINN unaffected.
There are additional new features as well as minor bug fixes:

New Features:

* A new calculation setting for non-bonded interaction cutoff for NAMD simulation input is introduced. In the previous version, filtering cutoff distance parameter specified for filtering cutoff distance and non-bonded cutoff for NAMD simulation input.

* gRINN now supports Charmm simulation input as well.

Minor bug fixes:

* A bug which cause multiple parameter files reading to fail for NAMD simulation input is fixed.

* A minor bug which caused incorrect protein structure display upon start of View Results interface in Mac OS version is fixed.

### v1.0.1 (2017/12/27)

Initial release of gRINN.

## Credits

gRINN was coded in Python 2.7.

Several open-source packages, including ProDy, MDTraj, PyQt5, matplotlib, seaborn, pandas, networkx and PyMol are distributed with gRINN. More details can be found in license.rst.

A full list of credits can be found in http://grinn.readthedocs.io/en/latest/credits.html.