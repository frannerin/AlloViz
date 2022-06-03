Introduction
============

Welcome to gRINN (get Residue Interaction Energies and Networks).

gRINN is a software tool for residue interaction-energy based investigation of protein Molecular Dynamics simulations.

Main functionality includes:

* Calculation of pairwise amino acid non-bonded interaction energies from NAMD or GROMACS generated Molecular Dynamics (MD) simulation trajectories by interoperating with NAMD/GMX simulation engines.

* Equal-time linear correlations between interaction energy time series.

* Custom residue selections & multiple processor usage.

* A visualization interface for:

	* Viewing pairwise interaction energies and their correlations alongside an embedded PyMol molecular viewer.
	* Protein Energy Network construction and visualization of simple residue-based local network metrics (Degrees, Betweenness Centrality and Closeness Centrality) 
	* Shortest path analysis.