"""Main AlloViz classes: `Protein` and `Delta`

The :class:`~AlloViz.Protein` class is AlloViz's main class, processing a structure and
trajectory(ies) input files and allowing to calculate, analyze and visualize the
allosteric communication networks with its associated methods.

The :class:`~AlloViz.Delta` class takes two :class:`~AlloViz.Protein` objects as input to
calculate the delta-network to highlight the differences in the allosteric communication
between two systems.
   
Notes
=====
This module is not meant to be used directly, as the classes are imported in the
namespace of AlloViz itself.

"""

import os, io, re

from multiprocess import Pool
import MDAnalysis as mda
import numpy as np

# from .Analysis import Whole, Incontact, Intercontact
#from .Filtering import Filtering#Whole, Incontact, Intercontact
# from .Visualization import Edges, Nodes
from . import Analysis
from . import Elements# import Edges, Nodes


from .utils import rgetattr, rhasattr
from . import utils
from . import trajutils

from .. import Wrappers


class _Store:
    r"""Private class that is only used to create empty class' attributes to which other
    attributes are added.
    """

    pass


class Protein:
    r"""AlloViz main class.

    Objects of the class process the input structure and trajectory files and allow to
    :meth:`~AlloViz.Protein.calculate` the different networks, which are stored in
    classes of the :mod:`AlloViz.Wrappers` module and can be filtered, analyzed and
    visualized with subsequent methods.

    Parameters
    ----------
    pdb : str
        Filename of the PDB structure to read, process and use.
    trajs : str or list
        Filename(s) of the MD trajectory (or trajectories) to read and use. File format
        must be recognized by MDAnalysis (e.g., xtc).
    GPCR : bool or int, default: False
        Use `True` if the structure is a GPCR, or pass the ID of a GPCRmd database
        dynamics entry, without specifying the `pdb` nor the `trajs` parameters, to
        automatically retrieve the files from the database and process them. They will
        be downloaded to `path`, which if left undefined will default to the GPCRmd ID.
    name : str, default: "protein"
        Name string to be used, e.g., for the title of the colorbar shown when
        representing a network with nglviewer.
    path : str, optional
        Path to store results in. It can exist already or not, and a new folder inside
        it called `data` will be created to store computation results and other
        files. If unspecified, it defaults to "." or the GPCRmd ID in case it is used.
    protein_sel : str, default: :attr:`AlloViz.Protein._protein_sel`
        MDAnalysis atom selection string to select the protein structure from the
        Universe (e.g., in case simulations in biological conditions are used, to avoid
        selecting extra chains, water molecules, ions...). It defaults to "(same segid as
        protein) and (not segid LIG) and (not chainid L)" and it can be extended using,
        e.g., :attr:`AlloViz.Protein._protein_sel` `+ " and {customselection}"`.

    Other Parameters
    ----------------
    psf : str, optional
        Optional kwarg of the filename of the .psf file corresponding to the pdb used.
        It is required alongside the `parameters` file to use the gRINN network
        construction method.
    parameters : str, optional
        Optional kwarg of the filename of the MD simulation force-field parameters file
        in NAMD format. It is required alongside the `psf` file to use the gRINN network
        construction method.
    special_res : dict, optional
        Optional kwarg of a dictionary containing a mapping of special residue 3/4-letter
        code(s) present in the structure to the corresponding standard 1-letter code(s).

    Attributes
    ----------
    pdb
    trajs
    GPCR
    protein : :class:`~MDAnalysis.core.universe.Universe`
        Universe of the processed pdb with only the selected `protein_sel` atoms.
    u : :class:`~MDAnalysis.core.universe.Universe`
        Universe of the processed pdb and trajectory files with only the `protein_sel`
        atoms.

    Raises
    ------
    FileNotFoundError
        If any of the files passed in the `pdb` and `traj` (and `psf` and `parameters`
        if provided) parameters cannot be accessed.

    See Also
    --------
    AlloViz.Delta : Class for calculation of the delta-network(s) between two Protein
                    objects.

    Examples
    --------
    >>> opioidGPCR = AlloViz.Protein(GPCR=169)
    >>> print(opioidGPCR.u)
    <Universe with 88651 atoms>
    """

    _protein_sel = "(same segid as protein) and (not segid LIG) and (not chainid L)"
    """
    Class attribute used to select only protein atoms from the input files with
    :external:ref:`MDAnalysis selection syntax <selection-commands-label>`. It defaults
    to "(same segid as protein) and (not segid LIG) and (not chainid L)" and it can be
    extended through the `protein_sel` parameter using, e.g.,
    :attr:`AlloViz.Protein._protein_sel` + " and {customselection}"
    """

    def __init__(
        self,
        pdb="",
        trajs=[],
        GPCR=False,
        name=None,
        path=None,
        protein_sel=None,
        **kwargs,
    ):
        self.GPCR = GPCR

        # If a GPCRmd ID is passed
        if not isinstance(self.GPCR, bool) and isinstance(self.GPCR, int):
            self.name = f"{self.GPCR}" if not name else name
            self._path = f"{self.GPCR}" if not path else path
            os.makedirs(self._path, exist_ok=True)

            # Download the files from GPCRmd
            if not any(
                [
                    re.search("(pdb$|psf$|xtc$|dcd$|parameters$)", file)
                    for file in os.listdir(self._path)
                ]
            ):
                trajutils.download_GPCRmd_files(self.GPCR, self._path)
            files = os.listdir(self._path)

            # Retrieve filenames from the files downloaded into self._path
            get_filename = (
                lambda ext: self._path
                + "/"
                + next(file for file in files if re.search(f"{ext}$", file))
            )
            self.pdb = get_filename("pdb")
            self.trajs = list(
                sorted(
                    f"{self._path}/{traj}"
                    for traj in files
                    if re.search("^(?!\.).*\.(xtc|dcd)$", traj)
                )
            )
            self.psf = get_filename("psf")
            self._paramf = get_filename("parameters")

        # If pdb and trajectory files are passed
        else:
            self.name = pdb.replace(".pdb", "").split("/")[-1] if not name else name
            self._path = "." if not path else path
            os.makedirs(self._path, exist_ok=True)

            self.pdb = pdb
            self.trajs = trajs if isinstance(trajs, list) else [trajs]

            # Check if a .psf and a force-field parameters files have been passed for gRINN network calculation
            passed_psf_params = "parameters" in kwargs and "psf" in kwargs
            if passed_psf_params:
                self.psf = kwargs["psf"]
                self._paramf = kwargs["parameters"]

            # Check if all the filehandles passed exist, and raise an error if not
            files_to_check = (
                self.trajs + [self.pdb]
                if not passed_psf_params
                else self.trajs + [self.pdb, self.psf, self._paramf]
            )
            files_exist = {file: os.path.isfile(file) for file in files_to_check}
            if any([not file_exist for file_exist in files_exist.values()]):
                raise FileNotFoundError(
                    f"Some of the files could not be found: {files_exist}"
                )

        self._protein_sel = Protein._protein_sel if not protein_sel else protein_sel

        # Set the names of the directories and files that will be creating when processing the input files
        self._datadir = f"{self._path}/data"
        os.makedirs(self._datadir, exist_ok=True)

        self._pdbf = f"{self._datadir}/protein.pdb"
        self._trajs = dict(
            [
                (num + 1, f"{self._datadir}/traj_{num+1}.xtc")
                for num in range(len(self.trajs))
            ]
        )
        self._psff = self._pdbf.replace("pdb", "psf") if hasattr(self, "psf") else None

        # Names of the directories and files of the future pdb and trajectory(ies) of the residues' Center Of Mass
        compath = f"{self._datadir}/COM_trajs"
        os.makedirs(compath, exist_ok=True)
        self._compdbf = f"{compath}/ca.pdb"
        self._comtrajs = {num: f"{compath}/{num}.xtc" for num in self._trajs}

        # If the processed filenames don't exist yet as files, process the input files; if special_res kwarg is used it will be passed
        if any(
            [
                not os.path.isfile(f)
                for f in list(self._trajs.values())
                + list(self._comtrajs.values())
                + [self._pdbf]
            ]
        ):
            self._pdbf, self._trajs, self._psff = trajutils.process_input(
                self.GPCR,
                self.pdb,
                self._protein_sel,
                self._pdbf,
                self._compdbf,
                rgetattr(self, "psf"),
                self._psff,
                self.trajs,
                self._trajs,
                self._comtrajs,
                # **kwargs,
                **{"special_res": kwargs["special_res"]} if "special_res" in kwargs else {},
            )

        # Set the protein/pdb and trajectory(ies) MDAnalysis' Universes with the processed files
        self.protein = mda.Universe(self._pdbf)
        """:class:`~MDAnalysis.core.universe.Universe`"""
        self.u = mda.Universe(self._pdbf, *list(self._trajs.values()))
        """:class:`~MDAnalysis.core.universe.Universe`"""

        # Bonded cysteines must be identified to remove them from non-covalent contacts calculations (i.e., PyInteraph2_Contacts)
        self._bonded_cys = trajutils.get_bonded_cys(self._pdbf)

        # _dihedral_residx returns a list of indices of all the protein's residues in the Universes that are not in the extremes of the chain(s)
        # These will be the residues for which Dihedral correlations can be calculated (residues in the extremes have some dihedrals missing)
        # MDEntropy_AlphaAngle uses information of residues i-1, i, i+1, i+2; so end=-2 is passed to also exclude the second-to-last residue of the chain(s)
        _res_arrays = np.split(
            self.protein.residues.resindices,
            np.where(np.diff(self.protein.residues.resnums) != 1)[0] + 1,
        )
        self._dihedral_residx = lambda end=-1: [
            elem for arr in _res_arrays for elem in arr[1:end]
        ]

        # _translate_ix returns a translated DataFrame Index (string) or MultiIndex (tuple) element according to the mapper dictionary passed
        # The translation function _translate_ix(mapper) is passed as the function used by the .map method of a DataFrame Index or MultiIndex to translate it
        self._translate_ix = (
            lambda mapper: lambda ix: tuple(mapper[_] for _ in ix)
            if isinstance(ix, tuple)
            else mapper[ix]
        )

    def __sub__(self, other):
        """
        The result of subtracting an "other" Protein object from "self" is a new object
        (delta) which has the attributes that the Proteins have in common (regarding
        network construction methods/packages, filtering schemes and elements-analysis),
        whose values arise from subtracting the corresponding dataframes (i.e.,
        subtracting the networks) to create a delta-network.
        """
        # Create an empty class (_Store dummy class) instance to store the results
        delta = _Store()

        # For each package that the two Proteins have in common; taking them from those attributes in self's __dict__ that are valid package names
        for pkg in (
            key for key in self.__dict__ if utils.pkgname(key, fail=False) and key in other.__dict__
        ):
            # Set an attribute with the package's name in the empty object
            setattr(delta, pkg, _Store())

            # For each filtering scheme (taken from self's __dict__ with a valid name) for which results have been analyze in "self" and that are also in "other"
            for filtering in (
                key
                for key in getattr(self, pkg).__dict__
                if any([f in key for f in utils.filteringsl])
                and key in getattr(other, pkg).__dict__
            ):                
                # Also set an attribute with the filtering name
                setattr(getattr(delta, pkg), filtering, _Store())

                # And finally for each Element (taken from self's __dict__ with a valid name) that is also in "other"
                for elem in (
                    key
                    for key in rgetattr(self, pkg, filtering).__dict__
                    if key.lower() in ["nodes", "edges"]
                    and key in rgetattr(other, pkg, filtering).__dict__
                ):
                    # Calculate the Elements' difference exploiting their custom subtraction special method definition
                    dif = rgetattr(self, pkg, filtering, elem) - rgetattr(
                        other, pkg, filtering, elem
                    )
                    # And save the result as a new Element object (to exploit its view method)
                    elemclass = eval(f"Elements.{elem.capitalize()}")
                    setattr(
                        rgetattr(delta, pkg, filtering),
                        elem,
                        elemclass(dif),
                    )
                    rgetattr(delta, pkg, filtering, elem)._parent = self._delta

        return delta.__dict__

    def calculate(self, pkgs="all", cores=1, **kwargs):
        r"""Calculate rwa edge weights of allosteric networks.

        Send the computation of the raw edge weights for the selected network
        construction methods.

        Parameters
        ----------
        pkgs : str or list, default: "all"
            Package(s)/Network construction method(s) for which to send raw edge weight
            computation. "all" sends the computation for all available methods within
            AlloViz (check :data:`AlloViz.AlloViz.utils.pkgsl`).
        cores : int, default: 1
            Number of cores to use for parallelization with a `multiprocess` Pool.
            Default value only uses 1 core with a custom :class:`AlloViz.utils.dummypool`
            that performs computations synchronously.

        Other Parameters
        ----------------
        taskcpus : int, optional
            Optional kwarg to specify the amount of cores that parallelizable network
            construction methods can use (i.e., AlloViz's method, getcontacts, dynetan,
            PyInteraph, MDEntropy and gRINN).
        stride : int, optional
            Optional kwarg to specify the striding to be done on the trajectory(ies) for
            computationally-expensive network construction methods: i.e., dynetan and
            AlloViz's own method. Default is no striding.
        namd : str, optional
            Optional kwarg pointing to the namd2 executable location; if the `namd`
            command is accessible through the CLI it is automatically retrieved with the
            `distutils` package.
        GetContacts_threshold : float, optional
            Optional kwarg to specify the minimum contact frequency (0-1, default 0)
            threshold for GetContacts results, which will be used to filter out contacts
            with a frequency (average) lower than it.
        chis : int, optional
            Optional kwarg to specify the number of side-chain chi dihedral angles (up to
            5) to combine when sending the calculation of a child of the Combined_Dihs
            Wrappers' base class that includes chi dihedrals in its calculation.
        MDEntropy_method : str, optional, {"knn", "grassberger", "chaowangjost"}
            Optional kwarg to specify the method to calculate the entropy of the
            variables for Mutual Information estimation when using one of the
            MDEntropy network construction methods (default: "grassberger").

        See Also
        --------
        AlloViz.Wrappers.Base.Base : Base class to launch and store calculation results.
        AlloViz.Protein.filter : Class method to filter the network raw edge
                                 weights with different criteria.
        AlloViz.Protein.analyze : Class method to analyze the calculated raw edge weights
                                  with graph theory-based methods.
        AlloViz.Protein.view : Class method to visualize the network on the protein
                               structure.

        Notes
        -----
        Calculation results are stored as new instance attributes with the same name as
        the selected packages/network construction methods. If the object is created
        providing more than one trajectory file, the average and standard error of the
        weights between the replicas are also calculated.

        Examples
        --------
        If we have a 6-core computer and want to use 2 cores to compute the values for
        each of the three trajectories of the protein (e.g., GPCRmd stores three replicas
        of each structure):

        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate("dynetan", cores=6, taskcpus=2)
        >>> print(opioidGPCR.dynetan.raw.shape)
        (41041, 5)
        """
        # Calculate for "all" packages or the ones passed as parameter (check that they are on the list of available packages and retrieve their case-sensitive names, else raise an Exception)
        pkgs = utils.make_list(pkgs, if_all=utils.pkgsl, apply=utils.pkgname)
        combined_dihs = [pkg for pkg in pkgs if "Dihs" in pkg]
        
        # Add the calculations of the individual dihedrals if a combined_dihedral has been passed
        bb = ["Phi", "Psi"]
        sc = [f"Chi{i+1}" for i in range(4)]
        if len(combined_dihs) > 0:
            for comb in combined_dihs:
                pkg = comb.split("_")[0]
                if pkg == "CARDS":
                    if "Sidechain" in comb or "Backbone" in comb:
                        pkg = comb.rsplit("_", 3)[0]
                    else:
                        pkg = comb.rsplit("_", 2)[0]
                dihs = bb if "Backbone" in comb else sc if "Sidechain" in comb else bb+sc
            pkgs += [f"{pkg}_{dih}" for dih in dihs]

        # Objects from the classes in the Wrappers module need to be passed a dictionary "d" containing all the attributes of the source Protein object and the passed kwargs
        d = self.__dict__.copy()
        d.update(kwargs)

        # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
        # Changing it inside the `utils` module allows to share the same one between modules
        if cores > 1:
            mypool = Pool(cores)
        else:
            mypool = utils.dummypool()
        utils.pool = mypool
        print(utils.pool)
            
        if any(["CARDS" in pkg for pkg in pkgs]):
            Wrappers.CARDS_w.CARDS(self, d)

        for pkg in set(pkgs) - set(combined_dihs):
            # Establish the corresponding Wrappers' class
            pkgclass = eval(f"Wrappers.{utils.pkgname(pkg)}")

            # Setting the class as a new attribute will initialize all calculations asynchronously (synchronously if a dummypool is used)
            if not hasattr(self, pkgclass.__name__):
                setattr(self, pkgclass.__name__, pkgclass(self, d))
        
        # if cores > 1:
        # Close the pool
        mypool.close()
        mypool.join()
        mypool = utils.dummypool()
        
        if len(combined_dihs) > 0:
            # Calculate now the combination of dihedrals, which is just a combination of the already-calculated data
            for pkg in combined_dihs:
                pkgclass = eval(f"Wrappers.{utils.pkgname(pkg)}")
                if not hasattr(self, pkgclass.__name__):
                    setattr(self, pkgclass.__name__, pkgclass(self, d))
        
        return getattr(self, pkgclass.__name__) if len(pkgs) == 1 else None
    
    
    
    def filter(self, pkgs="all", filterings="all", **kwargs):
        r"""Filter network edges
        
        Filter the networks according to the selected criteria to perform analyses on
        (all or) a subset of the edges. It calls
        :meth:`AlloViz.Wrappers.Base.Base.filter` and results are stored in instances
        of the :class:`AlloViz.AlloViz.Filtering.Filtering` class. The different filtering
        options are detailed in the :mod:`~AlloViz.AlloViz.Filtering` module.

        Parameters
        ----------
        pkgs : str or list, default: "all"
            Package(s)/Network construction method(s) for which to analyze their raw edge
            weights, which must be already calculated and their data saved as instance
            attribute. In this case, "all" sends the computation for all available
            methodsthat are already calculated and saved as instance attributes.
        filterings : str or list of strs and/or lists, default: "all"
            Filtering scheme(s) with which to filter the list of network edges before
            analysis. It can be a string, or a list of strings and/or lists. A list of
            lists (also with or without strings) is used to filter with a combination of
            criteria. All available (and combinable) filtering options are functions in
            the :mod:`~AlloViz.AlloViz.Filtering` module: 
            :func:`~AlloViz.AlloViz.Filtering.All`,
            :func:`~AlloViz.AlloViz.Filtering.GetContacts_edges`,
            :func:`~AlloViz.AlloViz.Filtering.No_Sequence_Neighbors`,
            :func:`~AlloViz.AlloViz.Filtering.GPCR_Interhelix`. The default "all"
            performs all the available filtering schemes (no combinations).
        
        Other Parameters
        ----------------
        GetContacts_threshold : float
            Optional kwarg that can be passed to specify the minimum contact frequency
            (0-1, default 0) threshold, which will be used to filter out contacts with a
            frequency (average) lower than it before analysis.
        Sequence_Neighbor_distance : int
            Optional kwarg that can be passed to specify the minimum number of sequence
            positions/distance between residues of a pair to retain in
            `No_Sequence_Neighbors` filtering, which defaults to 5.
        Interresidue_distance : int or float
            Optional kwarg that can be passed to specify the minimum number of angstroms
            that the CA atoms of residue pairs should have between each other in the initial
            PDB/structure (default 10 Å) to be considered spatially distant.

        See Also
        --------
        AlloViz.AlloViz.Filtering.Filtering : Filtering class.
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.analyze : Class method to analyze the calculated raw edge weights
                                  with graph theory-based methods.
        AlloViz.Protein.view : Class method to visualize the network on the protein
                               structure.
        
        Examples
        --------
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate(["getcontacts", "dynetan"], cores=6, taskcpus=2)
        >>> opioidGPCR.filter("dynetan", ["GetContacts_edges", ["GetContacts_edges", "GPCR_Interhelix"]])
        >>> opioidGPCR.dynetan.GetContacts_edges_GPCR_Interhelix
        <AlloViz.AlloViz.Filtering.Filtering at 0x7f892c3c0fa0>
        """
        # Filter "all" packages (all the available packages that have been previously calculated and are in __dict__)
        # or the ones passed as parameter (check that they are on the list of available packages and retrieve their case-sensitive names, else raise an Exception)
        pkgs = utils.make_list(pkgs, if_all = [key for key in utils.pkgsl if key in self.__dict__], apply = utils.pkgname)
        
        for pkgn in pkgs:
            pkg = rgetattr(self, pkgn)
            if not pkg:
                print(f"{pkgn} calculation results are needed first")
                continue
            result = pkg.filter(filterings, **kwargs)
            
        return result if (len(pkgs) == 1) else None
    
    def analyze(self, pkgs="all", filterings="all", elements="edges", metrics="all", normalize=True, cores=1, **kwargs):
        r"""Analyzed filtered network
        
        Analyze the selected (un)filtered networks with the passed elements-metrics. It
        calls :meth:`AlloViz.AlloViz.Analysis.analyze` and results are stored in
        instances of classes from the :mod:`AlloViz.AlloViz.Elements` module, which
        extend the :class:`pandas.DataFrame` class.

        Parameters
        ----------
        pkgs : str or list, default: "all"
            Package(s)/Network construction method(s) for which to analyze their raw edge
            weights, which must be already calculated and their data saved as instance
            attribute. In this case, "all" sends the computation for all available
            methods that are already calculated and saved as instance attributes.
        filterings : str or list, default: "all"
            Filtering scheme(s) for which to perform the analyses, which must exist
            already for the selected packages. "all" sends the computation for all
            available schemes that are already saved.
        elements : str or list, {"edges", "nodes"}
            Network elements for which to perform the analysis.
        metrics : str or list, default: "all"
            Network metrics to compute, which must be keys in the `nodes_dict` or
            `edges_dict` dictionaries. Default is "all" and it sends the computation for
            all the metrics defined in the corresponding dictionary of the selected
            elements in `element`.
        normalize : bool, default: True
            Passed to the NetworkX functions that calculate the metrics, to output
            normalized results or not.
        cores : int, default: 1
            Number of cores to use for parallelization with a `multiprocess` Pool.
            Default value only uses 1 core with a custom :class:`AlloViz.utils.dummypool`
            that performs computations synchronously.

        Other Parameters
        ----------------
        nodes_dict, edges_dict : dict
            Optional kwarg(s) of the dictionary(ies) that maps network metrics custom names
            (e.g., betweenness centrality, "btw") with their corresponding NetworkX
            function (e.g., "networkx.algorithms.centrality.betweenness_centrality").
            Functions strings must be written as if they were absolute imports, and must
            return a dictionary of edges or nodes, depending on the element dictionary in
            which they are. The keys of the dictionaries will be used to name the columns
            of the analyzed data that the functions produce. Defaults are
            :data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
            :data:`~AlloViz.AlloViz.Analysis.edges_dict`.

        See Also
        --------
        AlloViz.AlloViz.Analysis : Module with analysis functions.
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.filter : Class method to filter the network(s) raw edge
                                 weights with different criteria.
        AlloViz.Protein.view : Class method to visualize the network on the protein
                               structure.

        Examples
        --------
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate(["getcontacts", "dynetan"], cores=6, taskcpus=2)
        >>> opioidGPCR.filter("dynetan", "GetContacts_edges")
        >>> opioidGPCR.analyze("dynetan", "GetContacts_edges", "nodes", "btw")
        <AlloViz.AlloViz.Elements.Nodes at 0x7f892c3c0fa0>
        """
        # Analyze "all" packages (all the available packages that have been previously calculated and are in __dict__)
        # or the ones passed as parameter (check that they are on the list of available packages and retrieve their case-sensitive names, else raise an Exception)
        pkgs = utils.make_list(pkgs, if_all = [key for key in utils.pkgsl if key in self.__dict__], apply = utils.pkgname)
        
        # # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
        # # Changing it inside the `utils` module allows to share the same one between modules
        # if cores > 1:
        #     mypool = Pool(cores)
        # else:
        #     mypool = utils.dummypool()
        # utils.pool = mypool
        # print(utils.pool)       
        if cores > 1:
            mypool = Pool(cores)
        else:
            mypool = utils.dummypool()
        utils.pool = mypool
        print(utils.pool)   
        
        for pkg in pkgs:
            # Analyze for "all" filterings or the ones passed as parameter
            # "all": all filterings that have been previously done and are in the pkg's __dict__ and also contain at least one of the utils.filteringsl members (to include combinations)
            filterings = utils.make_list(
                filterings,
                if_all = [filt for filt in rgetattr(self, utils.pkgname(pkg), "__dict__") if any([f in filt for f in utils.filteringsl])]
            )
            
            for filtering in filterings:
                filtered = rgetattr(self, utils.pkgname(pkg), filtering)
                
                # if not filtered:
                #     print(f"{utils.pkgname(pkg)} {filtering} results are needed first")
                #     continue
                # elif filtered._filtdata.size == 0:
                #     print(f"{utils.pkgname(pkg)} {filtering} is not a connected network (or subnetwork)")
                #     continue
                # # result = 
                # Analysis.analyze(filtered, elements, metrics, normalize, **kwargs)
                filtered.analyze(elements, metrics, normalize, cores=1, **kwargs)
                
        # # Close the pool
        # mypool.close()
        # mypool.join()
        # mypool = utils.dummypool()
        # if cores > 1:
        # Close the pool
        mypool.close()
        mypool.join()
        mypool = utils.dummypool()
            
        #return #result if (len(pkgs) == 1 and len(filterings) == 1) else None
    

#     def analyze(
#         self,
#         pkg="all",
#         filterby="Whole",
#         element="edges",
#         metrics="all",
#         normalize=True,
#         cores=1,
#         **kwargs,
#     ):
#         r"""Analyze calculated edge weights with network analyses.

#         Perform analyses of the raw edge weights for the selected packages/network
#         construction methods, filtering schemes and network elements, calculating the
#         desired network metrics.

#         Parameters
#         ----------
#         pkg : str or list, default: "all"
#             Package(s)/Network construction method(s) for which to analyze their raw edge
#             weights, which must be already calculated and their data saved as instance
#             attribute. In this case, "all" sends the computation for all available
#             methodsthat are already calculated and saved as instance attributes.
#         filterby : str or list, {"Whole", "Incontact", "Intercontact"}
#             Filtering schemes with which to filter the list of network edges before
#             analysis. "Whole" is the default and means no filtering, and "Incontact" and
#             "Intercontact" (see Notes) require that raw data from the `getcontacts`
#             package has been calculated first.
#         element : str or list, {"edges", "nodes"}
#             Network elements for which to perform the analysis.
#         metrics : str or list, default: "all"
#             Network metrics to compute, which must be keys in the `nodes_dict` or
#             `edges_dict` dictionaries. Default is "all" and it sends the computation for
#             all the metrics defined in the corresponding dictionary of the selected
#             elements in `element`.
#         normalize : bool, default: True
#             Passed to the NetworkX functions that calculate the metrics, to output
#             normalized results or not.
#         cores : int, default: 1
#             Number of cores to use for parallelization with a `multiprocess` Pool.
#             Default value only uses 1 core with a custom :class:`AlloViz.utils.dummypool`
#             that performs computations synchronously.

#         Other Parameters
#         ----------------
#         nodes_dict, edges_dict : dict
#             Optional kwarg(s) of the dictionary(ies) that maps network metrics custom names
#             (e.g., betweenness centrality, "btw") with their corresponding NetworkX
#             function (e.g., "networkx.algorithms.centrality.betweenness_centrality").
#             Functions strings must be written as if they were absolute imports, and must
#             return a dictionary of edges or nodes, depending on the element dictionary in
#             which they are. The keys of the dictionaries will be used to name the columns
#             of the analyzed data that the functions produce. Defaults are
#             :data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
#             :data:`~AlloViz.AlloViz.Analysis.edges_dict`.
#         **kwargs
#             `GetContacts_threshold` kwarg can be passed to specify the minimum contact
#             frequency (0-1, default 0) threshold, which will be used to filter out
#             contacts with a frequency (average) lower than it before analysis. Make sure
#             to delete/have deleted all previous analysis attributes and files of any
#             network construction method. `Intercontact_distance` kwarg can be passed to
#             specify the minimum number of sequence positions/distance between residues of
#             a pair to retain in Intercontact filtering, which defaults to 5.


#         See Also
#         --------
#         AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
#                                     weights with different network construction methods.
#         AlloViz.Protein.view : Class method to visualize the network on the protein
#                                structure.

#         Notes
#         -----
#         Method returns nothing, but analysis results are stored as nested attributes
#         "inside" each of the packages' attributes of the :class:`~AlloViz.Protein` object,
#         first using the name of the filtering scheme (e.g., `.Package.filterby`) and
#         lastly with the analyzed network element name (e.g., `.Package.filterby.element`).
#         If the `Protein` object was created providing more than one trajectory file, the
#         analyses are performed both on the replicas' weights and the average, and an
#         average and standard error of the replicas' analysis results are also calculated.

#         "Incontact" filtering only retains edges of residue pairs in contact -those for
#         which GetContacts is able to compute contact frequencies-, and "Intercontact"
#         only keeps edges of pairs that are both in contact and apart in the sequence (more
#         than 5 positions away in the sequence).

#         Examples
#         --------
#         >>> opioidGPCR = AlloViz.Protein(GPCR=169)
#         >>> opioidGPCR.calculate(["getcontacts", "dynetan"], cores=6, taskcpus=2)
#         >>> opioidGPCR.analyze("dynetan", filterby="Intercontact", element=["edges", "nodes"], metrics="btw")
#         >>> print(opioidGPCR.dynetan.Intercontact.edges.df.shape)
#         (3410, 5)
#         """
#         # Calculate for "all" packages (all the available packages that have been previously calculated and are in __dict__)
#         # or the ones passed as parameter (check that they are on the list of available packages and retrieve their case-sensitive names, else raise an Exception)
#         pkgs = (
#             [pkg for pkg in self.__dict__ if pkg in utils.pkgsl]
#             if pkg == "all"
#             else [utils.pkgname(p) for p in pkg]
#             if isinstance(pkg, list)
#             else [utils.pkgname(pkg)]
#         )
#         # Calculate for the passed Filterings
#         filterbys = (
#             utils.filterbysl
#             if filterby == "all"
#             else filterby
#             if isinstance(filterby, list)
#             else [filterby]
#         )
#         # And for the passed Elements, and their corresponding network Metrics if they are defined in the dictionaries
#         elements = element if isinstance(element, list) else [element]
#         metrics = (
#             set(list(nodes_dict.keys()) + list(edges_dict.keys()))
#             if metrics == "all"
#             else metrics
#             if isinstance(metrics, list)
#             else [metrics]
#         )

#         # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
#         # Changing it inside the utils module allows to share the same one between modules
#         if cores > 1:
#             mypool = Pool(cores)
#         else:
#             mypool = utils.dummypool()
#         utils.pool = mypool
#         print(utils.pool)

#         # For each package to be analyzed
#         for pkgn in pkgs:
#             pkg = rgetattr(self, pkgn)
#             if not pkg:
#                 print(f"{pkgn} calculation results are needed first")
#                 continue

#             # Add (or retrieve) the passed Filterings attributes, and add the passed Elements-Metrics to them
#             for filterby in filterbys:
#                 # Establish the corresponding filterby/Analysis' class and check if it is already an attribute or set it otherwise
#                 anaclass = (
#                     eval(filterby.capitalize())
#                     if isinstance(filterby, str)
#                     else filterby
#                 )
#                 get_ana = lambda: rgetattr(pkg, anaclass.__name__)
#                 if not get_ana():
#                     setattr(
#                         pkg,
#                         anaclass.__name__,
#                         anaclass(pkg, elements, metrics, normalize, **kwargs),
#                     )
#                 else:
#                     # Add the passed Elements-Metrics
#                     get_ana().add_metrics(elements, metrics, normalize, **kwargs)

#         # Close the pool
#         mypool.close()
#         mypool.join()

    def view(
        self,
        pkg,
        metric,
        filtering="Whole",
        element="edges",
        num: int = 20,
        colors: list = ["orange", "turquoise"],
        nv=None,
    ):
        r"""Visualize the analyzed networks on the protein structure.

        Visualize the network corresponding to the selected package/network construction
        method, filtering scheme, network element and network metric or metrics, all of
        which must be previously calculated and analyzed. The number of (each of the)
        elements to show and the color scale can be specified. It calls
        :meth:`AlloViz.Elements.Element.view`.

        Parameters
        ----------
        pkg : str
            Package/Network construction method for which to show the network.
        metric : str
            Network metric for which to show the network.
        filtering : str, {"Whole", "Incontact", "Intercontact"}
            Filtering scheme for which to show the network.
        element : str or list, {"edges", "nodes"}
            Network element or elements to show on the protein structure representing the
            chosen package, filtering scheme and metric.
        num : int, default: 20
            Number of (each of the) network elements to show on the structure.
        colors : list, default: ["orange", "turquoise"]
            List of two colors to assign to the minimum and maximum values of the network
            to be represented, respectively. Middle value is assigned "white" and it will
            be the mean of the network values or 0 if the network has both negative and
            positive values.
        nv : :class:`nglview.NGLWidget`, optional
            A structure representation into which the shapes representing the chosen
            network elements will be added.

        Returns
        -------
        :class:`nglview.NGLWidget`

        See Also
        --------
        AlloViz.AlloViz.Elements : Module with classes for storage and visualization of
                                   filtered and analyzed networks.
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.filter : Class method to filter the network(s) raw edge
                                 weights with different criteria.
        AlloViz.Protein.analyze : Class method to visualize the network on the protein
                                    structure.

        Examples
        --------
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate("dynetan", cores=6, taskcpus=2)
        >>> opioidGPCR.filter("dynetan", "GetContacts_edges")
        >>> opioidGPCR.analyze("dynetan", "GetContacts_edges", "nodes", "btw")
        >>> opioidGPCR.view("dynetan", "btw", "GetContacts_edges", element=["edges", "nodes"])
        NGLWidget()
        """
        # Get the package name with the correct nomenclature or raise an error if it doesn't exist
        pkg = utils.pkgname(pkg)
        # Create a list of the passed elements
        elements = element if isinstance(element, list) else [element]

        # Define a function to retrieve the corresponding Element attribute of the passed elem and check that it exists
        get_element = lambda element: rgetattr(
            self, pkg, filtering.capitalize(), element.lower()
        )
        if isinstance(get_element(elements[0]), bool):
            raise Exception(
                f"{elements[0]} analysis with {filtering} filtering needs to be sent first."
            )

        # Retrieve the NGLWidget from the Visualization module's class method `view`
        nv = get_element(elements[0]).view(metric, num, colors, nv)

        # If both elements were passed, check that the second Element attribute exists and update the NGLWidget
        if len(element) == 2:
            if not get_element(elements[1]):
                raise Exception(
                    f"{elements[1]} analysis with {filtering} filtering needs to be sent first."
                )
            # The NGLWidget object `nv` is an optional argument of the `view` method and it is used to add the shapes to the passed NGLWidget
            nv = get_element(elements[1]).view(metric, num, colors, nv)

        # Returning of the NGLWidget is needed for its visualization in the notebook
        return nv


class Delta:
    r"""AlloViz class for calculating a delta-network.

    Used to calculate the delta-network between two :class:`~AlloViz.Protein` objects,
    using all the available combinations of packages/network construction methods,
    filterings, and elements-metrics that they have in common. A structural alignment of
    the structures is performed with PyMOL with the selected method to find the
    corresponding residues between the two structures for subtraction of edge weights.

    Parameters
    ----------
    refstate : :class:`AlloViz.Protein`
        Object of the :class:`AlloViz.Protein` class to use as reference (values of the
        other :class:`AlloViz.Protein` object will be subtracted from this one's values).
    state2 : :class:`AlloViz.Protein`
        Object of the :class:`AlloViz.Protein` class to compare with the reference one to
        build the delta-network (this one's values will be subtracted from the reference
        one's values).
    pymol_aln : {"super", "align", "cealign"}, default: "super"
        PyMOL's method to apply for the structural alignment of the two structures, used
        to retrieve the corresponding residues between the two for subtraction of edge
        weights. It is performed by :meth:`AlloViz.Delta._make_struct_aln`.

    Attributes
    ----------
    refstate
    state2

    See Also
    --------
    AlloViz.Protein : AlloViz main class for calculation of protein allosteric
                      communication networks.

    Notes
    -----
    For PyMOL's structural alignment, `align` starts from a sequence alignment and is
    optimal to align structures with identical sequences, but is deeply affected by
    sequence differences. `super` performs a sequence-independent dynamic programming
    structural alignment and it works well for structures with low sequence identity.
    `cealign` uses the Combinatorial Extension (CE) algorithm and it is preferred for
    structures with little to no sequence similarity (twilight zone). The methods are
    described in `PyMOL's documentation <https://pymolwiki.org/index.php/Main_Page>`_.

    Examples
    --------
    >>> activeB2AR = AlloViz.Protein(GPCR=117, name="Active B2AR")
    >>> inactiveB2AR = AlloViz.Protein(GPCR=160, name="Inactive B2AR")
    >>> for protein in [activeB2AR, inactiveB2AR]:
    >>>     protein.calculate("pytraj_CA")
    >>>     protein.analyze(metrics="btw")
    >>> delta = AlloViz.Delta(activeB2AR, inactiveB2AR)
    >>> print(delta.pytraj_CA.Whole.edges.df.shape)
    (40690, 5)
    """

    def __init__(self, refstate, state2, pymol_aln="super"):
        # Store the passed Proteins as attributes and as a private list
        self.refstate, self.state2 = refstate, state2
        self._states = [self.refstate, self.state2]
        # Add the present object as a new private attribute of each Protein object
        for state in self._states:
            setattr(state, "_delta", self)

        # Perform a PyMOL structural alignment with the passed method
        self._aln = self._make_struct_aln(pymol_aln)

        # Future: used source-target nodes for metrics calculation
        # nodes_subset = self._add_nodes_subset()
        # self.sources_subset, self.targets_subset = nodes_subset.values()

        # Calculate delta-network, simply by exploiting the classes' explicit __sub__ special methods
        delta = self.refstate - self.state2
        # Add the attributes of the calculated delta-network to the present object
        self.__dict__.update(delta)

    def _make_struct_aln(self, pymol_aln):
        r"""Perform a PyMOL structural alignment.

        Return a :class:`Bio.Align.MultipleSeqAlignment` as a result of parsing the
        ClustalW-formatted structural alignment provided by PyMOL. The alignment is also
        used to create an `_aln_mapper` attribute in each Protein of the Delta object
        that maps the alignment positions to its residues, which can/will be used by the
        Protein's private method `_translate_ix`.

        Parameters
        ----------
        pymol_aln : {"super", "align", "cealign"}
            PyMOL's method to apply for the structural alignment of the two structures,
            usedto retrieve the corresponding residues between the two for subtraction of
            edge weights.

        Returns
        -------
        :class:`Bio.Align.MultipleSeqAlignment`

        Notes
        -----
        Check `PyMol's documentation <https://pymolwiki.org/index.php/Main_Page>`_ for
        descriptions of the structural alignment methods and their fitness for identical
        or dissimilar structures.
        """
        import pymol2
        from Bio import AlignIO
        from Bio.SeqUtils import seq1

        # Perform PyMOL structural alignment: https://bioinformatics.stackexchange.com/questions/19105/perform-protein-structure-based-sequence-alignment-in-python
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(self.refstate._pdbf, "prot1")
            pymol.cmd.load(self.state2._pdbf, "prot2")
            eval(f"pymol.cmd.{pymol_aln}('prot1', 'prot2', object='aln')")

            with pymol.cmd.lockcm:
                aln = pymol2._cmd.get_seq_align_str(pymol.cmd._COb, "aln", -1, 0, -1)

        # Parse retrieved alignment string with Biopython's AlignIO
        with io.StringIO(aln) as f:
            alignment = AlignIO.read(f, "clustal")

        # Create the _aln_mapper attribute in each Protein passed to the Delta object to map the alignment's positions to its residues
        iterator = zip(
            self._states, alignment
        )  # looks like: ((refstate, alignment[0]), (state2, alignment[1]))
        for state, aln in iterator:
            # Get a list of residues in the Protein in the standardized name used in AlloViz
            aas = [f"{res.resname}:{res.resid}" for res in state.protein.residues]

            # Make a dictionary/mapper that maps each alignment position to a residue of the Protein
            ## The 'aln' object from the alignment has a `seq` attribute, which is a string with the same length as the alignment of the two sequences
            ## In each position of 'aln', there is either a one-letter code of a residue of the sequence or a "-" indicating a gap with respect to the other sequence
            ## So, each position of the alignment that is not a "-" corresponds to a residue of the Protein that will be in the 'aas' list
            ## 'aln_mapper' is generated by mapping the "popped first element of aas" to the corresponding alignment position (a pop and a position is skipped when it is "-")
            aln_mapper = dict(
                (aas.pop(0), pos)
                for pos, res in enumerate(aln.seq)
                if res != "-" and res == seq1(aas[0].split(":")[0])
            )
            # Also store the reverse mapping
            aln_mapper.update(dict(map(reversed, aln_mapper.items())))
            # And set the attribute in the Protein object
            setattr(state, "_aln_mapper", aln_mapper)

        return alignment

    #     def _add_nodes_subset(self):
    #         bfac = lambda atom: f"{atom.tempfactor:.2f}"
    #         is_wb = lambda atom: atom.name == "N" and -8.1 <= atom.tempfactor <= 8.1 and (bfac(atom) != "1.00" and bfac(atom) != "0.00")

    #         targets = [3.50, 6.30, 7.49, 7.50, 7.51, 7.52, 7.53]
    #         # In Ballesteros-Weinstein: Ionic lock 3.50 and 6.30; NPxxY 7.49-7.53

    #         def get_sources(states):
    #             resnum_to_wb = lambda nums, state: (bfac(atom)
    #                                                 for residue in nums
    #                                                 for atom in state.mdau.select_atoms(f"resnum {residue}")
    #                                                 if is_wb(atom))

    #             has_lig = [state if ("LIG" in (seg.segid for seg in state.mdau.segments)) else False for state in self.states]
    #             if any(has_lig):
    #                 state = next(item for item in has_lig if item != False)
    #                 aas = state.mdau.select_atoms("(same residue as around 4 segid LIG) and protein").residues
    #                 return list(resnum_to_wb(aas.resnums, state))
    #             else:
    #                 return [3.32] # Conserved position

    #         nodes_subset = {"sources": [val for val in get_sources(self.states) if val not in targets],
    #                         "targets": list(map(lambda x: f"{x:.2f}", targets))}

    #         format_res = lambda aa: f"A:{aa.resname}:{aa.resnum}"
    #         wb_to_aa = lambda wblist, state: list(map(format_res, (atom.residue
    #                                                                for residue in state.mdau.select_atoms("protein").residues
    #                                                                for atom in residue.atoms
    #                                                                if is_wb(atom) and bfac(atom) in wblist)))

    #         for state in self.states:
    #             state.sources_subset, state.targets_subset = wb_to_aa(nodes_subset["sources"], state), wb_to_aa(nodes_subset["targets"], state)

    #         return nodes_subset

    def view(
        self,
        pkg,
        metric,
        filtering="Whole",
        element="edges",
        num=20,
        colors=["orange", "turquoise"],
        nv=None,
    ):
        r"""Visualize the analyzed delta-networks on the protein structure.

        Visualize the delta-network corresponding to the selected package/network
        construction method, filtering scheme, network element and network metric or
        metrics, all of which must be already present in the delta-network object by
        having been previously calculated and analyzed in the respective Protein objects,
        and then correctly processed during the Delta object creation. Structure of the
        reference protein is used for displaying by default. The number of (each of the)
        elements to show and the color used for each Protein can be specified (first
        color will be used for negative values and thus for the second Protein, while the
        second color will be used for positive values and thus the reference Protein).

        Parameters
        ----------
        pkg : str, default: "all"
            Package/Network construction method for which to show the delta-network.
        metric : str, default: "all"
            Network metric for which to show the delta-network.
        filtering : str, {"Whole", "Incontact", "Intercontact"}
            Filtering scheme for which to show the delta-network.
        element : str or list, {"edges", "nodes"}
            Delta-network element or elements to show on the protein structure
            representing the chosen package, filtering scheme and metric.
        num : int, default: 20
            Number of (each of the) delta-network elements to show on the structure.
        colors : list, default: ["orange", "turquoise"]
            List of two colors to assign to the minimum and maximum values of the
            delta-network to be represented, respectively. Abstractly, the first color
            will be used for negative values and thus for elements that have higher value
            in the second Protein, while the second color will be used for positive
            values, and thus for the reference Protein. Dela-network values are
            interpolated setting 0 as the middle value, which is assigned "white".
        nv : :class:`nglview.NGLWidget`, optional
            A structure representation into which the shapes representing the chosen
            delta-network elements will be added.

        Returns
        -------
        :class:`nglview.NGLWidget`

        Examples
        --------
        >>> activeB2AR = AlloViz.Protein(GPCR=117, name="Active B2AR")
        >>> inactiveB2AR = AlloViz.Protein(GPCR=160, name="Inactive B2AR")
        >>> for protein in [activeB2AR, inactiveB2AR]:
        >>>     protein.calculate("pytraj_CA")
        >>>     protein.analyze(metrics="btw")
        >>> delta = AlloViz.Delta(activeB2AR, inactiveB2AR)
        >>> delta.view("pytraj_CA", "btw")
        NGLWidget()
        """
        # Function is the same one as the Protein class one but 'self' is passed to use Delta's attributes' data
        return Protein.view(self, pkg, metric, filterby, element, num, colors, nv)
