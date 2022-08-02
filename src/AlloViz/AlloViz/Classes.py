# \
# =========================================================
# Core object: Universe --- :mod:`MDAnalysis.core.universe`
# =========================================================
"""Main AlloViz classes: `Protein` and `Delta` --- :mod:`AlloViz`

The :class:`AlloViz.Protein` class is AlloViz's main class, processing a structure and
trajectory(ies) input files and allowing to calculate, analyze and visualize the
allosteric communication networks with its associated methods.

The :class:`AlloViz.Delta` class takes two :class:`AlloViz.Protein` objects as input to
calculate the delta-network to highlight the differences in the allosteric communication
between two systems.

Classes
=======
.. autoclass:: Protein
   :members:
.. autoclass:: Delta
   :members:
   
Notes
=====
This module is not meant to be used directly, as the classes are imported with the
namespace of AlloViz itself.

"""

import os, io, re

from multiprocess import Pool
import MDAnalysis as mda
import numpy as np

from .Analysis import Whole, Incontact, Intercontact
from .Visualization import Edges, Nodes

from .trajutils import ProteinBase
from .utils import rgetattr, rhasattr
from . import utils
#from . import trajutils

from .. import Wrappers






class _Store:
    pass

class Protein(ProteinBase):
    r"""AlloViz main class.

    Objects of the class store the initialization information about the protein
    structure and trajectory files, process it and allow to calculate the different
    networks and network analyses with the associated methods.

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
        e.g., :attr:`AlloViz.Protein._protein_sel` + " and {customselection}".

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
    protein : :class:`MDAnalysis.core.groups.AtomGroup`
        AtomGroup of the selected `protein_sel` string, taken from the pdb file.
    u : :class:`MDAnalysis.core.Universe`
        Universe of the pdb and trajectory files with only the `protein_sel` atoms.
    
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
    
    def __init__(self, pdb="", trajs=[], GPCR=False, name=None, path=None, protein_sel=None, **kwargs):
        self.GPCR = GPCR
        
        # If a GPCRmd ID is passed
        if not isinstance(self.GPCR, bool) and isinstance(self.GPCR, int):
            self.name = f"{self.GPCR}" if not name else name
            self._path = f"{self.GPCR}" if not path else path
            os.makedirs(self._path, exist_ok=True)
            
            # Download the files from GPCRmd
            if not any([re.search("(pdb$|psf$|xtc$|parameters$)", file) for file in os.listdir(self._path)]):
                self._download_GPCRmd_files()
            files = os.listdir(self._path)
            
            # Retrieve filenames from the files downloaded into self._path
            get_filename = lambda ext: self._path + '/' + next(file for file in files if re.search(f"{ext}$", file))
            self.pdb = get_filename("pdb")
            self.trajs = list(sorted( f"{self._path}/{traj}" for traj in files if re.search("^(?!\.).*\.xtc$", traj) ))
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
            files_to_check = self.trajs + [self.pdb] if not passed_psf_params else self.trajs + [self.pdb, self.psf, self._paramf]
            files_exist = {file: os.path.isfile(file) for file in files_to_check}
            if any([not file_exist for file_exist in files_exist.values()]):
                raise FileNotFoundError(f"Some of the files could not be found: {files_exist}")
            

        self._protein_sel = Protein._protein_sel if not protein_sel else protein_sel
        
        # Set the names of the directories and files that will be creating when processing the input files
        self._datadir = f"{self._path}/data"
        os.makedirs(self._datadir, exist_ok=True)
        
        self._pdbf = f"{self._datadir}/protein.pdb"
        self._trajs = dict( [(num+1, f"{self._datadir}/traj_{num+1}.xtc") for num in range(len(self.trajs))] )
        self._psff = self._pdbf.replace("pdb", "psf") if hasattr(self, "psf") else None
        
        # Names of the directories and files of the future pdb and trajectory(ies) of the residues' Center Of Mass
        compath = f"{self._datadir}/COM_trajs"
        os.makedirs(compath, exist_ok=True)
        self._compdbf = f"{compath}/ca.pdb"
        self._comtrajs = {num: f"{compath}/{num}.xtc" for num in self._trajs}
        
        # If the processed filenames don't exist yet as files, process the input files; if special_res kwarg is used it will be passed
        if any([not os.path.isfile(f) for f in list(self._trajs.values()) + list(self._comtrajs.values()) + [self._pdbf]]):
            self._process_input(**kwargs)
        
        # Set the protein/pdb and trajectory(ies) MDAnalysis' Universes with the processed files
        self.protein = mda.Universe(self._pdbf)
        self.u = mda.Universe(self._pdbf, *list(self._trajs.values()))
        

        # Bonded cysteines must be identified to remove them from non-covalent contacts calculations (i.e., PyInteraph2_Contacts)
        self._bonded_cys = self._get_bonded_cys()
        
        # _dihedral_residx returns a list of indices of all the protein's residues in the Universes that are not in the extremes of the chain(s)
        # These will be the residues for which Dihedral correlations can be calculated (residues in the extremes have some dihedrals missing)
        # MDEntropy_AlphaAngle uses information of residues i-1, i, i+1, i+2; so end=-2 is passed to also exclude the second-to-last residue of the chain(s)
        _res_arrays = np.split(self.protein.residues.resindices, np.where(np.diff(self.protein.residues.resnums) != 1)[0]+1)
        self._dihedral_residx = lambda end=-1: [elem for arr in _res_arrays for elem in arr[1:end]]
        
        # _translate_ix returns a translated DataFrame Index (string) or MultiIndex (tuple) element according to the mapper dictionary passed
        # The translation function _translate_ix(mapper) is passed as the function used by the .map method of a DataFrame Index or MultiIndex to translate it
        self._translate_ix = lambda mapper: lambda ix: tuple(mapper[_] for _ in ix) if isinstance(ix, tuple) else mapper[ix]
    


    
    def __sub__(self, other):
        delta = _Store()
        
        for pkg in (key for key in self.__dict__ if key.lower() in [x.lower() for x in utils.pkgsl] and key in other.__dict__):
            setattr(delta, pkg, _Store())
            for filterby in (key for key in getattr(self, pkg).__dict__ if key.lower() in utils.filterbysl and key in getattr(other, pkg).__dict__):
                setattr(getattr(delta, pkg), filterby, _Store())
                for elem in (key for key in rgetattr(self, pkg, filterby).__dict__ if key.lower() in ["nodes", "edges"] and key in rgetattr(other, pkg, filterby).__dict__):
                    dif = rgetattr(self, pkg, filterby, elem) - rgetattr(other, pkg, filterby, elem)
                    elemclass = eval(elem.capitalize())
                    setattr(rgetattr(delta, pkg, filterby), elem, elemclass(self._delta, dif))
                    
        return delta.__dict__
        
        
        
    
    def calculate(self, pkg="all", cores=1, **kwargs):
        r"""Calculate edge weights of allosteric networks.

        Send the computation of the raw edge weights for the selected network
        construction methods.

        Parameters
        ----------
        pkg : str or list, default: "all"
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
        namd : str, optional
            Optional kwarg pointing to the namd2 executable location; if the `namd`
            command is accessible through the CLI it is automatically retrieved with the
            `distutils` package.
        GetContacts_threshold : float, optional
            Optional kwarg to specify the minimum contact frequency (0-1, default 0)
            threshold for GetContacts results, which will be used to filter out contacts
            with a frequency (average) lower than it.

        See Also
        --------
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
        pkgs = utils.pkgsl if pkg=="all" else [utils.pkgname(p) for p in pkg] if isinstance(pkg, list) else [utils.pkgname(pkg)]
        
        # Objects from the classes in the Wrappers module need to be passed a dictionary "d" containing all the attributes of the source Protein object and the passed kwargs
        d = self.__dict__.copy()
        d.update(kwargs)
        
        # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
        # Changing it inside the utils module allows to share the same one between modules
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for pkg in pkgs:
            # Establish the corresponding Wrappers' class
            pkgclass = eval(f"Wrappers.{pkg}")
            
            # Setting the class as a new attribute will initialize all calculations asynchronously (synchronously if a dummypool is used)
            if not hasattr(self, pkgclass.__name__):
                setattr(self, pkgclass.__name__, pkgclass(self, d))
        
        # Close and reset the pool
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
    
    
    
    def analyze(self, pkg="all", filterby="Whole", element="edges", metrics="all", normalize=True, cores=1,
                nodes_dict = {'btw': 'networkx.algorithms.centrality.betweenness_centrality',
                              'cfb': 'networkx.algorithms.centrality.current_flow_betweenness_centrality'},
                edges_dict = {'btw': 'networkx.algorithms.centrality.edge_betweenness_centrality',
                              'cfb': 'networkx.algorithms.centrality.edge_current_flow_betweenness_centrality'},
                **kwargs):
        r"""Analyze calculated edge weights with network analyses.
        
        Perform analyses of the raw edge weights for the selected packages/network
        construction methods, filtering schemes and network elements, calculating the
        desired network metrics.

        Parameters
        ----------
        pkg : str or list, default: "all"
            Package(s)/Network construction method(s) for which to analyze their raw edge
            weights, which must be already calculated and their data saved as instance
            attribute. In this case, "all" sends the computation for all available
            methodsthat are already calculated and saved as instance attributes.
        filterby : str or list, {"Whole", "Incontact", "Intercontact"}
            Filtering schemes with which to filter the list of network edges before
            analysis. "Whole" is the default and means no filtering, and "Incontact" and
            "Intercontact" (see Notes) require that raw data from the `getcontacts`
            package has been calculated first.
        element : str or list, {"edges", "nodes"}
            Network elements for which to perform the analysis.
        metrics : str or list, default: "all"
            Network metrics to compute, which must be keys in the `nodes_dict` or
            `edges_dict` dictionaries. Default is "all" and it sends the computation for
            all the metrics defined in the corresponding dictionary of the selected
            elements in `element`.

        Returns
        -------
        None

        Other Parameters
        ----------------
        normalize : bool, default: True
            Passed to the NetworkX functions that calculate the metrics, to output
            normalized results or not.
        cores : int, default: 1
            Number of cores to use for parallelization with a `multiprocess` Pool.
            Default value only uses 1 core with a custom :class:`AlloViz.utils.dummypool`
            that performs computations synchronously.
        nodes_dict, edges_dict : dict
            Dictionary that maps network metrics custom names (e.g., betweenness
            centrality, "btw") with their corresponding NetworkX function (e.g.,
            "networkx.algorithms.centrality.betweenness_centrality"). Functions strings
            must be written as if they were absolute imports, and must return a
            dictionary of edges or nodes, depending on the element dictionary in which
            they are. The keys of the dictionaries will be used to name the columns of
            the analyzed data that the functions produce.
        **kwargs
            `GetContacts_threshold` kwarg can be passed to specify the minimum contact
            frequency (0-1, default 0) threshold, which will be used to filter out
            contacts with a frequency (average) lower than it before analysis. Make sure
            to delete/have deleted all previous analysis attributes and files of any 
            network construction method. `Intercontact_distance` kwarg can be passed to
            specify the minimum number of sequence positions/distance between residues of
            a pair to retain in Intercontact filtering, which defaults to 5.
            

        See Also
        --------
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.view : Class method to visualize the network on the protein
                               structure.

        Notes
        -----
        Method returns nothing, but analysis results are stored as nested attributes
        "inside" each of the packages' attributes of the :class:`AlloViz.Protein` object,
        first using the name of the filtering scheme (e.g., `.Package.filterby`) and
        lastly with the analyzed network element name (e.g., `.Package.filterby.element`). 
        If the `Protein` object was created providing more than one trajectory file, the
        analyses are performed both on the replicas' weights and the average, and an
        average and standard error of the replicas' analysis results are also calculated.
        
        "Incontact" filtering only retains edges of residue pairs in contact -those for
        which GetContacts is able to compute contact frequencies-, and "Intercontact"
        only keeps edges of pairs that are both in contact and apart in the sequence (more
        than 5 positions away in the sequence).
        
        Examples
        --------
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate(["getcontacts", "dynetan"], cores=6, taskcpus=2)
        >>> opioidGPCR.analyze("dynetan", filterby="Intercontact", element=["edges", "nodes"], metrics="btw")
        >>> print(opioidGPCR.dynetan.Intercontact.edges.df.shape)
        (3410, 5)
        """
        # Calculate for "all" packages (all the available packages that have been previously calculated and are in __dict__)
        # or the ones passed as parameter (check that they are on the list of available packages and retrieve their case-sensitive names, else raise an Exception)
        pkgs = [pkg for pkg in self.__dict__ if pkg in utils.pkgsl] if pkg=="all" else [utils.pkgname(p) for p in pkg] if isinstance(pkg, list) else [utils.pkgname(pkg)]
        # Calculate for the passed Filterings
        filterbys = utils.filterbysl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        # And for the passed Elements, and their corresponding network Metrics if they are defined in the dictionaries
        elements = element if isinstance(element, list) else [element]
        metrics = set(list(nodes_dict.keys()) + list(edges_dict.keys())) if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
        # Changing it inside the utils module allows to share the same one between modules
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        # For each package to be analyzed
        for pkgn in pkgs:
            pkg = rgetattr(self, pkgn)
            if not pkg:
                print(f"{pkgn} calculation results are needed first")
                continue
            
            # Add (or retrieve) the passed Filterings attributes, and add the passed Elements-Metrics to them
            for filterby in filterbys:
                # Establish the corresponding filterby/Analysis' class and check if it is already an attribute or set it otherwise
                anaclass = eval(filterby.capitalize()) if isinstance(filterby, str) else filterby
                get_ana = lambda: rgetattr(pkg, anaclass.__name__)
                if not get_ana():
                     setattr(pkg, anaclass.__name__, anaclass(pkg, **kwargs))
                
                # Add the passed Elements-Metrics
                get_ana().add_metrics(elements, metrics, normalize, nodes_dict, edges_dict)
        
        # Close and reset the pool
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
        
    
    
        
    def view(self, pkg, metric, filterby="Whole", element="edges", num:int=20, colors:list=["orange", "turquoise"]):
        r"""Visualize the analyzed networks on the protein structure.
        
        Visualize the network corresponding to the selected package/network construction
        method, filtering scheme, network element and network metric or metrics, all of
        which must be previously calculated and analyzed. The number of (each of the) 
        elements to show and the color scale can be specified.

        Parameters
        ----------
        pkg : str, default: "all"
            Package/Network construction method for which to show the network.
        metric : str, default: "all"
            Network metric for which to show the network.
        filterby : str, {"Whole", "Incontact", "Intercontact"}
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
        
        Returns
        -------
        NGLWidget()

        See Also
        --------
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.analyze : Class method to visualize the network on the protein
                                    structure.

        Notes
        -----
        ``view`` method is also available in the Elements' attributes, e.g.,
        `AlloViz.Protein.Package.Analysis.Element.view` 
        (see :mth:`AlloViz.Visualization.Element.view`).
        
        Examples
        --------
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate("dynetan", cores=6, taskcpus=2)
        >>> opioidGPCR.analyze("dynetan", element=["edges", "nodes"], metrics="btw")
        >>> opioidGPCR.view("dynetan", "btw", element=["edges", "nodes"])
        NGLWidget()
        """
        # Get the package name with the correct nomenclature or raise an error if it doesn't exist
        pkg = utils.pkgname(pkg)
        # Create a list of the passed elements
        elements = element if isinstance(element, list) else [element]
        
        # Define a function to retrieve the corresponding Element attribute of the passed elem and check that it exists
        get_element = lambda element: rgetattr(self, pkg, filterby.capitalize(), element.lower())
        if not get_element(elements[0]):
            raise Exception(f"{elements[0]} analysis with {filterby} filtering needs to be sent first.")
        
        # Retrieve the NGLWidget from the Visualization module's class method `view`
        nv = get_element(elements[0]).view(metric, num, colors)
        
        # If both elements were passed, check that the second Element attribute exists and update the NGLWidget 
        if len(element) == 2:
            if not get_element(elements[1]):
                raise Exception(f"{elements[1]} analysis with {filterby} filtering needs to be sent first.")
            # The NGLWidget object `nv` is an optional argument of the `view` method and it is used to add the shapes to the passed NGLWidget
            nv = get_element(elements[1]).view(metric, num, colors, nv)
        
        # Returning of the NGLWidget is needed for its visualization in the notebook
        return nv


    
    
    
class Delta:
    def __init__(self, refstate, state2, pymol_aln="super"):
        # https://www.researchgate.net/post/Which_command_should_I_use_to_do_stuctural_comparison_on_a_protein_binding_with_different_ligand           
        self.state1, self.state2 = refstate, state2
        self._states = [self.state1, self.state2]
        for state in self._states:
            setattr(state, "_delta", self)            
            
        self._aln = self._make_struct_aln(pymol_aln)
        
        # nodes_subset = self._add_nodes_subset()
        # self.sources_subset, self.targets_subset = nodes_subset.values()
        
        
        delta = self.state1 - self.state2
        self.__dict__.update(delta)
    
    
    def _make_struct_aln(self, pymol_aln):
        import pymol2
        from Bio import AlignIO
        from Bio.SeqUtils import seq1


        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(self.state1._pdbf, 'prot1')
            pymol.cmd.load(self.state2._pdbf, 'prot2')
            eval(f"pymol.cmd.{pymol_aln}('prot1', 'prot2', object='aln')")

            with pymol.cmd.lockcm:
                aln = pymol2._cmd.get_seq_align_str(pymol.cmd._COb, 'aln',
                        -1, 0, -1)
        
        with io.StringIO(aln) as f:
            alignment = AlignIO.read(f, "clustal")
        
        for state, aln in zip(self._states, alignment):
            aas = [f"{res.resname}:{res.resid}" for res in state.protein.residues]
            aln_mapper = dict( (aas.pop(0), pos) for pos, res in enumerate(aln.seq) if res != "-" and res == seq1(aas[0].split(":")[0]) )#, custom_map = state._res_dict
            aln_mapper.update( dict(map(reversed, aln_mapper.items())) )
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
    
    
    
    def view(self, pkg, metric, filterby="Whole", element:list=["edges"], num=20, colors=["orange", "turquoise"]):
        return Protein.view(self, pkg, metric, filterby, element, num, colors)