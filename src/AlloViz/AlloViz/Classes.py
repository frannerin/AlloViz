import os, io, re

from multiprocess import Pool


from .Analysis import Whole, Incontact, Intercontact
from .Visualization import Edges, Nodes

from .utils import rgetattr, rhasattr, capitalize
from . import utils
from . import trajutils

from .. import Wrappers






class _Store:
    pass

class Protein:
    r"""AlloViz main class.

    Objects of the class store the initialization information about the protein
    structure and trajectory files, process it and allow to calculate the different
    networks and network analyses with associated methods.
    
    Last updated on 3/6

    Parameters
    ----------
    pdb : str
        Filename of the PDB structure to read, process and use.
    trajs : str or list
        Filename(s) of the MD trajectory (or trajectories) to read and use. File format
        must be recognized by MDAnalysis (e.g., xtc), and other network construction
        packages may have stricter file format requirements.
    GPCR : bool or int, default: False
        Use `True` if the structure is a GPCR, or pass the ID of a GPCRmd database
        dynamics entry, without specifying the `pdb` nor the `trajs` parameters, to 
        automatically retrieve the files from the database and process them. They will
        be downloaded to `path`, which if left undefined will default to the GPCRmd ID.
    name : str, default: `pdb`
        Name string of the protein/class instance to be used for visualization, e.g., 
        for the title of the colorbar shown when representing a network with nglviewer.
    path : str, optional
        Path to store results in. It can exist already or not, and a new folder inside
        it called `data` will be also created to store computation results and other
        files. If unspecified, it defaults to "." or the GPCRmd ID in case it is used.

    Other Parameters
    ----------------
    protein_sel : str, default: "same segid as protein"
        MDAnalysis atom selection string to select the protein structure from the 
        Universe (e.g., in case simulation in biological conditions are used, with
        water molecules, ions...). Some network construction methods use the source
        files directly, and might end up providing data with a different size than that
        obtained with this selection, which will cause errors.
    psf : str, optional
        Filename of the .psf file corresponding to the pdb used. It is required by gRINN.
    parameters : str, optional
        Filename of the MD simulation force-field parameters file in NAMD format. It is
        required by gRINN.
    **kwargs
        A dictionary can be passed with the `special_res` keyword argument, containing
        a mapping of special residue 3/4-letter code(s) present in the structure to the
        corresponding 1-letter code(s). Most times it is not necessary, and it may also
        cause problems with network construction methods that do not interpret these 
        special residues as part of the protein, which will return data with size
        different to what is expected causing errors.
    
    Attributes
    ----------
    pdb
    trajs
    path
    GPCR
    pdbu : MDAnalysis.core.Universe
        Universe of the pdb file provided.
    prot : MDAnalysis.core.groups.AtomGroup
        AtomGroup of the selected `protein_sel` string, taken from the `pdbu` Universe.
    mdau : MDAnalysis.core.Universe
        Universe of the pdb and trajectory file(s) provided.
    
    Raises
    ------
    FileNotFoundError
        If any of the files passed in `pdb` and `traj` (and `psf` and `parameters` if 
        provided) parameters cannot be accessed.
        
    See Also
    --------
    AlloViz.Delta : Class for calculation of the delta-network(s) between two Protein
                    objects.

    Examples
    --------
    >>> opioidGPCR = AlloViz.Protein(GPCR=169)
    >>> print(opioidGPCR.mdau)
    <Universe with 88651 atoms>
    """
#     Notes
#     -----
#     Notes about the implementation algorithm (if needed).

#     This can have multiple paragraphs.

#     You may include some math:

#     .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

#     And even use a Greek symbol like :math:`\omega` inline.

#     References
#     ----------
#     Cite the relevant literature, e.g. [1]_.  You may also cite these
#     references in the notes section above.

#     .. [1] O. McNoleg, "The integration of GIS, remote sensing,
#        expert systems and adaptive co-kriging for environmental habitat
#        modelling of the Highland Haggis using object-oriented, fuzzy-logic
#        and neural-network techniques," Computers & Geosciences, vol. 22,
#        pp. 585-588, 1996.
    
    def __init__(self, pdb="", trajs=[], GPCR=False, name=None, path="", protein_sel="same segid as protein", psf=None, parameters=None, **kwargs):
        self.GPCR = GPCR
        
        if not isinstance(self.GPCR, bool) and isinstance(self.GPCR, int):
            self._gpcrmdid = self.GPCR
            self._path = f"{self.GPCR}" if len(path) == 0 else path
            os.makedirs(self._path, exist_ok=True)
            if not any([re.search("(pdb$|psf$|xtc$|parameters$)", file) for file in os.listdir(self._path)]):
                self._download_files()
                
            files = os.listdir(self._path)
            self._pdbf = self._path + '/' + next(file for file in files if re.search("pdb$", file))
            self._psff = self._path + '/' + next(file for file in files if re.search("psf$", file))
            self._paramf = self._path + '/' + next(file for file in files if re.search("parameters$", file))
            self._trajs = dict(enumerate(sorted( f"{self._path}/{traj}" for traj in files if re.search("^(?!\.).*\.xtc$", traj) ), 1))
        
        else:
            psf_params = parameters is not None and psf is not None
            
            self._pdbf = pdb
            self._trajs = {1: trajs} if isinstance(trajs, str) else dict(enumerate(trajs, 1))
            
            files_to_check = list(self._trajs.values()) + [self._pdbf] if not psf_params else list(self._trajs.values()) + [self._pdbf, psf, parameters]
            files_exist = {file: os.path.isfile(file) for file in files_to_check}
            if any([not file_exist for file_exist in files_exist.values()]):
                raise FileNotFoundError(f"Some of the files could not be found: {files_exist}")
            
            self._path = "." if len(path) == 0 else path
            os.makedirs(self._path, exist_ok=True)
            if psf_params:
                self._psff = psf
                self._paramf = parameters
        
        # Temporary patch
        self.pdb = self._pdbf
        self.trajs = list(self._trajs.values())
        self.path = self._path
        
        self.name = name if name is not None else self._pdbf
        self._protein_sel = protein_sel
                
        self._datadir = f"{self._path}/data"
        os.makedirs(self._datadir, exist_ok=True)
        
        self._res_dict = self._get_res_dict(**kwargs)
        self.pdbu, self.prot, self.u = self._get_mdau(self._protein_sel)
        self._dihedral_residx = self._get_dihedral_residx(self.prot)
        self._translate_ix = lambda mapper: lambda ix: tuple(mapper[_] for _ in ix) if isinstance(ix, tuple) else mapper[ix]
        
        
        
    _download_files = trajutils.download_files
    _get_mdau = trajutils.get_mdau
    _add_comtrajs = trajutils.add_comtrajs
    _make_dcds = trajutils.make_dcds
    _get_bonded_cys = trajutils.get_bonded_cys
    _get_res_dict = trajutils.get_res_dict
    _get_dihedral_residx = staticmethod(trajutils.get_dihedral_residx)
    
    
    
    
    def __sub__(self, other):
        delta = _Store()
        
        for pkg in (key for key in self.__dict__ if key.lower() in [x.lower() for x in utils.pkgsl] and key in other.__dict__):
            setattr(delta, pkg, _Store())
            for filterby in (key for key in getattr(self, pkg).__dict__ if key.lower() in ["whole", "incontact", "intercontact", "raw"] and key in getattr(other, pkg).__dict__): #if not re.search("(^_|raw)", key)
                setattr(getattr(delta, pkg), filterby, _Store())
                for elem in (key for key in rgetattr(self, pkg, filterby).__dict__ if key.lower() in ["nodes", "edges"] and key in rgetattr(other, pkg, filterby).__dict__): 
                # setattr(getattr(delta, pkg), filterby, _Store())
                # for norm in (key for key in rgetattr(self, pkg, filterby).__dict__ if key in ["norm", "no_norm"] and key in rgetattr(self, pkg, filterby).__dict__):
                    # dif = rgetattr(self, pkg, filterby, norm) - rgetattr(other, pkg, filterby, norm)
                    # setattr(rgetattr(delta, pkg, filterby), norm, Edges(self._pair, dif))
                    dif = rgetattr(self, pkg, filterby, elem) - rgetattr(other, pkg, filterby, elem)
                    elemclass = eval(elem.capitalize())
                    setattr(rgetattr(delta, pkg, filterby), elem, elemclass(self._delta, dif))
                    
        return delta.__dict__
        
        
        
    
    def calculate(self, pkg="all", cores=1, **kwargs): # ow
        r"""Calculate edge weights of allosteric networks.

        Send the computation of the raw edge weights for the selected network
        construction methods.

        Parameters
        ----------
        pkg : str or list, default: "all"
            Package(s)/Network construction method(s) for which to send raw edge weight
            computation. "all" sends the computation for all available methods within
            AlloViz (check `AlloViz.AlloViz.utils.pkgsl`).
        cores : int, default: 1
            Number of cores to use for parallelization with a `multiprocess` Pool.
            Default value only uses 1 core with a custom `dummypool` (check
            `AlloViz.utils.dummypool`) that perform computations synchronously.

        Returns
        -------
        None

        Other Parameters
        ----------------
        **kwargs
            `taskcpus` keyword argument can be passed to specify the amount of cores
            that parallelizable network construction methods can use (i.e., AlloViz's
            method, getcontacts, dynetan, PyInteraph, MDEntropy and gRINN). `namd`
            keyword argument can be passed to point to the namd2 executable location; if
            the `namd` command is accessible through the CLI it is automatically
            retrieved with the `distutils` package.

        See Also
        --------
        AlloViz.Protein.analyze : Class method to analyze the calculated raw edge weights
                                  with graph theory-based methods.
        AlloViz.Protein.visualize : Class method to visualize the network on the protein
                                    structure.
                                    
        Notes
        -----
        Method returns nothing, but calculation results are stored as new instance
        attributes with the same name as the selected packages/network construction
        methods. If the object is created providing more than one trajectory file,
        the average and standard error of the weights between the replicas are also
        calculated.

        Examples
        --------        
        If we have a 6-core computer and want to use 2 cores to compute the values for
        each of the three trajectories of the protein (GPCRmd stores three replicas of
        each structure):
        
        >>> opioidGPCR = AlloViz.Protein(GPCR=169)
        >>> opioidGPCR.calculate("dynetan", cores=6, taskcpus=2)
        >>> print(opioidGPCR.dynetan.raw.shape)
        (41041, 5)
        """
        pkgs = utils.pkgsl if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        
        if any(["COM" in pkg for pkg in pkgs]):
            self._add_comtrajs(self._protein_sel)
            
        if any([re.search("(carma|grinn)", pkg.lower()) for pkg in pkgs]):
            self._make_dcds()
        
        d = self.__dict__.copy()
        d.update(kwargs)
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        # taskcpus = kwargs.pop("taskcpus") if "taskcpus" in kwargs else cores
        
        for pkg in pkgs: #self._set_pkgclass(self, pkg, d) #**kwargs)
            try:
                pkg_ix = [pkg.lower() for pkg in utils.pkgsl].index(pkg.lower())
            except:
                raise Exception(f"{pkg} isn't a valid name of an AlloViz network construction method.")
    
            pkgclass = eval(f"Wrappers.{utils.pkgsl[pkg_ix]}")
            
            if not hasattr(self, pkgclass.__name__):
                setattr(self, pkgclass.__name__, pkgclass(self, d))
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
    
    
    
    def analyze(self, pkg="all", filterby="Whole", element="edges", metrics="all", normalize=True, cores=1,
                nodes_dict = {'btw': 'networkx.algorithms.centrality.betweenness_centrality',
                              'cfb': 'networkx.algorithms.centrality.current_flow_betweenness_centrality'},
                edges_dict = {'btw': 'networkx.algorithms.centrality.edge_betweenness_centrality',
                              'cfb': 'networkx.algorithms.centrality.edge_current_flow_betweenness_centrality'}):
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
            Default value only uses 1 core with a custom `dummypool` (check
            `AlloViz.utils.dummypool`) that perform computations synchronously.
        nodes_dict, edges_dict : dict
            Dictionary that maps network metrics custom names (e.g., betweenness
            centrality, "btw") with their corresponding NetworkX function (e.g.,
            "networkx.algorithms.centrality.betweenness_centrality"). Functions strings
            must be written as if they were absolute imports, and must return a
            dictionary of edges or nodes, depending on the element dictionary in which
            they are. The keys of the dictionaries will be used to name the columns of
            the analyzed data that the functions produce.

        See Also
        --------
        AlloViz.Protein.calculate : Class method to calculate the network(s) raw edge
                                    weights with different network construction methods.
        AlloViz.Protein.visualize : Class method to visualize the network on the protein
                                    structure.

        Notes
        -----
        Method returns nothing, but analysis results are stored as nested attributes
        "inside" each of the packages' attributes of the `Protein` object, first using
        the name of the filtering scheme (e.g., `.Package.filterby`) and lastly with the
        analyzed network element name (e.g., `.Package.filterby.element`). If the
        `Protein` object was created providing more than one trajectory file, the
        analyses are performed both on the replicas' weights and the average, and an
        average and standard error of the replicas' analysis results are also calculated.
        
        "Incontact" filtering only retains edges of residue pairs in contact -those for
        which getcontacts is able to compute contact frequencies-, and "Intercontact"
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
        pkgs = [pkg for pkg in self.__dict__ if pkg in utils.pkgsl] if pkg=="all" else pkg if isinstance(pkg, list) else [pkg]
        filterbys = utils.filterbysl if filterby=="all" else filterby if isinstance(filterby, list) else [filterby]
        elements = element if isinstance(element, list) else [element]
        metrics = set(list(nodes_dict.keys()) + list(edges_dict.keys())) if metrics=="all" else metrics if isinstance(metrics, list) else [metrics]
        
        utils.pool = utils.dummypool()
        if cores>1:
            mypool = Pool(cores)
            utils.pool = mypool
        print(utils.pool)
        
        for pkgn in pkgs:
            # pkgclass = eval(f"Wrappers.{pkg}")#eval(f"Pkgs.{capitalize(pkg)}") if isinstance(pkg, str) else pkg
            pkg = rgetattr(self, pkgn)
            if not pkg:
                print(f"{pkgn} calculation results are needed first")
                continue
                
            for filterby in filterbys:
                anaclass = eval(f"{capitalize(filterby)}") if isinstance(filterby, str) else filterby
                filterby = lambda: rgetattr(pkg, anaclass.__name__)
                if not filterby():
                     setattr(pkg, anaclass.__name__, anaclass(pkg))
                filterby()._add_metrics(elements, metrics, normalize, nodes_dict, edges_dict)
        
        if cores>1:
            mypool.close()
            mypool.join()
            utils.pool = utils.dummypool()
    
        
    
    
        
    def view(self, pkg, metric, filterby="incontact", element:list=["edges"], num=20, colors=["orange", "turquoise"]):        
        get_element = lambda element: rgetattr(self, capitalize(pkg), capitalize(filterby), element.lower())
        
        nv = get_element(element[0]).view(metric, num, colors)
        
        if len(element) == 2:
            nv = get_element(element[1]).view(metric, num, colors, nv)
            
        return nv


    
    
    
class Delta:
    def __init__(self, refstate, state2, **kwargs):
        self.state1, self.state2 = refstate, state2
        self._states = [self.state1, self.state2]
        for state in self._states:
            setattr(state, "_delta", self)
            
            # translate_ix = lambda ix: tuple(to_aln_pos[_] for _ in ix) if isinstance(ix, tuple) else to_aln_pos[ix]
            
            
        self._aln = self._make_struct_aln(**kwargs)
        
        # nodes_subset = self._add_nodes_subset()
        # self.sources_subset, self.targets_subset = nodes_subset.values()
        
        
        delta = self.state1 - self.state2
        self.__dict__.update(delta)
    
    
    
    def _make_struct_aln(self, aln_method="TMalign_pair"):
        # ****** Pairwise Structural Alignment Methods:
        # --------------------------------------------
        # align_pdbpair        built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # lalign_pdbpair       built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # extern_pdbpair       built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # thread_pair          built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # fugue_pair           http://mizuguchilab.org/fugue/                             [pg:        fugueali is NOT Installed][]
        # pdb_pair             built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # sap_pair             https://mathbio.crick.ac.uk/wiki/Software#SAP              [pg:             sap is  Installed][/gpcr/users/frann/networks/tcoffee/bin/sap]
        # sara_pair            built_in                                                   [pg:        t_coffee is  Installed][built_in]
        # daliweb_pair         built_in                                                   [pg:     dalilite.pl is  Installed][built_in]
        # dali_pair            built_in                                                   [pg:     dalilite.pl is  Installed][built_in]
        # mustang_pair         http://lcb.infotech.monash.edu.au/mustang/                 [pg:         mustang is  Installed][/gpcr/users/frann/networks/tcoffee/bin/mustang]
        # TMalign_pair         http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f      [pg:         TMalign is  Installed][/gpcr/users/frann/networks/tcoffee/bin/TMalign]
        
        from subprocess import Popen, PIPE
        from Bio import AlignIO
        from Bio.SeqUtils import seq1
        from distutils.spawn import find_executable

        path = find_executable('t_coffee').replace('/t_coffee', '')
        env = {"HOME": path, "PATH": os.environ["PATH"] + f":{path}"}
        tcoffee = Popen(f"t_coffee -pdb={self.state1._protpdb},{self.state2._protpdb} -method {aln_method} -outfile no -template_file no -output clustalw -align".split(" "),
                        stdout=PIPE, stderr=PIPE, encoding="utf-8", env=env)
        # tmalign = Popen(f"TMalign {self.state1._protpdb} {self.state2._protpdb}".split(" "),# -method {aln_method} -outfile no -template_file no -output clustalw -align".split(" "),
        #                 stdout=PIPE, stderr=PIPE, encoding="utf-8")#, env=env)
        
        stderr = tcoffee.stderr.read()
        # stderr = tmalign.stderr.read()
        if "error" not in stderr.lower():
            alignment = AlignIO.read(io.StringIO(tcoffee.stdout.read()), "clustal")
#             aln = tmalign.stdout.readlines()[-4:]
#             fasta = f""">{self.state1._protpdb}
# {aln[0].strip()}
# >{self.state2._protpdb}
# {aln[2].strip()}
# """
#             alignment = AlignIO.read(io.StringIO(fasta), "fasta")
        else:
            print(stderr)
            raise Exception("Structural alignment of the structures couldn't be performed.")
        
        for state, aln in zip(self._states, alignment):
            aas = [f"{res.resname}:{res.resid}" for res in state.prot.residues]
            aln_mapper = dict( (aas.pop(0), pos) for pos, res in enumerate(aln.seq) if res != "-" and res == seq1(aas[0].split(":")[0], custom_map = state._res_dict) )
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
    
    
    
    def view(self, pkg, metric, filterby="incontact", element:list=["edges"], num=20, colors=["orange", "turquoise"]):
        # norm = self._get_norm_str(normalize)
        # if not rhasattr(self.data, filterby, norm, pkg):
        #     self.analyze(pkg, normalize=normalize, filterby=filterby)
        
        get_element = lambda element: rgetattr(self, capitalize(pkg), capitalize(filterby), element.lower())
        
        nv = get_element(element[0]).view(metric, num, colors)
        
        if len(element) == 2:
            nv = get_element(element[1]).view(metric, num, colors, nv)
            
        return nv