"""Module with classes to filter and store networks

Functions of the module (and/or combinations of them) are used for filtering the passed
raw networks with the main class :class:`~AlloViz.AlloViz.Filtering.Filtering`.

"""

import os
import time

import pandas
import numpy as np
import networkx
from multiprocess import Pool

from . import Analysis
from . import utils
from .utils import rgetattr, rhasattr


def All(pkg, data, **kwargs):
    """No filtering
    """
    return data


def GetContacts_edges(pkg, data, GetContacts_threshold, **kwargs):
    """Retain only edges found by GetContacts

    It retains only the raw edges for which `GetContacts` has been able to calculate a
    contact frequency value, i.e., residue pairs that are in physicochemical contact.

    `GetContacts` data is retrieved from the passed `pkg`'s
    :class:`~AlloViz.Protein` object stored in it and thus is needed to have been
    calculated previously. It may (or not) have been calculated using the
    `GetContacts_threshold` kwarg, filtered afterwards with
    :meth:`AlloViz.Wrappers.GetContacts.GetContacts.filter_contacts`, or the
    `GetContacts_threshold` kwarg can also be passed to the analysis to filter the
    `GetContacts` data in this point.

    Other Parameters
    ----------------
    **kwargs
        `GetContacts_threshold` kwarg can be passed to specify the minimum contact
        frequency (0-1, default 0) threshold, which will be used to filter out
        contacts with a frequency (average) lower than it before analysis.
    """
    # Get GetContacts raw data if its calculation is available, else raise an Exception
    try:
        gc = pkg.protein.GetContacts.raw
    except:
        raise Exception("GetContacts results are needed first")
    
    # If GetContacts_threshold kwarg is passed, use it to filter the present data without affecting the saved GetContacts' Protein attribute
    gc = gc[gc["weight"] >= GetContacts_threshold]
        
    # Return the data filtered retaining only the indices that are in GetContacts data
    return data.filter(gc.index, axis=0)


def No_Sequence_Neighbors(pkg, data, Sequence_Neighbor_distance, **kwargs):
    """Filter out residue pairs too close in the sequence

    It only retains edges between residue pairs that are minimum a certain number of
    positions away in the protein sequence (default: 5). It can be assumed that residues
    that are close in the protein sequence will be considerably correlated and/or
    interacting, because of their intrinsic short spatial distance. Specially, it is
    known that a turn of an alpha-helix is completed almost every 4 residues, and that
    the residues i and i+4 (and i-5) in an alpha-helix interact with a backbone hydrogen
    bond. Thus, to take out these interactions and also the contacts due to closeness of
    residues i and i+5, the default value is 5 sequence positions.

    Other Parameters
    ----------------
    **kwargs
        `Sequence_Neighbor_distance` kwarg can be passed to specify the minimum number of
        sequence positions/distance between residues of a pair to retain in Intercontact
        filtering, which defaults to 5.
    """
    # Define a function to retrieve the residue number from the residue nomenclature used (chainID:)RES:resnum
    resnum = lambda res: int(res.rsplit(":")[-1])
    # Return the indices whose residue pair satisfies the Intercontact_distance
    # The absolute value of the subtraction between the two residue numbers must be greater than the parameter to indicate that they are at a greater distance in the sequence than it
    indices = [
        idx
        for idx in data.index
        if abs(resnum(idx[0]) - resnum(idx[1])) > Sequence_Neighbor_distance #Intercontact_distance
    ]
    # Return the filtered data
    return data.filter(indices, axis=0)
    
    
def GPCR_Interhelix(pkg, data, **kwargs):
    """Retain only edges between different TMs(/ECLs/ICLs) of GPCRs

    It can only be used for GPCR systems (`GPCR==True`) for which GPCRdb's generic
    numbering could be retrieved, and it only retains edges between residue pairs that
    have generic numbering and are in different transmembrane helices (and/or
    intra-cellular or extracellular loops) according to it.
    """
    if not pkg.protein.GPCR:
        raise Exception("GPCR_Interhelix filtering is only available for GPCR systems.")
    
    # Create a dictionary mapping the residue numbers to their corresponding TM or ICL/ECL (or 0)
    mapper = dict(zip(
        pkg.protein.protein.select_atoms("name CA").resnums,
        np.floor(pkg.protein.protein.select_atoms("name CA").tempfactors)
    ))
    
    # Define a function to retrieve the residue number from the residue nomenclature used (chainID:)RES:resnum
    resnum = lambda res: int(res.rsplit(":")[-1])
    def are_interhelix(idx):
        # Make a list of the residue pair's TMs
        TMs = [mapper[resnum(idx[0])], mapper[resnum(idx[1])]]
        # Then, check that (1) none of the retrieved TMs are 0 (residue with no generic numbering) and that (2) they are different TMs/ICLs/ECLs 
        return all([
                    all(TMs),
                    len(set(TMs)) > 1
                   ])

    indices = [idx for idx in data.index if are_interhelix(idx)]
    # Return the filtered data
    return data.filter(indices, axis=0)


def Spatially_distant(pkg, data, Interresidue_distance, **kwargs):
    """Retain only edges between spatially distant residue pairs
    
    It only retains edges between residue pairs whose CA atoms are minimum a certain
    number of angstroms away from each other in the initial PDB/structure (default 10 Å).
    The relationship found between these residues can be considered purely allosteric, as
    they are spatially distant and have no direct communication but can be found to be
    interacting/correlated...
    
    Other Parameters
    ----------------
    **kwargs
        'Interresidue_distance' kwarg can be passed to specify the minimum number of
        angstroms that the CA atoms of residue pairs should have between each other in
        the initial PDB/structure (default 10 Å) to be considered spatially distant.
    """
    from MDAnalysis.analysis import distances
    
    # https://userguide.mdanalysis.org/1.1.1/examples/analysis/distances_and_contacts/distances_within_selection.html
    # Create a triangular matrix with all inter-residue distances
    CAs = pkg.protein.protein.select_atoms("name CA")
    n_ca = len(CAs)
    self_distances = distances.self_distance_array(CAs.positions)
    sq_dist_arr = np.zeros((n_ca, n_ca))
    triu = np.triu_indices_from(sq_dist_arr, k=1)
    sq_dist_arr[triu] = self_distances
    
    # Transform the matrix into a pandas DataFrame
    resnames = [f"{aa.resname}:{aa.resid}" for aa in pkg.protein.protein.residues]
    df = pandas.DataFrame(sq_dist_arr, columns=resnames, index=resnames)
    df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
    df = pandas.DataFrame({"dist": df.stack()})
        
    indices = df[df["dist"] >= Interresidue_distance].index
    # Return the filtered data
    return data.filter(indices, axis=0)



class NoNetworkException(Exception): pass


class Filtering:
    r"""Class for network filtering

    Instances of this class are used to be added as attributes to instances of
    Wrappers' classes (as attributes of a :class:`AlloViz.Protein` object) with the
    purpose of storing (un)filtered networks according to different criteria.

    Raw network edges stored as a DataFrame in the passed `pkg` object are filtered with
    the corresponding function or combination of functions of the
    :mod:`~AlloViz.AlloViz.Filtering` module and converted to individual NetworkX'
    :external:ref:`Graphs <graph>` with :meth:`~AlloViz.AlloViz.Analysis.Analysis._get_G`
    to be stored as private attributes.

    Parameters
    ----------
    pkg : instance of a class from :mod:`AlloViz.Wrappers`
        The object is used to retrieve the raw network edges for network analysis. It
        also contains/gives access to the corresponding:class:`~AlloViz.Protein` object
        and its information.
    filtering : str or list
        Filtering scheme(s) with which to filter the list of network edges. A list of
        strings is used to filter with a combination of criteria. All available (and 
        combinable) filtering options are functions in the
        :mod:`~AlloViz.AlloViz.Filtering` module: 
        :func:`~AlloViz.AlloViz.Filtering.All`,
        :func:`~AlloViz.AlloViz.Filtering.GetContacts_edges`,
        :func:`~AlloViz.AlloViz.Filtering.No_Sequence_Neighbors`,
        :func:`~AlloViz.AlloViz.Filtering.GPCR_Interhelix`,
        :func:`~AlloViz.AlloViz.Filtering.Spatially_distant`.
    name : str
        Name of the filtering scheme. It will be the same name as the passed filtering
        option if it is a single string, or the names of the filtering options joined by
        "_" if a lsit has been passed. It is used to name the pkg's attribute in which
        the class instance is saved.

    Other Parameters
    ----------------
    GetContacts_threshold : float
        Optional kwarg that can be passed to specify the minimum contact
        frequency (0-1, default 0) threshold, which will be used to filter out
        contacts with a frequency (average) lower than it before analysis.
    Sequence_Neighbor_distance : int
        Optional kwarg that can be passed to specify the minimum number of sequence
        positions/distance between residues of a pair to retain in No_Sequence_Neighbors
        filtering, which defaults to 5.
    Interresidue_distance : int or float
        Optional kwarg that can be passed to specify the minimum number of angstroms
        that the CA atoms of residue pairs should have between each other in the initial
        PDB/structure (default 10 Å) to be considered spatially distant.
         
    Attributes
    ----------
    graphs : dict of external:ref:`Graph <graph>` objects

    See Also
    --------
    AlloViz.Protein.filter : Class method to filter the network(s) raw edge weights with
                             different criteria.
    AlloViz.Wrappers.Base.Base.filter : Pkg's method to filter the network(s) raw edge
                                        weights with different criteria.
    """

    def __init__(self, pkg, filtering, name, *, GetContacts_threshold=0, Sequence_Neighbor_distance=5, Interresidue_distance=10):
        self._pkg = pkg
        self._name = name

        # Establish the paths and filenames
        self._path = f"{self._pkg.protein._datadir}/{self._pkg._name}/{self._name}"
        os.makedirs(self._path, exist_ok=True)
        self._datapq = lambda element, metric: f"{self._path}/{element}_{metric}.pq"

        # Get the filtered data according to the filtering scheme(s)
        filterings = (
            filtering
            if isinstance(filtering, list)
            else [filtering]
        )
        
        data = self._pkg.raw
        for filt in filterings:
            filtfunc = eval(filt)
            data = filtfunc(self._pkg, data, 
                            GetContacts_threshold=GetContacts_threshold, 
                            Sequence_Neighbor_distance=Sequence_Neighbor_distance,
                            Interresidue_distance=Interresidue_distance)
        # Drop all-0 rows (not taking into account weight_std column if it's present)
        self._filtdata = data.drop(data[(data.drop(columns=[c for c in data.columns if "std" in c]) == 0).all(axis=1)].index, axis=0)
        
        self._graph_distances = -np.log(abs(self._filtdata) + 10E-10) + 10E-10
        # un-transform standard error columns
        self._graph_distances.loc[:,["std" in c for c in self._graph_distances.columns]] = data.loc[:,["std" in c for c in data.columns]]
        
        # For each column in the filtered data that is not a standard error, create an analyzable NetworkX's graph and save it
        self.graphs = {}
        for col in [c for c in self._graph_distances if "std" not in c]:
            try:
                # The approach in lit. is to use -log10(|corr|) as edge weights/distances in the network for analyses
                # e.g., https://www.pnas.org/doi/full/10.1073/pnas.0810961106            
                self.graphs[col] = self._get_G(self._graph_distances[col])
            except NoNetworkException as e:
                print(e)
                
                

    def _get_G(self, column):
        r"""Return a DataFrame's column as a Graph

        Transform a column from a :class:`pandas.DataFrame` (passed as a
        :class:`pandas.Series`) into a NetworkX' :class:`~networkx.Graph`; retaining only the
        largest component if the passed information results in an unconnected network.

        Parameters
        ----------
        column : :class:`pandas.Series`
        """
        # Transform the column into a 3-column dataframe with the nodes' name and the value as weight. Drop 0s, NAs and use absolute value
        weights = (
            column[column != 0].dropna().abs().rename("weight").reset_index()
        )  # btw calculations fail with 0 value weights and cfb prob with negative
        # Create the NetworkX's Graph and a list of its components
        network = networkx.from_pandas_edgelist(weights, "level_0", "level_1", "weight")
        components = list(networkx.connected_components(network))
        
        if len(components) == 0:
            raise NoNetworkException(
                " ".join((
                    f"EXCEPTION! No connected components in network ({network.number_of_nodes()} nodes):",
                    self._pkg._name,
                    self._name,
                    column.name,
                ))
            )
        else:
            # Check that the largest connected component has the same size as the total number of nodes, else select the largest component to return as the Graph to be analyzed
            largest_component = max(components, key=len)

            if len(largest_component) < network.number_of_nodes():
                print(
                    f"WARNING! Unconnected network ({network.number_of_nodes()} nodes):",
                    self._pkg._name,
                    self._name,
                    column.name,
                    "\n",
                    f"Largest network component will be used for analysis. Sizes (number of nodes) of all components: {[len(comp) for comp in components]}",
                )
                network = network.subgraph(largest_component)

            return network
    
    def analyze(self, elements="edges", metrics="all", cores=1, nodes_dict=Analysis.nodes_dict, edges_dict=Analysis.edges_dict, **kwargs):
        r"""Analyze the filtered network
        
        Send the analyses of the filtered network for the passed combinations of
        elements-metrics. The individual :external:ref:`Graphs <graph>` saved as private
        attributes upon object initialization are analyzed independently with the
        :mod:`~AlloViz.AlloViz.Analysis` module (this function calls
        :func:`AlloViz.AlloViz.Analysis.analyze`). Results are stored as new instances
        of classes from the :mod:`AlloViz.AlloViz.Elements` module, which extend the
        :class:`pandas.DataFrame` class.

        Parameters
        ----------
        elements : str or list, {"edges", "nodes"}
            Network element for which to perform the analysis.
        metrics : str or list, default: "all"
            Network metrics to compute, which must be keys in the `nodes_dict` or
            `edges_dict` dictionaries. Default is "all" and it sends the computation for
            all the metrics defined in the corresponding dictionary of the selected
            elements in `element`.
        cores : int, default: 1
            Number of cores to use for parallelization with a `multiprocess` Pool.
            Default value only uses 1 core with a custom :class:`AlloViz.utils.dummypool`
            that performs computations synchronously.

        Other Parameters
        ----------------
        nodes_dict, edges_dict : dict
            Optional kwarg(s) of the dictionary(ies) that maps network metrics custom 
            names (e.g., betweenness centrality, "btw") with their corresponding NetworkX
            function (e.g., "networkx.algorithms.centrality.betweenness_centrality").
            Functions strings must be written as if they were absolute imports, and must
            return a dictionary of edges or nodes, depending on the element dictionary in
            which they are. The keys of the dictionaries will be used to name the columns
            of the analyzed data that the functions produce. Defaults are
            :data:`~AlloViz.AlloViz.Analysis.nodes_dict` and
            :data:`~AlloViz.AlloViz.Analysis.edges_dict`.
        **kwargs
            Other optional keyword arguments that will be passed to the NetworkX analysis
            function(s) that is(are) used on the method call in case they need extra
            parameters.
        """
        # Depending on the desired cores, use a dummypool (synchronous calculations) or a `multiprocess` Pool
        # Changing it inside the `utils` module allows to share the same one between modules
        if cores > 1:
            mypool = Pool(cores)
        else:
            mypool = utils.dummypool()
        utils.pool = mypool
        
        if self._filtdata.size == 0:
            print(f"{self._pkg._name} {self._name} is not a connected network (or subnetwork)")
        else:
            Analysis.analyze(self, elements, metrics, nodes_dict, edges_dict, **kwargs)
        
        # Close the pool
        utils.pool.close()
        utils.pool.join()
        utils.pool = utils.dummypool()