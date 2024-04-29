"""Base module for package/software wrapping for AlloViz

Main :class:`~AlloViz.Wrappers.Base.Base` class stores the information about the
:class:`~AlloViz.Protein` object and defines the general private methods that most
child classes use to launch, manage and store the network calculations. It exploits the
'__new__' special method to establish the object's attributes and that can be extended in
the child classes and '__init__' is used to launch calculations, which can be extended or
overriden by child classes.

In this module there are additional classes that override or extend 
:class:`~AlloViz.Wrappers.Base.Base`'s private methods and that are used to be combined
to create child classes with multiple inheritance for specific uses.

"""

import os, time

import pandas
import numpy as np

from contextlib import redirect_stdout, redirect_stderr
from importlib import import_module
from lazyasd import LazyObject

from ..AlloViz.Filtering import Filtering
from ..AlloViz.Elements import Edges
from ..AlloViz.utils import get_pool, rgetattr, rhasattr
from ..AlloViz.info import citations
from ..AlloViz import utils




def lazy_import(key, val):
    """Return a lazily-imported module to store in a variable
    
    It uses `lazyasd <https://github.com/xonsh/lazyasd>`_.

    Parameters
    ----------
    key : str
        The variable name with which the module is going to be accessed. In 
        `import numpy as np` "key" would be `np`.
    val : str
        Absolute import of the module.
    """
    # If the import is inside AlloViz's 'Packages' folder (src/Packages in the repo and AlloViz/Packages in the installation), this extra "package" argument is needed
    extra_arg = {"package": 'AlloViz'} if 'Packages' in val else {}
    # Return the LazyObject that is going to be stored in a variable.
    return LazyObject(lambda: import_module(val, **extra_arg), globals(), key)







class Base:
    """Base class for network calculation
    
    This class uses the '__new__' special method to establish the object's attributes,
    which can be extended and/or modified in the child classes' '__new__'. The '__init__'
    special method is then used to launch calculations, and it can also be extended or
    overriden by child classes.

    It also defines the private methods :meth:`~AlloViz.Wrappers.Base.Base._calculate` to
    launch calculations, which exploits the child classes' custom '_computation' private
    method so that the different packages can be run with them using the standardized
    information stored in the objects attributes; and
    :meth:`~AlloViz.Wrappers.Base.Base._save_pq` to save the results of the calculations
    in parquet format. Both can also be extended or overriden by child classes.
    
    :meth:`~AlloViz.Wrappers.Base.Base.filter` allows to filter the calculated raw edges
    for posterior analysis.
    
    Parameters
    ----------
    protein : :class:`~AlloViz.Protein` object
        Protein object from which to extract the information for calculation.    
    d : dict
        Dictionary containing all the elements of the protein's `__dict__` and also any
        keyword arguments passed for calculation. This is needed for successful
        parallelization, as pickling objects of Base's child classes means pickling the
        Protein object and the attributes that are added to it during its `__init__` are
        lost.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls)
        
        new._name = new.__class__.__name__
        new.protein = protein
        new._d = d
        
        new._pdbf = d["_pdbf"]
        new._trajs = d["_trajs"]
        
        new._path = f"{d['_datadir']}/{new._name}/raw"
        os.makedirs(new._path, exist_ok=True)
        new._rawpq = lambda xtc: f"{new._path}/{xtc if isinstance(xtc, int) else xtc.rsplit('.', 1)[0]}.pq"
        
        return new
    
    
    def __getnewargs__(self):
        # This is necessary for successful parallelization/pickling
        return self.protein, self._d
        
        
    def __init__(self, *args):
        # Define the list of .pq files that we expect are going to be saved (or be retrieved) and a function to check which of them already exist
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        
        # If any of the .pq files don't exist, send the calculation for those that don't; then continue to read the saved data
        if any(no_exist(pqs)):
            for xtc in (xtc for xtc in self._trajs if no_exist(pqs)[xtc-1]):
                self._calculate(xtc)
        
        
        # Function to wait for the calculations to finish in the background; returns the .pq files to be read and added as attributes when they do
        def wait_calculate(pqs):
            while any(no_exist(pqs)):
                time.sleep(5)
            return pqs
        
        # Function to read the newly calculated (or retrieved) data files, concatenate them, calculate averages and std if necessary and return it to be added as raw data attribute
        def get_raw(pqs):
            print(f"adding raw data of {self._name} for {self._pdbf}: ", pqs)
            # List of the corresponding .pq files to concatenate
            flist = map(lambda pq: pandas.read_parquet(pq), pqs)
            df = pandas.concat(flist, axis=1)
            
            # Perform min-max normalization (to range 0-1) on data coming from AlloViz or MDEntropy related methods (their MI values are low) or PyInteraph2_Atomic_Contacts_Strength
            if "AlloViz" in self._name or "MDEntropy" in self._name or "PyInteraph2_Atomic_Contacts_Strength" in self._name:
                df = (df-df.min())/(df.max()-df.min())
            
            # If there are more than 1 trajectory, calculate the average and standard error of the trajectories' raw data
            if len(self._trajs) > 1:
                cols = [f"{num}" for num in self._trajs]
                df["weight"] = df[cols].fillna(0).mean(axis=1)
                df["weight_std"] = df[cols].fillna(0).std(axis=1)
            else:
                # If there is only 1 trajectory, simply use its raw data as the "weight" variable
                df.rename(columns={"1": "weight"}, inplace=True)
                
            return Edges(df, parent=self.protein)
        
        # Function to call get_raw to read and process the raw data and add it as the "raw" attribute to self
        add_raw = lambda pqs: setattr(self, "raw", get_raw(pqs))
        # Wait asynchronously for analysis to end and then add the data calling add_raw
        get_pool().apply_async(wait_calculate,
                               args=(pqs,),
                               callback=add_raw)
        
        
        pkg = [name for name in citations if name in self.__class__.__name__]
        if len(pkg)==1 and "AlloViz" not in pkg:
            print(f"Please, make sure to correctly cite the package used to compute the network: {pkg[0]} ({citations[pkg[0]]})")
        
        
    
    def _calculate(self, xtc, *args):
        """Send the calculation for a single trajectory file
        
        Send the calculation with the specific class' `_computation` private method and
        capture the standard output and standard error of the calculation into a log
        file. The callback function to which the reuslts of the computation are passed is
        :meth:`~AlloViz.Wrappers.Base.Base._save_pq` to save the returned calculated
        results. \*args parameter is not used.
        
        Parameters
        ----------
        xtc : int
            Trajectory number from the processed trajectory dictionary of the Protein
            object.
        """
        # Function to send the _computation and capture stdout and stderr to a log file
        def send_and_log(xtc, *args):
            with open(f"{self._path}/{self._name}.log", "a+") as f:
                with redirect_stdout(f), redirect_stderr(f):
                    return self._computation(xtc, *args)
        
        # Send the computation of the trajectory file to the pool with _save_pq as callback
        get_pool().apply_async(send_and_log,
                         args=(xtc, *args),
                         callback=self._save_pq)
        
    
    
    
    def _save_pq(self, args):
        """Save the calculation results in a parquet file
        
        It is used as callback function after the calculation with a specific class'
        `_computation` private method to save the calculated data of a single trajectory
        file in tabular format in a parquet (.pq) file.
        
        Parameters
        ----------
        args : list
            List of elements returned by the `_computation` function. It will always
            contain a `corr` matrix (numpy.ndarray) and the `xtc` trajectory number. It
            may also include a residue list `resl` with the residue indices out of all
            the residues in the Protein object for which values were calculated. This is
            used, e.g., for calculations involving dihedral angles, in which residues in
            the extremes will not have some dihedrals and thus values won't be calculated
            for them.
        """
        corr, xtc, *resl = args
        
        # If resl is passed, check that corr has the appropriate shape or select the values using resl if not
        if len(resl) != 0:
            resl = resl[0]
            if corr.shape != (len(resl), len(resl)):
                corr = corr[np.ix_(resl, resl)]
        # Else, simply make resl a slice that will select data of the same size as corr (i.e., all data)
        elif len(resl) == 0:
            resl = slice(0, corr.shape[0])
        
        # Make a list of residue names selecting the appropriate residues with resl
        resnames = [f"{aa.resname}:{aa.resid}" for aa in self._d["protein"].residues[resl]]
        
        # Transform the square, symmetric corr matrix into a DataFrame and select only the upper triangle (without the diagonal: k=1) and stack it into a tabular format for saving
        df = pandas.DataFrame(corr, columns=resnames, index=resnames)
        df = df.where( np.triu(np.ones(df.shape), k=1).astype(bool) )
        df = pandas.DataFrame({f"{xtc}": df.stack()})
        # if not len(df[f"{xtc}"].unique()) == 1:
        df.to_parquet(self._rawpq(xtc))
    
    
    
    
    def filter(self, filterings="All", **kwargs):
        r"""Filter network edges
        
        Filter the networks according to the selected criteria to perform analyses on
        (all or) a subset of the edges. It calls
        :meth:`AlloViz.Wrappers.Base.Base.filter` and results are stored in instances
        of the :class:`AlloViz.AlloViz.Filtering.Filtering` class. The different filtering
        options are detailed in the :mod:`~AlloViz.AlloViz.Filtering` module.

        Parameters
        ----------
        filterings : str or list of strs and/or lists, default: "all"
            Filtering scheme(s) with which to filter the list of network edges before
            analysis. It can be a string, or a list of strings and/or lists: a list of
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
        >>> opioidGPCR.dynetan.filter(["GetContacts_edges", ["GetContacts_edges", "GPCR_Interhelix"]])
        >>> opioidGPCR.dynetan.GetContacts_edges_GPCR_Interhelix
        <AlloViz.AlloViz.Filtering.Filtering at 0x7f892c3c0fa0>
        """

        for filt in filterings:
            # Name used to store as attribute will be that of the filtering scheme chosen or the combination's names joined by "_"
            name = filt if isinstance(filt, str) else "_".join(filt)
            if not rgetattr(self, name):
                setattr(
                    self,
                    name,
                    Filtering(self, filt, name, **kwargs),
                )
        
        return rgetattr(self, name) if len(filterings) == 1 else None
        

    
    
    
# class Use_COM(Base):
#     """Class for using the COM's structure file and trajectories
    
#     Classes that inherit this class use the residues' COM structure and trajectory(ies)
#     files for calculations instead of the whole protein's.
#     """
    
#     def __new__(cls, protein, d):
#         new = super().__new__(cls, protein, d)
        
#         new._pdbf = new._d["_compdbf"]
#         new._trajs = new._d["_comtrajs"]
        
#         return new
    
    
    


class Multicore(Base):
    """Class for multi-core packages calculations
    
    This class defines additional attributes 'taskcpus' and '_empties' in '__new__' for
    packages that are able to perform multi-core calculations by themselves (besides
    AlloViz's parallelization of the calculations for each trajectory file). 'taskcpus'
    is the number of cores the package should use per trajectory and '_empties' the
    number of empty jobs that should be sent to the same Pool that the task/_computation
    is being sent to to "occupy" the number of cores that the package is going to use in
    reality, as sending it to the Pool would only take 1 of its workers. It extends the
    :meth:`~AlloViz.Wrappers.Base.Multicore._calculate` private method to do it.
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        # If taskcpus is not specified, set it to half of the available cores in the system
        if "taskcpus" not in new._d:
            new.taskcpus = int(np.ceil(os.cpu_count()/2))
        else:
            new.taskcpus = new._d["taskcpus"]
        
        # The actual _computation job will only occupy 1 worker of the Pool, so _empties should be the rest of desired taskcpus
        new._empties = new.taskcpus-1
            
        return new
    
    
    
    def _calculate(self, xtc, *args):
        """Extend the function to send empty jobs and finish them when the main job saves the data
        
        """
        super()._calculate(xtc, *args)
        
        
        def calculate_empty(pq):
            # print("sleeping", pq, os.getpid())
            while not os.path.isfile(pq):
                time.sleep(5)
            return
        
        for _ in range(self._empties): get_pool().apply_async(calculate_empty, args=(self._rawpq(xtc),))
    
    
    



        
class Combined_Dihs(Base):
    """Class for combination of dihedral angle data
    
    This class' child classes are used to combine the information from multiple dihedral
    angles by overriding the `_calculate` and `_save_pq` private methods. It checks
    that the calculations of the desired dihedral angles for the current trajectory
    number (`xtc`) are available, or else raises an error, and saves the data with the
    child class' `_save_pq` specific private method.
    """
    
    def _calculate(self, xtc):
        # Child classes' names are: correlationplus_Backbone_Dihs_Avg, AlloViz_Sidechain_Dihs_Max, etc... 
        pkg = self._name.rsplit("_", 2)[0] if ("Sidechain" in self._name or "Backbone" in self._name) else self._name.rsplit("_", 1)[0]
        
        # If any of the dihedral calculations don't exist, raise error
        attrs_exist = {f"{pkg}_{Dih}": rhasattr(self, "protein", f"{pkg}_{Dih}") for Dih in self._dihs}
        if any([not attr_exist for attr_exist in attrs_exist.values()]):
            raise Exception(
                f"Individual dihedrals calculations are needed first: {attrs_exist}"
            )
            
        # Function to get the name of the files that we aim to retrieve
        get_rawpq = lambda Dih: rgetattr(self, "protein", f"{pkg}_{Dih}", "_rawpq")(xtc)
        pqs = [get_rawpq(Dih) for Dih in self._dihs]
        
        self._save_pq(pqs, xtc)

        
        
        
class Combined_Dihs_Avg(Combined_Dihs):
    """Class for combination of dihedral angle data by averaging
    
    This class' child classes are used to combine the information from multiple dihedral
    angles by averaging, with its specific
    :meth:`~AlloViz.Wrappers.Base.Combined_Dihs_Avg._save_pq` private method.
    """
    
    def _save_pq(self, pqs, xtc):
        # Make a list of the trajectories' dihedrals computations results (absolute numbers)
        dfs = [pandas.read_parquet(pq)[f"{xtc}"].abs() for pq in pqs]
        # Concatenate them, drop edges for which none of the columns have a value (sanity check) and calculate the rowwise means
        avgs = pandas.concat(dfs, axis=1).dropna(how="all").mean(axis=1)
        # Save the averages as the final data
        pandas.DataFrame({f"{xtc}": avgs}).to_parquet(self._rawpq(xtc))

        
        
# class Combined_Dihs_Max(Combined_Dihs):
#     """Class for combination of dihedral angle data by taking the maximum value
    
#     This class' child classes are used to combine the information from multiple dihedral
#     angles by taking for each edge the maximum edge weight of all the dihedral networks
#     that are combined, with its specific
#     :meth:`~AlloViz.Wrappers.Base.Combined_Dihs_Max._save_pq` private method.
#     """
    
#     def _save_pq(self, pqs, xtc):
#         # Make a list of the trajectories' dihedrals computations results (absolute numbers)
#         dfs = [pandas.read_parquet(pq)[f"{xtc}"].abs() for pq in pqs]
#         # Concatenate them, drop edges for which none of the columns have a value (sanity check) and retrieve the rowwise max values
#         maxs = pandas.concat(dfs, axis=1).dropna(how="all").max(axis=1)
#         # Save the max values as the final data
#         pandas.DataFrame({f"{xtc}": maxs}).to_parquet(self._rawpq(xtc))