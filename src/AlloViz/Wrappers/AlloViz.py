"""AlloViz's own network construction method wrapper

It calculates the Mutual Information (MI) generalized correlation of the residues'
backbone dihedral angles (Phi, Psi and Omega) and also their average.

"""

import numpy as np

from .Base import lazy_import, Multicore, Combined_Dihs

imports = {
"_npeet_lnc": ".Packages.NPEET_LNC.lnc",
"_mda_dihedrals": "MDAnalysis.analysis.dihedrals",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
def _calculate_row(in_data, values, nres):
    """"""
    res, others = in_data
    row = np.zeros((nres,))
    
    # Populate the row with the MI values between the present residue and each of the others
    for other in others:
        row[other] = _npeet_lnc.MI.mi_LNC([values[res], values[other]])
        
    return res, row




    
class AlloViz(Multicore):
    """AlloViz network construction method base class

    :meth:`~AlloViz.Wrappers.AlloViz.AlloViz._computation` uses
    :mod:`MDAnalysis.analysis.dihedrals` to extract data and 
    `NPEET_LNC <https://github.com/ViktorvdValk/NPEET_LNC>`_ for MI calculation.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        if "stride" in d:
            new.stride = d["stride"]
        return new
    
    def _computation(self, xtc):
        """"""
        # Establish the protein atoms and the indices of the residues included in the calculation (from the Protein's _dihedral_residx)
        prot = self._d["u"].atoms
        selected_res = self._d["_dihedral_residx"]()
        
        # Establish a function to retrieve the AtomGroups forming a certain dihedral angle in a certain residue, and get the selections for all the selected residues
        select_dih = lambda res: eval(f"res.{self._dih.lower()}_selection()")
        selected = [select_dih(res) for res in prot.residues[selected_res]]
        
        # Depending on the trajectory number, establish the start and stop frames to get dihedral angle values time series from the Universe
        offset = 0
        # If there is more than 1 trajectory file, take into account the number of frames of all the trajectories previous to the present file
        if hasattr(self._d["u"].trajectory, "readers"):
            get_frames = lambda numreader: self._d["u"].trajectory.readers[numreader].n_frames
            for numreader in range(xtc-1):
                offset += get_frames(numreader)
        # Else, the "stop" frames (offset+get_frames) is just the number of frames of the trajectory
        else:
            get_frames = lambda _: self._d["u"].trajectory.n_frames

        step = {"step": self.stride} if hasattr(self, "stride") else {}
        values = _mda_dihedrals.Dihedral(selected).run(start=offset, stop=offset+get_frames(xtc-1), **step).results.angles.transpose()
        
        # Create a numpy array to store results
        corr = np.zeros((len(selected_res), len(selected_res)))
        print("prot.n_residues:", prot.n_residues, "len(selected_res)", len(selected_res), "dihedrals array shape:", values.shape, "empty corr array shape:", corr.shape)
        
        # Calculate MIs: https://stackoverflow.com/q/72783941
        in_data = []
        iterator = set(range(len(selected_res)))
        # For each residue (row in the matrix), add a tuple to in_data with its id and the id of the rest of residues with a "higher" id
        for res1 in iterator:
            # "others" includes all the other residues except itself and residues "below" it, as the MI between them and the present one have already been calculated (triangular matrix)
            others = list(iterator - set(range(res1)) - {res1})
            # Sanity check to avoid adding the last residue, which will have no residues "above" it
            if len(others) > 0:
                in_data.append((res1, others))
        
        # Iterate over "in_data" (one calculation per residue: one calculation per row with _calculate_row) and store each row data in "results"
        from concurrent.futures import ProcessPoolExecutor as Pool
        from functools import partial
        with Pool(self.taskcpus) as p:
            results = list(p.map(partial(_calculate_row, values=values, nres=len(selected_res)), in_data))
            p.shutdown()
        
        # Populate corr with results
        for result in results:
            res, row = result
            corr[res, :] = row
            
            
        return corr, xtc, selected_res
    
    
    
    
    
class AlloViz_Psi(AlloViz):
    """AlloViz network construction method's MI of Psi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "psi"
        return new

    
    
    
class AlloViz_Phi(AlloViz):
    """AlloViz network construction method's MI of Phi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "phi"
        return new

        

class AlloViz_Omega(AlloViz):
    """AlloViz network construction method's MI of Omega backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "omega"
        return new

        
        
        
    
class AlloViz_Dihs(Combined_Dihs, AlloViz):
    """AlloViz network construction method's of the combination of the backbone dihedrals' MIs
    """
    _Phi = AlloViz_Phi
    _Psi = AlloViz_Psi
    _Omega = AlloViz_Omega
    pass