"""AlloViz's own network construction method wrapper

It calculates the Mutual Information (MI) generalized correlation of the residues'
dihedral angles and also their combinations.

"""

import numpy as np

from .Base import lazy_import, Multicore, Combined_Dihs_Avg, Combined_Dihs_Max

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




    
class AlloViz_Base(Multicore):
    """AlloViz network construction method base class

    :meth:`~AlloViz.Wrappers.AlloViz.AlloViz_w._computation` uses
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
        # Establish the protein atoms from the Universe with the trajectories
        prot = self._d["u"].atoms
        
        # Establish the atoms that form the desired dihedral (if they are not the default ones from the selection function and are provided) 
        # and the selection function itself, and then retrieve the selected AtomGroups
        dih_atoms = self._dih_atoms if hasattr(self, "_dih_atoms") else {}
        # select_dih = lambda res: eval(f"res.{self._dih.lower()}_selection(**dih_atoms)")
        # selected = [select_dih(res) for res in prot.residues]
        # selected = eval(f"prot.residues.{self._dih.lower()}_selections(**dih_atoms)")
        selected = []
        for res in prot.residues:
            ag = eval(f"res.{self._dih.lower()}_selection(**dih_atoms)")
            selected.append(ag)
        
        # Save the residue indices of the protein for which an AtomGroup forming the desired dihedral could be found, and also make a list of the AtomGroups without Nones
        selected_res, selected = zip(*[(i, ag) for i, ag in enumerate(selected) if ag])
        # selected_res =  d["_dihedral_residx"]() # Used before for backbone dihedrals; from the Protein's _dihedral_residx
        
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
        
        # Establish the striding if it was passed for calculation
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
            
            
        return corr, xtc, list(selected_res)
    
    
    
    
    
class AlloViz_Psi(AlloViz_Base):
    """AlloViz network construction method's MI of Psi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "psi"
        return new

    
    
    
class AlloViz_Phi(AlloViz_Base):
    """AlloViz network construction method's MI of Phi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "phi"
        return new

        

class AlloViz_Omega(AlloViz_Base):
    """AlloViz network construction method's MI of Omega backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "omega"
        return new

        
        
        
    
class AlloViz_Backbone_Dihs(Combined_Dihs_Avg, AlloViz_Base):
    """AlloViz network construction method's of the combination of the backbone
    dihedrals' MIs by averaging
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dihs = ["Phi", "Psi"]#, "Omega"]
        return new
    
class AlloViz_Backbone_Dihs_Max(Combined_Dihs_Max, AlloViz_Base):
    """AlloViz network construction method's of the combination of the backbone
    dihedrals' MIs by taking the max value
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dihs = ["Phi", "Psi"]#, "Omega"]
        return new


    




class AlloViz_Chi1(AlloViz_Base):
    """AlloViz network construction method's MI of Chi1 side-chain dihedral
    
    GLY and ALA residues don't have the Chi1 side-chain dihedral.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "chi1"
        return new
    
    
    
class AlloViz_Chi2(AlloViz_Base):
    """AlloViz network construction method's MI of Chi2 side-chain dihedral
    
    In addition to GLY and ALA, CYS, SER, THR and VAL residues don't have the Chi2
    dihedral.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "chi1"
        # dihedral atomnames from: http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
        new._dih_atoms = {"n_name": "CA",
                          "ca_name": "CB",
                          "cb_name": "CG CG1",
                          "cg_name": "CD OD1 ND1 CD1 SD"}
        return new

                          
                          
class AlloViz_Chi3(AlloViz_Base):
    """AlloViz network construction method's MI of Chi3 side-chain dihedral
    
    Chi3 side-chain dihedral is only available for ARG, GLN, GLU, LYS and MET residues.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "chi1"
        new._dih_atoms = {"n_name": "CB",
                          "ca_name": "CG",
                          "cb_name": "CD SD",
                          "cg_name": "NE OE1 CE"}
        return new

    
class AlloViz_Chi4(AlloViz_Base):
    """AlloViz network construction method's MI of Chi4 side-chain dihedral
    
    Chi4 side-chain dihedral is only available for ARG and LYS residues.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "chi1"
        new._dih_atoms = {"n_name": "CG",
                          "ca_name": "CD",
                          "cb_name": "CE NE",
                          "cg_name": "CZ NZ"}
        return new
    
    
class AlloViz_Chi5(AlloViz_Base):
    """AlloViz network construction method's MI of Chi5 side-chain dihedral
    
    Chi5 side-chain dihedral is only available for ARG residues.
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "chi1"
        new._dih_atoms = {"n_name": "CD",
                          "ca_name": "NE",
                          "cb_name": "CZ",
                          "cg_name": "NH1"}
        return new
    
    
    
    
class AlloViz_Sidechain_Dihs(Combined_Dihs_Avg, AlloViz_Base):
    """AlloViz network construction method's of the combination of the side-chain
    dihedrals' MIs by averaging
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = [f"Chi{i+1}" for i in range(new.chis)]
        return new
    
    
class AlloViz_Sidechain_Dihs_Max(Combined_Dihs_Max, AlloViz_Base):
    """AlloViz network construction method's of the combination of the side-chain
    dihedrals' MIs by taking the max value
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = [f"Chi{i+1}" for i in range(new.chis)]
        return new
    
    
    
    
    
class AlloViz_Dihs(Combined_Dihs_Avg, AlloViz_Base):
    """AlloViz network construction method's of the combination of all backbone and
    side-chain dihedrals' MIs by averaging
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = ["Phi", "Psi"] + [f"Chi{i+1}" for i in range(new.chis)]#, "Omega"]
        return new
    
    
class AlloViz_Dihs_Max(Combined_Dihs_Max, AlloViz_Base):
    """AlloViz network construction method's of the combination of all backbone and
    side-chain dihedrals' MIs by taking the max value
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = ["Phi", "Psi"] + [f"Chi{i+1}" for i in range(new.chis)]#, "Omega"]
        return new