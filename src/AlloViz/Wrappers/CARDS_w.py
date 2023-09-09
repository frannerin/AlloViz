"""CARDS network construction method wrapper

It calculates the Mutual Information (MI) generalized correlation of the residues'
dihedral angles, and also the MI of their disorder, the MI between their values and their
disorder (pure disorder or disorder-mediated) and the holistic MI (MI of the values plus
of the disorder), and also their combinations.

"""

import pandas, os, time
import numpy as np

from .Base import lazy_import, Base, Multicore, Combined_Dihs_Avg

imports = {
"_enspara": "enspara.apps.collect_cards",
"_cards": "enspara.cards",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
class CARDS(Multicore):
    """CARDS network construction method base class for calculations
    """
    
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        new.buffer = d["CARDS_buffer"] if "CARDS_buffer" in d else 15
        # buffer_width: int, default=15
        # The width of the no-man's land between rotameric bins. Angles
        # in this range are not used in the calculation.
        
        return new
    
    
    def __init__(self, *args):
        super().__init__(*args)
        
        # Wait for the calculations to finish before finishing class instance creation; calcs are needed to proceed with other CARDS Wrappers
        pqs = [self._rawpq(xtc) for xtc in self._trajs]
        no_exist = lambda pqs: [not os.path.isfile(pq) for pq in pqs]
        while any(no_exist(pqs)):
            time.sleep(5)
    
    
    def _computation(self, xtc):
        """"""
        class args: pass
        args.trajectories = [[self._trajs[xtc]]]
        args.topology = [self._pdbf]
        
        result = _cards.cards(_enspara.load_trajs(args), self.buffer, self.taskcpus)
        
        return result, xtc
        
        
        
    def _save_pq(self, args):
        result, xtc = args
        ss_mi, dd_mi, sd_mi, ds_mi, inds = result        
        disorder_mediated = sd_mi + ds_mi + dd_mi
        holistic = disorder_mediated + ss_mi
        # inds has shape (n_dihedrals, 4 atoms) and the MI matrices have shapes (n_dihs, n_dihs); i.e., indexing inds with [i, :] translates to indexing MIs with [i, i]
        names = ("MI",  "Disorder", "Disorder_mediated", "Holistic")
        MIs = list(zip((ss_mi, dd_mi, disorder_mediated, holistic), names))
        
        # Translate the array of atom ids that form the dihedrals into an array of the protein's MDAnalysis' AtomGroups 
        atoms = self._d["protein"].atoms[inds]
        # Translate it into an array of the atomnames and get the unique atomnames combinations to identify the dihedrals
        atomnames = atoms.names
        unqs = set(map(tuple, atomnames))
        # Also translate it into an array of the residue index of each of the atoms of each dihedral to assign each dihedral to a residue index (the majority one)
        resids = list(map(list, (atom.resindex for atom in atoms)))
        resids = np.array([max(l, key = l.count) for l in resids])
        
        # Initialize a list for each of the dihedrals to store the atomname combinations found previously that define each of them
        # Based on AlloViz.AlloViz.info.dihedrals_atoms dictionary, and mostly can be done looking at the first atom that defines the dihedral
        phi, psi, chi1, chi2, chi3, chi4 = ([] for i in range(6))
        for unq in unqs:
            i = unq[0]
            if i == "C": phi.append(unq)
            elif i == "N":
                if unq[-1] == "N": psi.append(unq)
                elif "G" in unq[-1]: chi1.append(unq)
            elif i == "CA": chi2.append(unq)
            elif i == "CB": chi3.append(unq)
            elif i == "CG": chi4.append(unq)
            else: raise Exception(unq)
        
        # Define a Base class without an init to define children for each of the dih-MI combination to use its _save_pq without launching the __init__ and thus _calculate
        class CARDS_Base(Base):
            def __init__(self, *args):
                pass
            
        # For each dihedral, retrieve the indices (e.g., of inds' first dimension) of the data that correspond to it according to the combinations of atomnames retrieved
        for dih in ["phi", "psi", "chi1", "chi2", "chi3", "chi4"]:
            indices = np.zeros(atomnames.shape[0], dtype=bool) # All false 1D array
            for dih_atoms in eval(dih): # For each combination of atomnames that can define the dihedral
                ixs = np.all(atomnames == dih_atoms, axis=1)
                indices = np.logical_or(indices, ixs) # Combines the Trues of both arrays into one (union)
            # And then save the subset of each of the MI matrices for each dihedral
            for mi, name in MIs:
                # Define the class name with the dihedral and the MI matrix name, construct an object of the class and use the _save_pq method
                classname = f"CARDS_{name}_{dih.capitalize()}"
                exec(f"class {classname}(CARDS_Base): pass")
                dihclass = eval(classname)
                args = mi[np.ix_(indices, indices)], xtc, resids[indices]
                dihclass(self.protein, self._d)._save_pq(args)
        
        # Save an empty pq file so that the wait_calculate launched in the class' __init__ exits
        pandas.DataFrame(columns=[f"{num}" for num in self._trajs]).to_parquet(self._rawpq(xtc))

        
        
        
        
for dih in ["phi", "psi", "chi1", "chi2", "chi3", "chi4"]:
    for name in ("MI",  "Disorder", "Disorder_mediated", "Holistic"):
        exec(f'''
class CARDS_{name}_{dih.capitalize()}(Base): 
    """CARDS network construction method's {name + " MI" if name != "MI" else name} of {dih.capitalize()} dihedral
    """
    pass
''')
        

        
        
for name in ("MI",  "Disorder", "Disorder_mediated", "Holistic"):
    exec(f'''
class CARDS_{name}_Backbone_Dihs(Combined_Dihs_Avg):
    """CARDS network construction method's of the combination of the backbone
    dihedrals' {name + " MIs" if name != "MI" else "MIs"} by averaging.
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dihs = ["Phi", "Psi"]
        return new
''')
        
        
        
for name in ("MI",  "Disorder", "Disorder_mediated", "Holistic"):
    exec(f'''
class CARDS_{name}_Sidechain_Dihs(Combined_Dihs_Avg):
    """CARDS network construction method's of the combination of the side-chain
    dihedrals' {name + " MIs" if name != "MI" else "MIs"} by averaging.
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = [f"Chi{{i+1}}" for i in range(new.chis)]
        return new
''')
        
        
        
for name in ("MI",  "Disorder", "Disorder_mediated", "Holistic"):
    exec(f'''
class CARDS_{name}_Dihs(Combined_Dihs_Avg):
    """CARDS network construction method's of the combination of all backbone and
    side-chain dihedrals' {name + " MIs" if name != "MI" else "MIs"} by averaging.
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new.chis = d["chis"] if "chis" in d else 4
        new._dihs = ["Phi", "Psi"] + [f"Chi{{i+1}}" for i in range(new.chis)]
        return new
''')