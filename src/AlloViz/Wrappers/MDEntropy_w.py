"""MDEntropy wrapper

It calculates the Mutual Information (MI) of the combination of the backbone's dihedral angles, of
the Alpha Angles (torsion angle defined by the i-1, i, i+1 and i+2 residues' CA atoms), and of the
Contacts.

"""

import numpy as np

from .Base import lazy_import, Multicore

imports = {
"_mdtraj": "mdtraj",
"_mdentropy": ".Packages.mdentropy.mdentropy.metrics",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    


class MDEntropy_Base(Multicore):
    """MDEntropy base class
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        new.method = d["MDEntropy_method"] if "MDEntropy_method" in d else "grassberger"
        
        # # Opción si las necesidades de memoria incrementan con taskcpus
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3
        
        # # Opción si las necesidades altas de memoria sólo son con el inicio del cálculo
        # extra_taskcpus = int((self.taskcpus/4 - 1) * 4) if self.taskcpus>=4 else 0
        # taskcpus = 1 + extra_taskcpus # Minimum 4*2000 of memory (this taskcpus is like that to use 4 cpus-2000 mem in .sh files)
        # empties = 3 - extra_taskcpus if extra_taskcpus<=3 else 0

        # Opción de 1 empty por tasckpu, pasando las extras a empties; 1 taskcpu requiere 3,2G; hay una necesidad mayor de memoria al principio pero i pretend i do not see
        half = int(np.floor(new.taskcpus/2))
        new.taskcpus = half if new.taskcpus > 1 else 1
        new._empties = half

        return new
        
        
        
    def _computation(self, xtc):#pdb, traj, xtc, pq, taskcpus):
        """"""
        mytraj = _mdtraj.load(self._trajs[xtc], top=self._pdbf) # hopefully mdtraj is loaded from the Classes module
        mi = self._function(threads=self.taskcpus, normed=True, method=self.method, **self._types) # n_bins=3, method='knn', normed=True
        corr = mi.partial_transform(traj=mytraj, shuffle=0, verbose=True)
        return corr, xtc, self._resl
    
    
    

    
class MDEntropy_Contacts(MDEntropy_Base):
    """MDEntropy's MI of Contacts
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._function = _mdentropy.ContactMutualInformation
        new._types = {}
        new._resl = []
        return new

    
        
        
        
        
class MDEntropy_Backbone_Dihs(MDEntropy_Base):
    """MDEntropy's MI of the backbone's dihedrals
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        new._function = _mdentropy.DihedralMutualInformation
        new._types = {"types": ["phi", "psi"]}#, "omega"]}
        new._resl = new._d["_dihedral_residx"]()
        return new
    
    
class MDEntropy_Phi(MDEntropy_Backbone_Dihs):
    """MDEntropy's MI of the phi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        new._types = {"types": ["phi"]}
        return new
    
    
class MDEntropy_Psi(MDEntropy_Backbone_Dihs):
    """MDEntropy's MI of the psi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        new._types = {"types": ["psi"]}
        return new
    
    
    
    
    

class MDEntropy_AlphaAngle(MDEntropy_Base):
    """MDEntropy's MI of the Alpha Angles
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        new._function = _mdentropy.AlphaAngleMutualInformation
        new._types = {"types": ["alpha"]}
        new._resl = new._d["_dihedral_residx"](-2)
        return new
    
    
    def _computation(self, xtc):
        """"""
        corr, xtc, dihedral_resl = super()._computation(xtc)
        
        length = corr.shape[0]
        to_insert = {0, length+1, length+2}
        new_length = length + len(to_insert)
        new_indices = set(range(new_length))
        
        new_corr = np.zeros((new_length, new_length), dtype=corr.dtype)
        corr_indices = np.array(list(new_indices - to_insert))
        new_corr[corr_indices.reshape(-1,1), corr_indices] = corr
        
        return new_corr, xtc, dihedral_resl