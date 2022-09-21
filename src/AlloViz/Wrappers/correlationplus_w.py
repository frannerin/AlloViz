"""correlationplus wrapper

It calculates the Pearson's correlation and the Linear Mutual Information (MI) of the residues'
CA atoms and COMs, and also the Pearson's correlation of the residues' backbone dihedral angles 
(Phi, Psi and Omega) and their average.

"""

from .Base import lazy_import, Base, Use_COM, Combined_Dihs_Avg, Combined_Dihs_Max

imports = {
"_corrplus": ".Packages.correlationplus.correlationplus.calculate",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
class correlationplus_CA_Pear(Base):
    """correlationplus' Pearson's correlation of CA atoms
    """
    def _computation(self, xtc):
        corr = _corrplus.calcMDnDCC(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc

        
        
        
class correlationplus_CA_LMI(correlationplus_CA_Pear):
    """correlationplus' LMI of CA atoms
    """    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _corrplus.calcMD_LMI(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc
        
        
        
class correlationplus_COM_Pear(Use_COM, correlationplus_CA_Pear):
    """correlationplus' Pearson's correlation of residues' COM
    """
    pass
        

        
        
class correlationplus_COM_LMI(Use_COM, correlationplus_CA_LMI):
    """correlationplus' LMI of residues' COM
    """    
    pass

        
        
        
class correlationplus_Psi(Base):
    """correlationplus' Pearson's correlation of Psi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "psi"
        return new
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        """"""
        corr = _corrplus.calcMDsingleDihedralCC(self._pdbf, self._trajs[xtc], dihedralType = self._dih, saveMatrix = False) # outputs a n_res x n_res matrix nevertheless
        return corr, xtc, self._d["_dihedral_residx"]()
    
    
    
class correlationplus_Phi(correlationplus_Psi):
    """correlationplus' Pearson's correlation of Phi backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "phi"
        return new

        

class correlationplus_Omega(correlationplus_Psi):
    """correlationplus' Pearson's correlation of Omega backbone dihedral
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "omega"
        return new

        
        
class correlationplus_Backbone_Dihs_Avg(Combined_Dihs_Avg, correlationplus_CA_Pear):
    """correlationplus' combination of the backbone dihedrals' Pearson's correlations by
    averaging
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dihs = ["Phi", "Psi", "Omega"]
        return new
    
    
class correlationplus_Backbone_Dihs_Max(Combined_Dihs_Max, correlationplus_CA_Pear):
    """correlationplus' combination of the backbone dihedrals' Pearson's correlations by
    taking the max value
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dihs = ["Phi", "Psi", "Omega"]
        return new