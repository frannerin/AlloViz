from .Base import Base, Use_COM, Combined_Dihs

from ..AlloViz.utils import lazy_import

imports = {
"_corrplus": "..Packages.correlationplus.correlationplus.calculate",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
class correlationplus_CA_Pear(Base):    
    
    def _computation(self, xtc):
        corr = _corrplus.calcMDnDCC(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc

        
        
        
class correlationplus_CA_LMI(correlationplus_CA_Pear):
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _corrplus.calcMD_LMI(self._pdbf, self._trajs[xtc], saveMatrix = False)
        return corr, xtc
        
        
        
class correlationplus_COM_Pear(Use_COM, correlationplus_CA_Pear):
    pass
        

        
        
class correlationplus_COM_LMI(Use_COM, correlationplus_CA_LMI):
    pass

        
        
        
class correlationplus_Psi(Base):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "psi"
        return new
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        corr = _corrplus.calcMDsingleDihedralCC(self._pdbf, self._trajs[xtc], dihedralType = self._dih, saveMatrix = False) # outputs a n_res x n_res matrix nevertheless
        return corr, xtc, self._d["_dihedral_residx"]()#[1, -1]
    
    
    
class correlationplus_Phi(correlationplus_Psi):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "phi"
        return new

        

class correlationplus_Omega(correlationplus_Psi):
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._dih = "omega"
        return new

        
        
class correlationplus_Dihs(Combined_Dihs, correlationplus_CA_Pear):
    pass