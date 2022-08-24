"""PyInteraph2 wrapper

It calculates interaction energies and contact frequencies.

"""

from .Base import lazy_import, Base

imports = {
"_pyinteraph": "pyinteraph.main",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    
    

class PyInteraph2_Base(Base):
    """PyInteraph2 base class
    """    
    def _computation(self, xtc):#pdb, traj, xtc, pq, CLIargs):
        """"""
        corr = _pyinteraph.main(f"-s {self._pdbf} -t {self._trajs[xtc]} {self._CLIargs}".split())        
        return corr, xtc
    
    
    
class PyInteraph2_Contacts(PyInteraph2_Base):
    """PyInteraph2's contact frequencies
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        
        #reslist = set(new._d["u"].select_atoms(d["_protein_sel"]).residues.resnames)
        new._CLIargs = f"-m --cmpsn-graph dummy"# --cmpsn-residues {','.join(reslist)}"
        
        #new._bonded_cys_indices = new.protein._bonded_cys#(new._d["_pdbf"])
        
        return new
                
        
    def _computation(self, xtc):#pdb, traj, xtc, pq, CLIargs):
        """"""
        corr, xtc = super()._computation(xtc)#corr = _pyinteraph.main(f"-s {self._pdbf} -t {self._trajs[xtc]} {self._CLIargs}".split()) / 100
        corr = corr / 100
        
        bonded_cys_indices = self._d["_bonded_cys"]
        
        if len(bonded_cys_indices) > 1:
            for res1 in bonded_cys_indices:
                for res2 in bonded_cys_indices:
                    if res1 != res2:
                        corr[res1, res2] = 0
        
        return corr, xtc

    
    
class PyInteraph2_Energy(PyInteraph2_Base):
    """PyInteraph2's interaction energies
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)
        new._CLIargs = "-p --kbp-graph dummy"
        return new