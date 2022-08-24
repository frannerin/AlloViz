"""MD-TASK wrapper

It calculates the Pearson's correlation of the residues' CA atoms.

"""

from .Base import lazy_import, Multicore

imports = {
"_mdtask": ".Packages.MD-TASK.mdtask.calc_correlation",
}

for key, val in imports.items():
    exec(f"{key} = lazy_import(*{key, val})")
    
    
    

class MDTASK(Multicore):
    """MD-TASK's Pearson's correlation of CA atoms
    """
    def __new__(cls, protein, d):
        new = super().__new__(cls, protein, d)        
        new._empties = 2
        return new
    
    
    def _computation(self, xtc):#pdb, traj, xtc, pq):
        """"""
        corr = _mdtask.correlate(_mdtask.parse_traj(traj = self._trajs[xtc], topology = self._pdbf))
        return corr, xtc