from Pkgs import *
import pytraj




class PytrajCA(Correlationpkg):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)
    
    
    
    def _computation(self, pdb, traj, xtc, pq):
        top = pytraj.load_topology(pdb)
        traj = pytraj.load(traj, top, mask = f'@{mask}')
        corr = pytraj.matrix.correl(traj, f'@{mask}')
        return corr, xtc, pq

    
    
class PytrajCB(PytrajCA):
    def __init__(self, state):
        super().__init__(state)
    
    
    def _calculate(self, xtc):
        pool, pdb, traj, pq = super()._calculate(xtc)
        mask = self.__class__.__name__[-2:]
        self.selection = "protein and not resname GLY"
        
        pool.apply_async(self._computation,
                         args=(pdb, traj, mask, xtc, pq),
                         callback=self._save_pq)
